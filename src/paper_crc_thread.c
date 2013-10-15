/* paper_crc_thread.c
 *
 * Routine to read packets from network and validate their CRCs (if present).
 * This does not process the data nor pass it on to any other thread.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <stdint.h>
#include <string.h>
#include <pthread.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <sys/types.h>
#include <errno.h>
#include <zlib.h>

#include <xgpu.h>

#include "hashpipe.h"
#include "paper_databuf.h"

typedef struct {
    uint64_t mcnt;
    int      fid;	// Fengine ID
    int      xid;	// Xengine ID
} packet_header_t;

// The fields of a block_info_t structure hold (at least) two different kinds
// of data.  Some fields hold data that persist over many packets while other
// fields hold data that are only applicable to the current packet (or the
// previous packet).
typedef struct {
    int initialized;
    int32_t  self_xid;
    uint64_t mcnt_start;
    uint64_t mcnt_offset;
    uint64_t mcnt_prior;
    int out_of_seq_cnt;
    int block_i;
    // The m,x,f fields hold three of the five dimensional indices for
    // the first data word of the current packet (i.e. t=0 and c=0).
    int m; // formerly known as sub_block_i
    int f;
    int block_active[N_INPUT_BLOCKS];
} block_info_t;

static hashpipe_status_t *st_p;

static inline void get_header (struct hashpipe_udp_packet *p, packet_header_t * pkt_header)
{
    uint64_t raw_header;
    raw_header = be64toh(*(unsigned long long *)p->data);
    pkt_header->mcnt        = raw_header >> 16;
    pkt_header->xid         = raw_header        & 0x00000000000000FF;
    pkt_header->fid         = (raw_header >> 8) & 0x00000000000000FF;
}

#define ELAPSED_NS(start,stop) \
  (((int64_t)stop.tv_sec-start.tv_sec)*1000*1000*1000+(stop.tv_nsec-start.tv_nsec))

static void *run(hashpipe_thread_args_t * args)
{
    // Local aliases to shorten access to args fields
    // Our output buffer happens to be a paper_input_databuf
    hashpipe_status_t st = args->st;
    const char * status_key = args->thread_desc->skey;

    st_p = &st;	// allow global (this source file) access to the status buffer

    // Get inital value for crc32 function
    uint32_t init_crc = crc32(0,0,0);

    // Flag that holds off the crc thread
    int holdoff = 1;

    // Force ourself into the hold off state
    hashpipe_status_lock_safe(&st);
    hputi4(st.buf, "CRCHOLD", 1);
    hashpipe_status_unlock_safe(&st);

    while(holdoff) {
	// We're not in any hurry to startup
	sleep(1);
	hashpipe_status_lock_safe(&st);
	// Look for CRCHOLD value
	hgeti4(st.buf, "CRCHOLD", &holdoff);
	if(!holdoff) {
	    // Done holding, so delete the key
	    hdel(st.buf, "CRCHOLD");
	}
	hashpipe_status_unlock_safe(&st);
    }

    /* Read network params */
    struct hashpipe_udp_params up = {
	.bindhost = "0.0.0.0",
	.bindport = 8511,
	.packet_size = 8200
    };
    hashpipe_status_lock_safe(&st);
    // Get info from status buffer if present (no change if not present)
    hgets(st.buf, "BINDHOST", 80, up.bindhost);
    hgeti4(st.buf, "BINDPORT", &up.bindport);
    // Store bind host/port info etc in status buffer
    hputs(st.buf, "BINDHOST", up.bindhost);
    hputi4(st.buf, "BINDPORT", up.bindport);
    hputu4(st.buf, "CRCPKOK", 0);
    hputu4(st.buf, "CRCPKERR", 0);
    hputs(st.buf, status_key, "running");
    hashpipe_status_unlock_safe(&st);

    struct hashpipe_udp_packet p;

    /* Give all the threads a chance to start before opening network socket */
    sleep(1);


    /* Set up UDP socket */
    int rv = hashpipe_udp_init(&up);
    if (rv!=HASHPIPE_OK) {
        hashpipe_error("paper_crc_thread",
                "Error opening UDP socket.");
        pthread_exit(NULL);
    }
    pthread_cleanup_push((void *)hashpipe_udp_close, &up);

    /* Main loop */
    uint64_t packet_count = 0;
    uint64_t good_count = 0;
    uint64_t error_count = 0;
    uint64_t elapsed_wait_ns = 0;
    uint64_t elapsed_recv_ns = 0;
    uint64_t elapsed_proc_ns = 0;
    float ns_per_wait = 0.0;
    float ns_per_recv = 0.0;
    float ns_per_proc = 0.0;
    struct timespec start, stop;
    struct timespec recv_start, recv_stop;
    packet_header_t hdr;

    while (run_threads()) {

        /* Read packet */
	clock_gettime(CLOCK_MONOTONIC, &recv_start);
	do {
	    clock_gettime(CLOCK_MONOTONIC, &start);
	    p.packet_size = recv(up.sock, p.data, HASHPIPE_MAX_PACKET_SIZE, 0);
	    clock_gettime(CLOCK_MONOTONIC, &recv_stop);
	} while (p.packet_size == -1 && (errno == EAGAIN || errno == EWOULDBLOCK) && run_threads());

	// Break out of loop if stopping
	if(!run_threads()) break;

	// Increment packet count
	packet_count++;

	// Check CRC
        if(crc32(init_crc, (/*const?*/ uint8_t *)p.data, p.packet_size) == 0xffffffff) {
	    // CRC OK! Increment good counter
	    good_count++;
	} else {
	    // CRC error!  Increment error counter
	    error_count++;

	    // Log message
	    get_header(&p, &hdr);
	    hashpipe_warn("paper_crc", "CRC error mcnt %llu ; fid %u ; xid %u",
		    hdr.mcnt, hdr.fid, hdr.xid);
	}

	clock_gettime(CLOCK_MONOTONIC, &stop);
	elapsed_wait_ns += ELAPSED_NS(recv_start, start);
	elapsed_recv_ns += ELAPSED_NS(start, recv_stop);
	elapsed_proc_ns += ELAPSED_NS(recv_stop, stop);

        if(packet_count % 1000 == 0) {
	    // Compute stats
	    get_header(&p, &hdr);
            ns_per_wait = (float)elapsed_wait_ns / packet_count;
            ns_per_recv = (float)elapsed_recv_ns / packet_count;
            ns_per_proc = (float)elapsed_proc_ns / packet_count;

            // Update status
            hashpipe_status_lock_busywait_safe(&st);
            hputu8(st.buf, "CRCMCNT", hdr.mcnt);
	    // Gbps = bits_per_packet / ns_per_packet
	    // (N_BYTES_PER_PACKET excludes header, so +8 for the header)
            hputr4(st.buf, "CRCGBPS", 8*(N_BYTES_PER_PACKET+8)/(ns_per_recv+ns_per_proc));
            hputr4(st.buf, "CRCWATNS", ns_per_wait);
            hputr4(st.buf, "CRCRECNS", ns_per_recv);
            hputr4(st.buf, "CRCPRCNS", ns_per_proc);
	    // TODO Provide some way to recognize request to zero out the
	    // CRCERR and CRCOK fields.
	    hputu8(st.buf, "CRCPKOK",  good_count);
	    hputu8(st.buf, "CRCPKERR", error_count);
            hashpipe_status_unlock_safe(&st);

	    // Start new average
	    elapsed_wait_ns = 0;
	    elapsed_recv_ns = 0;
	    elapsed_proc_ns = 0;
	    packet_count = 0;
        }

        /* Will exit if thread has been cancelled */
        pthread_testcancel();
    }

    /* Have to close all push's */
    pthread_cleanup_pop(1); /* Closes push(hashpipe_udp_close) */

    return NULL;
}

static hashpipe_thread_desc_t crc_thread = {
    name: "paper_crc_thread",
    skey: "CRCSTAT",
    init: NULL,
    run:  run,
    ibuf_desc: {NULL},
    obuf_desc: {NULL}
};

static __attribute__((constructor)) void ctor()
{
  register_hashpipe_thread(&crc_thread);
}

// vi: set ts=8 sw=4 noet :
