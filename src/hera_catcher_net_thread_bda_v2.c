/* hera_catcher_net_thread.c
 *
 * Routine to read packets sent to the catcher from 
 * network and put them into shared memory blocks. 
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

#include <xgpu.h>

#include "hashpipe.h"
#include "paper_databuf.h"


#define DEBUG_NET

#ifndef MIN
#define MIN(a,b) ((a) < (b) ? (a) : (b))
#endif

#ifndef MAX
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#endif

#define PKTSOCK_BYTES_PER_FRAME       (4864)
#define PKTSOCK_FRAMES_PER_BLOCK      (128)
#define PKTSOCK_NBLOCKS               (5000)
#define PKTSOCK_NFRAMES               (PKTSOCK_FRAMES_PER_BLOCK * PKTSOCK_NBLOCKS)

#define ELAPSED_NS(start,stop) \
  (((int64_t)stop.tv_sec-start.tv_sec)*1000*1000*1000+(stop.tv_nsec-start.tv_nsec))

#define MAX_OUT_OF_SEQ_PKTS           (4096)

typedef struct{
    uint64_t mcnt;        // timestamp of the packet
    uint32_t bcnt;        // unique id based on order of baselines sent to catcher
    uint32_t offset;      // channel within one x-eng 
    uint16_t ant0;
    uint16_t ant1;
    uint16_t xeng_id;	  // For time demux and starting channel 
    uint16_t payload_len; // should be 4096
} packet_header_t;

// The fields of a block_info_t structure hold meta-data 
// about the contents of each block
typedef struct {
    int initialized;
    uint32_t bcnt_start; 
    int block_curr;
    int block_next;
    long out_of_seq_cnt;
    long block_packet_counter[CATCHER_N_BLOCKS];
    char flags[CATCHER_N_BLOCKS][PACKETS_PER_BLOCK];
    char baselines[CATCHER_N_BLOCKS][BASELINES_PER_BLOCK];
} block_info_t;

static hashpipe_status_t *st_p;
static const char * status_key;

// Get physical block number given the bcnt
static inline int block_for_bcnt(uint32_t bcnt){
    return (bcnt/BASELINES_PER_BLOCK) % CATCHER_N_BLOCKS;
}

// Initialize a block by clearing its "good data" flag and saving the 
// bcnt of the first baseline in this block. The bcnt should be a multiple 
// of BASELINES_PER_BLOCK.
static inline void initialize_block(hera_catcher_input_databuf_t * db, uint64_t bcnt){
  int block_i = block_for_bcnt(bcnt);
  db->block[block_i].header.bcnt[0]   = bcnt;
  db->block[block_i].header.good_data = 0;
}

/* Initialize block_info */
// This function must be called once and only once per block_info structure!
// Subsequent calls are no-ops.
static inline void initialize_block_info(block_info_t * binfo, uint32_t bcnt){
    int i;

    // If this block_info structure has already been initialized
    if(binfo->initialized) {
	  return;
    }

    for(i = 0; i < CATCHER_N_BLOCKS; i++) {
      binfo->block_packet_counter[i] = 0;
      memset(binfo->flags[i], 0, PACKETS_PER_BLOCK*sizeof(char));
      memset(binfo->baselines[i], 0, BASELINES_PER_BLOCK*sizeof(char));
    }

    // The number of the current block is set by the first packet received
    binfo->bcnt_start     = bcnt;
    binfo->block_curr     = block_for_bcnt(bcnt);
    binfo->block_next     = (binfo->block_curr+1) % CATCHER_N_BLOCKS;
    binfo->out_of_seq_cnt = 0;
    binfo->initialized    = 1;
}

/* Get packet header */
static inline void get_header(unsigned char *p_frame, packet_header_t *pkt_header){
   packet_header_t *packet_header_raw = (packet_header_t *)PKT_UDP_DATA(p_frame);
   pkt_header->mcnt        = be64toh(packet_header_raw->mcnt);
   pkt_header->bcnt        = be32toh(packet_header_raw->bcnt);
   pkt_header->offset      = be32toh(packet_header_raw->offset);
   pkt_header->ant0        = be16toh(packet_header_raw->ant0);
   pkt_header->ant1        = be16toh(packet_header_raw->ant1);
   pkt_header->xeng_id     = be16toh(packet_header_raw->xeng_id);
   pkt_header->payload_len = be16toh(packet_header_raw->payload_len);
}

/* Set hashpipe block to filled */
// This sets the "current" block to be marked as filled.
// Returns bcnt of the block being marked filled.
static uint64_t set_block_filled(hera_catcher_input_databuf_t *db, block_info_t *binfo){
  static int last_filled = -1; 
  if (last_filled==-1){ 
    last_filled = (binfo->block_curr-1)% CATCHER_N_BLOCKS;
  }

  uint64_t block_missed_pkt_cnt;
  uint32_t missed_pkt_cnt;
  uint32_t block_i = binfo->block_curr;

  // Validate that we're filling blocks in the proper sequence
  last_filled = (last_filled+1) % CATCHER_N_BLOCKS;
  if(last_filled != block_i) {
    printf("block %d being marked filled, but expected block %d!\n", block_i, last_filled);
  }

  db->block[block_i].header.good_data = 1;

  // Set the block as filled
  if(hera_catcher_input_databuf_set_filled(db, block_i) != HASHPIPE_OK){
    hashpipe_error(__FUNCTION__, "error waiting for databuf filled call");
    pthread_exit(NULL);
  }

  // Calculate missing packets.
  block_missed_pkt_cnt = PACKETS_PER_BLOCK - binfo->block_packet_counter[block_i];

  // Update status buffer
  hashpipe_status_lock_busywait_safe(st_p);
  hputu4(st_p->buf, "NETBKOUT", block_i);
  if(block_missed_pkt_cnt){
    fprintf(stderr, "Expected %lu packets, Got %lu\n", PACKETS_PER_BLOCK, binfo->block_packet_counter[block_i]);
    // Increment MISSEDPK by number of missed packets for this block
    hgetu4(st_p->buf, "MISSEDPK", &missed_pkt_cnt);
    missed_pkt_cnt += block_missed_pkt_cnt;
    hputu4(st_p->buf, "MISSEDPK", missed_pkt_cnt);
    //  fprintf(stderr, "got %d packets instead of %d\n",
    //	    binfo->block_packet_counter[block_i], N_PACKETS_PER_BLOCK);
  }
  hashpipe_status_unlock_safe(st_p);

  return db->block[block_i].header.bcnt[0];
}


/* Process packet */
// This function returns -1 unless the given packet causes a block to be marked
// as filled in which case this function returns the marked block's first bcnt.
// Any return value other than -1 will be stored in the status memory as
// NETMCNT, so it is important that values other than -1 are returned rarely
// (i.e. when marking a block as filled)!!!
static inline uint64_t process_packet(
  hera_catcher_input_databuf_t *db, unsigned char *p_frame){

  static block_info_t binfo;
  packet_header_t pkt_header;
  const uint32_t *payload_p;
  uint32_t *dest_p;
  int i;
  int pkt_block_i;
  int b, x, t, o;
  uint32_t pkt_offset;
  uint32_t pkt_bcnt_dist;
  uint64_t bcnt = -1;

  // Parse packet header
  get_header(p_frame, &pkt_header);
  pkt_block_i = block_for_bcnt(pkt_header.bcnt);

  // Lazy init binfo
  if(!binfo.initialized){
    fprintf(stderr,"Initializing binfo..!\n");
    initialize_block_info(&binfo, pkt_header.bcnt);

    fprintf(stderr,"Initializing the first blocks..\n");
    // Initialize the newly acquired block
    initialize_block(db, pkt_header.bcnt); 
    initialize_block(db, pkt_header.bcnt+BASELINES_PER_BLOCK); 
    initialize_block(db, pkt_header.bcnt+2*BASELINES_PER_BLOCK); 
  }

  //fprintf(stderr, "curr:%d\tnext:%d\t",binfo.block_curr, binfo.block_next);
  //fprintf(stderr, "bcnt:%d\tblock_id:%d\n",pkt_header.bcnt, pkt_block_i);
  //fprintf(stderr, "xeng:%d\n",pkt_header.xeng_id);

  // Packet bcnt distance (how far away is this packet's bcnt from the
  // current bcnt).  Positive distance for pcnt mcnts > current mcnt.
  pkt_bcnt_dist = pkt_header.bcnt - binfo.bcnt_start; 

  // We expect packets for the current block, the next block and the block after.
  if (0 <= pkt_bcnt_dist && pkt_bcnt_dist < 3*BASELINES_PER_BLOCK){
    // If the packet is for the block after the next block (i.e. current 
    // block + 2 blocks), mark the current block as filled.
    if (pkt_bcnt_dist >= 2*BASELINES_PER_BLOCK){
       
       bcnt = set_block_filled(db, &binfo);

       // Print
       fprintf(stderr,"Filled Block: %d with bcnt: %lu\n",binfo.block_curr, bcnt);

       // Update binfo
       binfo.bcnt_start += BASELINES_PER_BLOCK;
       binfo.block_packet_counter[binfo.block_curr] = 0;
       memset(binfo.flags[binfo.block_curr],     0, PACKETS_PER_BLOCK*sizeof(char));
       memset(binfo.baselines[binfo.block_curr], 0, BASELINES_PER_BLOCK*sizeof(char));

       binfo.block_curr = binfo.block_next;
       binfo.block_next = (binfo.block_next+1) % CATCHER_N_BLOCKS; 

       // Wait (hopefully not long!) to acquire the block after next.
       if(hera_catcher_input_databuf_busywait_free(db, pkt_block_i) != HASHPIPE_OK) {
       if (errno == EINTR) {
           // Interrupted by signal, return -1
           hashpipe_error(__FUNCTION__, "interrupted by signal waiting for free databuf");
           pthread_exit(NULL);
           return -1; // We're exiting so return value is kind of moot
       } else {
           hashpipe_error(__FUNCTION__, "error waiting for free databuf");
           pthread_exit(NULL);
           return -1; // We're exiting so return value is kind of moot
       }
       }

       // Initialize the newly acquired block
       initialize_block(db, bcnt+2*BASELINES_PER_BLOCK); 

       // Reset out-of-seq counter
       binfo.out_of_seq_cnt = 0;
    }

    // Evaluate the location in the buffer to which to copy the packet data   
    b = pkt_header.bcnt % BASELINES_PER_BLOCK;
    x = pkt_header.xeng_id % N_XENGINES_PER_TIME;
    t = (pkt_header.mcnt/Nt) % TIME_DEMUX;  //Nt = 2
    o = pkt_header.offset;
    pkt_offset = hera_catcher_input_databuf_pkt_offset(b, t, x, o);
    //fprintf(stderr, "offset: %d\n", pkt_offset);
    //fprintf(stderr, "bcnt-loc:%d\txeng:%d\ttime:%d\tpktoffset:%d\n",b,x,t,o);
    
    // Copy data into buffer with byte swap
    payload_p = (uint32_t *)(PKT_UDP_DATA(p_frame) + sizeof(packet_header_t)); 
    dest_p    = (uint32_t *)(db->block[pkt_block_i].data + (pkt_header.payload_len * pkt_offset/sizeof(uint32_t)));

    for(i=0; i<(pkt_header.payload_len>>2); i++){
      dest_p[i] = be32toh(payload_p[i]);
    }
    //fprintf(stderr,"bcnt:%d\t block_id:%d\t Pkt cntr: %lu\n",b, pkt_block_i, binfo.block_packet_counter[pkt_block_i]);

    // If this is the first packet of this baseline, update header
    if(!binfo.baselines[pkt_block_i][b]){
      db->block[pkt_block_i].header.mcnt[b] = pkt_header.mcnt;
      db->block[pkt_block_i].header.ant_pair_0[b] = pkt_header.ant0;
      db->block[pkt_block_i].header.ant_pair_1[b] = pkt_header.ant1;
      db->block[pkt_block_i].header.bcnt[b] = pkt_header.bcnt;
      binfo.baselines[pkt_block_i][b] = 1;
    }
    // Update binfo
    binfo.flags[pkt_block_i][pkt_offset] = 1;
    binfo.block_packet_counter[pkt_block_i]++;

    // Check for duplicate packets
    // if(binfo.flags[pkt_block_i][pkt_offset]){
    //   // This slot is already filled
    //   //fprintf(stderr, "Packet repeated!!\n");
    //   binfo.out_of_seq_cnt++;
    //   return -1;
    // }

  }else{
    binfo.out_of_seq_cnt++;
    return -1;
  }
  
  //if (binfo.out_of_seq_cnt > MAX_OUT_OF_SEQ_PKTS){
  //  fprintf(stderr,"Lots of out of sequence packets!");
  //}
  return bcnt;  
}


static int init(hashpipe_thread_args_t *args){
  /* Read network params */
  char bindhost[80];
  int bindport = CATCHER_PORT;
  status_key = args->thread_desc->skey;

  strcpy(bindhost, "0.0.0.0");

  hashpipe_status_t st = args->st;

  hashpipe_status_lock_safe(&st);
  // Record version
  hputs(st.buf, "GIT_VER", GIT_VERSION);
  // Get info from status buffer if present (no change if not present)
  hgets(st.buf, "BINDHOST", 80, bindhost);
  hgeti4(st.buf, "BINDPORT", &bindport);
  // Store bind host/port info etc in status buffer
  hputs(st.buf, "BINDHOST", bindhost);
  hputi4(st.buf, "BINDPORT", bindport);
  hputu4(st.buf, "MISSEDFE", 0);
  hputu4(st.buf, "MISSEDPK", 0);
  hashpipe_status_unlock_safe(&st);

#ifndef TIMING_TEST
    /* Set up pktsock */
    struct hashpipe_pktsock *p_ps = (struct hashpipe_pktsock *)
    malloc(sizeof(struct hashpipe_pktsock));

    if(!p_ps) {
        perror(__FUNCTION__);
        return -1;
    }
  
    // Make frame_size be a divisor of block size so that frames will be
    // contiguous in mapped mempory.  block_size must also be a multiple of
    // page_size.  Easiest way is to oversize the frames to be 16384 bytes, which
    // is bigger than we need, but keeps things easy.
    p_ps->frame_size = PKTSOCK_BYTES_PER_FRAME;
    // total number of frames
    p_ps->nframes = PKTSOCK_NFRAMES;
    // number of blocks
    p_ps->nblocks = PKTSOCK_NBLOCKS;
  
    int rv = hashpipe_pktsock_open(p_ps, bindhost, PACKET_RX_RING);
    if (rv!=HASHPIPE_OK) {
        hashpipe_error("hera_pktsock_thread", "Error opening pktsock.");
        pthread_exit(NULL);
    }
  
    // Store packet socket pointer in args
    args->user_data = p_ps;
#endif

  // Success!
  return 0;
}

static void *run(hashpipe_thread_args_t * args){
  // Local aliases to shorten access to args fields
  // Our output buffer happens to be a paper_input_databuf
  hera_catcher_input_databuf_t *db = (hera_catcher_input_databuf_t *)args->obuf;
  hashpipe_status_t st = args->st;
  const char *status_key = args->thread_desc->skey;

  st_p = &st;	// allow global (this source file) access to the status buffer

  // Flag that holds off the net thread
  int holdoff = 0;

  // Force ourself into the hold off state
  //fprintf(stderr, "Setting CNETHOLD state to 1.Waiting for someone to set it to 0\n");
  //hashpipe_status_lock_safe(&st);
  //hputi4(st.buf, "CNETHOLD", 1);
  //hputs(st.buf, status_key, "holding");
  //hashpipe_status_unlock_safe(&st);

  while(holdoff) {
    sleep(1);
    hashpipe_status_lock_safe(&st);
    hgeti4(st.buf, "CNETHOLD", &holdoff);
    if(!holdoff) {
      // Done holding, so delete the key
      hdel(st.buf, "CNETHOLD");
      hputs(st.buf, status_key, "starting");
      fprintf(stderr, "Starting...\n");
    }
    hashpipe_status_unlock_safe(&st);
  }

  fprintf(stderr,"Waiting to acquire three blocks to start!\n");

  // Acquire first two blocks to start
  if(hera_catcher_input_databuf_busywait_free(db, 0) != HASHPIPE_OK){
    if (errno == EINTR){
      // Interrupted by signal, return -1
      hashpipe_error(__FUNCTION__, "interrupted by signal waiting for free databuf");
      pthread_exit(NULL);
    }else{
        hashpipe_error(__FUNCTION__, "error waiting for free databuf");
        pthread_exit(NULL);
    }
  }if(hera_catcher_input_databuf_busywait_free(db, 1) != HASHPIPE_OK){
    if (errno == EINTR){
      // Interrupted by signal, return -1
      hashpipe_error(__FUNCTION__, "interrupted by signal waiting for free databuf");
      pthread_exit(NULL);
    }else{
        hashpipe_error(__FUNCTION__, "error waiting for free databuf");
        pthread_exit(NULL);
    }
  }if(hera_catcher_input_databuf_busywait_free(db, 2) != HASHPIPE_OK){
    if (errno == EINTR){
      // Interrupted by signal, return -1
      hashpipe_error(__FUNCTION__, "interrupted by signal waiting for free databuf");
      pthread_exit(NULL);
    }else{
        hashpipe_error(__FUNCTION__, "error waiting for free databuf");
        pthread_exit(NULL);
    }
  }


  /* Read network params */
  int bindport = CATCHER_PORT;
  // (N_BYTES_PER_PACKET excludes header)
  size_t expected_packet_size = OUTPUT_BYTES_PER_PACKET + sizeof(packet_header_t);

#ifndef TIMING_TEST
    /* Get pktsock from args*/
    struct hashpipe_pktsock * p_ps = (struct hashpipe_pktsock*)args->user_data;
    pthread_cleanup_push(free, p_ps);
    pthread_cleanup_push((void (*)(void *))hashpipe_pktsock_close, p_ps);

    // Drop all packets to date
    unsigned char *p_frame;
    while((p_frame=hashpipe_pktsock_recv_frame_nonblock(p_ps))) {
	hashpipe_pktsock_release_frame(p_frame);
    }

    hashpipe_status_lock_safe(&st);
    // Get info from status buffer
    hgeti4(st.buf, "BINDPORT", &bindport);
    hputu4(st.buf, "MISSEDPK", 0);
    hputs(st.buf, status_key, "running");
    hashpipe_status_unlock_safe(&st);
#endif


  /* Main loop */
  uint64_t packet_count = 0;
  uint64_t wait_ns = 0; // ns for most recent wait
  uint64_t recv_ns = 0; // ns for most recent recv
  uint64_t proc_ns = 0; // ns for most recent proc
  uint64_t min_wait_ns = 99999; // min ns per single wait
  uint64_t min_recv_ns = 99999; // min ns per single recv
  uint64_t min_proc_ns = 99999; // min ns per single proc
  uint64_t max_wait_ns = 0;     // max ns per single wait
  uint64_t max_recv_ns = 0;     // max ns per single recv
  uint64_t max_proc_ns = 0;     // max ns per single proc
  uint64_t elapsed_wait_ns = 0; // cumulative wait time per block
  uint64_t elapsed_recv_ns = 0; // cumulative recv time per block
  uint64_t elapsed_proc_ns = 0; // cumulative proc time per block
  uint64_t status_ns = 0; // User to fetch ns values from status buffer
  float ns_per_wait = 0.0; // Average ns per wait over 1 block
  float ns_per_recv = 0.0; // Average ns per recv over 1 block
  float ns_per_proc = 0.0; // Average ns per proc over 1 block
  unsigned int pktsock_pkts = 0;  // Stats counter from socket packet
  unsigned int pktsock_drops = 0; // Stats counter from socket packet
  uint64_t pktsock_pkts_total = 0;  // Stats total for socket packet
  uint64_t pktsock_drops_total = 0; // Stats total for socket packet
  struct timespec start, stop;
  struct timespec recv_start, recv_stop;


  fprintf(stderr,"Starting to collect packets..\n");
  fprintf(stderr,"Channels per catcher packet:%ld\n",CHAN_PER_CATCHER_PKT);
  fprintf(stderr,"Packets per block: %ld\n",PACKETS_PER_BLOCK);

  while (run_threads()) {
    // Read packet
    clock_gettime(CLOCK_MONOTONIC, &recv_start);
    do{
      clock_gettime(CLOCK_MONOTONIC, &start);
      p_frame = hashpipe_pktsock_recv_udp_frame_nonblock(p_ps, bindport);
      clock_gettime(CLOCK_MONOTONIC, &recv_stop);
    } while (!p_frame && run_threads());

    // Make sure received packet size matches expected packet size.
    // -8 for the UDP header (not the *App* header!)
    int packet_size = PKT_UDP_SIZE(p_frame) - 8; 
    if (expected_packet_size != packet_size) {
      // Log warning and ignore wrongly sized packet
      #ifdef DEBUG_NET
        hashpipe_warn("hera_pktsock_thread", "Invalid pkt size (%d)", packet_size);
      #endif
      hashpipe_pktsock_release_frame(p_frame);
      continue;
    }

    packet_count++;
    
    // Copy packet into any blocks where it belongs.
    const uint64_t bcnt = process_packet((hera_catcher_input_databuf_t *)db, p_frame);
    // Release frame back to kernel
    hashpipe_pktsock_release_frame(p_frame);

    clock_gettime(CLOCK_MONOTONIC, &stop);
    wait_ns = ELAPSED_NS(recv_start, start);
    recv_ns = ELAPSED_NS(start, recv_stop);
    proc_ns = ELAPSED_NS(recv_stop, stop);
    elapsed_wait_ns += wait_ns;
    elapsed_recv_ns += recv_ns;
    elapsed_proc_ns += proc_ns;
    // Update min max values
    min_wait_ns = MIN(wait_ns, min_wait_ns);
    min_recv_ns = MIN(recv_ns, min_recv_ns);
    min_proc_ns = MIN(proc_ns, min_proc_ns);
    max_wait_ns = MAX(wait_ns, max_wait_ns);
    max_recv_ns = MAX(recv_ns, max_recv_ns);
    max_proc_ns = MAX(proc_ns, max_proc_ns);

    if(bcnt != -1) {
      // Update status
      ns_per_wait = (float)elapsed_wait_ns / packet_count;
      ns_per_recv = (float)elapsed_recv_ns / packet_count;
      ns_per_proc = (float)elapsed_proc_ns / packet_count;

      // Get stats from packet socket
      hashpipe_pktsock_stats(p_ps, &pktsock_pkts, &pktsock_drops);

      hashpipe_status_lock_busywait_safe(&st);

      hputu8(st.buf, "NETBCNT", bcnt);
      // Gbps = bits_per_packet / ns_per_packet
      // (N_BYTES_PER_PACKET excludes header, so +8 for the header)
      hputr4(st.buf, "NETGBPS", 8*(N_BYTES_PER_PACKET+8)/(ns_per_recv+ns_per_proc));
      hputr4(st.buf, "NETWATNS", ns_per_wait);
      hputr4(st.buf, "NETRECNS", ns_per_recv);
      hputr4(st.buf, "NETPRCNS", ns_per_proc);

      // Get and put min and max values.  The "get-then-put" allows the
      // user to reset the min max values in the status buffer.
      hgeti8(st.buf, "NETWATMN", (long int *)&status_ns);
      status_ns = MIN(min_wait_ns, status_ns);
          hputi8(st.buf, "NETWATMN", status_ns);

          hgeti8(st.buf, "NETRECMN", (long int *)&status_ns);
      status_ns = MIN(min_recv_ns, status_ns);
          hputi8(st.buf, "NETRECMN", status_ns);

          hgeti8(st.buf, "NETPRCMN", (long int *)&status_ns);
      status_ns = MIN(min_proc_ns, status_ns);
          hputi8(st.buf, "NETPRCMN", status_ns);

          hgeti8(st.buf, "NETWATMX", (long int *)&status_ns);
      status_ns = MAX(max_wait_ns, status_ns);
          hputi8(st.buf, "NETWATMX", status_ns);

          hgeti8(st.buf, "NETRECMX", (long int *)&status_ns);
      status_ns = MAX(max_recv_ns, status_ns);
          hputi8(st.buf, "NETRECMX", status_ns);

          hgeti8(st.buf, "NETPRCMX", (long int *)&status_ns);
      status_ns = MAX(max_proc_ns, status_ns);
          hputi8(st.buf, "NETPRCMX", status_ns);

          hputu8(st.buf, "NETPKTS",  pktsock_pkts);
          hputu8(st.buf, "NETDROPS", pktsock_drops);

          hgetu8(st.buf, "NETPKTTL", (long unsigned int*)&pktsock_pkts_total);
          hgetu8(st.buf, "NETDRPTL", (long unsigned int*)&pktsock_drops_total);
          hputu8(st.buf, "NETPKTTL", pktsock_pkts_total + pktsock_pkts);
          hputu8(st.buf, "NETDRPTL", pktsock_drops_total + pktsock_drops);

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

#ifndef TIMING_TEST
  /* Have to close all push's */
  pthread_cleanup_pop(1); /* Closes push(hashpipe_pktsock_close) */
  pthread_cleanup_pop(1); /* Closes push(hashpipe_udp_close) */
#endif
  return NULL;
}

static hashpipe_thread_desc_t hera_catcher_net_thread = {
    name: "hera_catcher_net_thread",
    skey: "CNETSTAT",
    init: init,
    run:  run,
    ibuf_desc: {NULL},
    obuf_desc: {hera_catcher_input_databuf_create}
};

static __attribute__((constructor)) void ctor()
{
  register_hashpipe_thread(&hera_catcher_net_thread);
}

// vi: set ts=8 sw=4 noet :
