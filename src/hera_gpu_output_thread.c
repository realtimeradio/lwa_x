/* Takes the output from the GPU and 
 * averages over time for the shorter 
 * baselines according to a config file.
 * Output is stored in CASPER ordered 
 * shared memory segments.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <string.h>
#include <pthread.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <sys/types.h>
#include <errno.h>
#include <endian.h>

#include <xgpu.h>

#include "hashpipe.h"
#include "paper_databuf.h"

#define LOG(x)           ((int)(log((x))/log(2))) // yields msb loc
#define CHECK_PWR2(x)    (!((x)&((x)-1)))

// Macros for generating values for the pkthdr_t fields
#define TIMESTAMP(x)      (htobe64((uint64_t)(x)))
#define BASELINE_ID(x)    (htobe32((uint32_t)(x)))
#define OFFSET(x)         (htobe32((uint32_t)(x)))
#define XENG_ID(x)        (htobe16((uint16_t)(x)))
#define PAYLOAD_LEN(x)    (htobe16((uint16_t)(x)))
#define ANTENNA(x)        (htobe16((uint16_t)(x)))

#define CONVERT(x)        (htobe32((x)))

typedef int32_t pktdata_t;
static XGPUInfo xgpu_info;

// Structure for packet header
typedef struct pkthdr {
  uint64_t timestamp;
  uint32_t baseline_id;
  uint32_t offset;
  uint16_t ant0;
  uint16_t ant1;
  uint16_t xeng_id;
  uint16_t payload_len;
} pkthdr_t;

// Structure of a packet
typedef struct struct_pkt {
  pkthdr_t hdr;
  pktdata_t data[OUTPUT_BYTES_PER_PACKET/sizeof(pktdata_t)];
} pkt_t;

// PACKET_DELAY_NS is number of nanoseconds to delay between packets.  This is
// to prevent overflowing the network interface's TX queue.
// Keep in mind the total throughput of the network, and also the number
// of x-engine instances which will be running on each host.

// Note that the below calculations account only for the packet's payload (i.e.
// it considers the size of the packet header to be negligible).
//
// 1000 megabit per second =   1 nanosecond per bit
//  100 megabit per second =  10 nanosecond per bit
//   10 megabit per second = 100 nanosecond per bit

// Set to 200 Mbps -- OK for two instances per node.
// With 16 nodes, amounts to 6.4 Gbps of data
#define PACKET_DELAY_NS (4 * XENG_CHAN_SUM * 8*OUTPUT_BYTES_PER_PACKET)

// bytes_per_dump depends on xgpu_info.triLength
static uint64_t bytes_per_dump = 0;
// packets_per_dump is bytes_per_dump / OUTPUT_BYTES_PER_PACKET
static unsigned int packets_per_dump = 0;

// Open and connect a UDP socket to the given host and port.  Note that port is
// a string and can be a stringified integer (e.g. "7148") or a service name
// (e.g. "ntp").  Returns -1 on error, otherwise a valid descriptor for an open
// and connected socket.
static int
open_udp_socket(const char *host, const char *port)
{
    struct addrinfo hints;
    struct addrinfo *result, *rp;
    int s, sfd=-1;

    // Obtain address(es) matching host/port

    memset(&hints, 0, sizeof(struct addrinfo));
    hints.ai_family = AF_INET;
    hints.ai_socktype = SOCK_DGRAM;
    hints.ai_flags = 0;
    hints.ai_protocol = 0;

    s = getaddrinfo(host, port, &hints, &result);
    if (s != 0) {
        hashpipe_error("getaddrinfo", gai_strerror(s));
        return -1;
    }

    // getaddrinfo() returns a list of address structures.
    // Try each address until we successfully connect(2).
    // If socket(2) (or connect(2)) fails, we (close the socket
    // and) try the next address.

    for (rp = result; rp != NULL; rp = rp->ai_next) {
        sfd = socket(rp->ai_family, rp->ai_socktype, rp->ai_protocol);

        if (sfd == -1) {
            continue;
        }

        if (connect(sfd, rp->ai_addr, rp->ai_addrlen) != -1) {
            break; // Success
        }

        close(sfd);
        sfd = -1;
    }

    freeaddrinfo(result);

#if 0
    // Print send buffer size
    int bufsize;
    unsigned int bufsizesize = sizeof(bufsize);
    getsockopt(sfd, SOL_SOCKET, SO_SNDBUF, &bufsize, &bufsizesize);
    printf("send buffer size is %d\n", bufsize);
#endif

    return sfd;
}


// Computes the triangular index of an (i,j) pair as shown here...
// NB: Output is valid only if i >= j.
//
//      i=0  1  2  3  4..
//     +---------------
// j=0 | 00 01 03 06 10
//   1 |    02 04 07 11
//   2 |       05 08 12
//   3 |          09 13
//   4 |             14
//   :
static inline off_t tri_index(const int i, const int j)
{
  return (i * (i+1))/2 + j;
}

// regtile_chan_length is the number of complex cross products per channel for
// the xGPU register tile order correlator output format with N_INPUTS.
// NB: N_INPUTS = N_STATION * N_POL
#define REGTILE_CHAN_LENGTH (4 * 4 * N_INPUTS/4 * (N_INPUTS/4+1) / 2)
static inline size_t regtile_chan_length()
{
  // Four cross products for each quadrant of 4 input x 4 input tile
  return REGTILE_CHAN_LENGTH;
}

// Returns index into the GPU's register tile ordered output buffer for the
// real component of the cross product of inputs in0 and in1.  Note that in0
// and in1 are input indexes (i.e. 0 based) and often represent antenna and
// polarization by passing (2*ant_idx+pol_idx) as the input number (NB: ant_idx
// and pol_idx are also 0 based).  Return value is valid if in1 >= in0.  The
// corresponding imaginary component is located xgpu_info.matLength words after
// the real component.
static off_t regtile_index(const int in0, const int in1)
{
  const int a0 = in0 >> 1;
  const int a1 = in1 >> 1;
  const int p0 = in0 & 1;
  const int p1 = in1 & 1;
  const int num_words_per_cell = 4;

  // Index within a quadrant
  const int quadrant_index = tri_index(a1/2, a0/2);
  // Quadrant for this input pair
  const int quadrant = 2*(a0&1) + (a1&1);
  // Size of quadrant
  const int quadrant_size = (xgpu_info.nstation/2 + 1) * xgpu_info.nstation/4;
  // Index of cell (in units of cells)
  const int cell_index = quadrant*quadrant_size + quadrant_index;
  //printf("%s: in0=%d, in1=%d, a0=%d, a1=%d, cell_index=%d\n", __FUNCTION__, in0, in1, a0, a1, cell_index);
  // Pol offset
  const int pol_offset = 2*p1 + p0;
  // Word index (in units of words (i.e. floats) of real component
  const int index = (cell_index * num_words_per_cell) + pol_offset;
  return index;
}

// Returns index into a CASPER ordered buffer for the real component of the
// cross product of inputs in0 and in1.  Note that in0 and in1 are input
// indexes (i.e. 0 based) and often represent antenna and polarization by
// passing (2*ant_idx+pol_idx) as the input number (NB: ant_idx ad pol_idx are
// also 0 based).  Return value is valid if in1 >= in0.  The corresponding
// imaginary component is located in the word immediately following the real
// component.  A casper ordered buffer consists of four complex values for each
// pair of input pairs.  Thus, the number of complex values in a casper ordered
// buffer are: 4 * (N/2 * (N/2 + 1)) / 2 = N * (N/2 + 1)

// For integration buffers, change the ordering to baselines x channels x stokes
// to make it easier to send the packets after integration.

/*  Each baseline (i.e, ant pair) needs a unique index for encoding the  
 *  baseline both while packetization and while building the integration
 *  buffers. The baseline index function below is an adaption of the 
 *  CASPER index without accounting for the size of each cell or the 
 *  stokes parameters.
 */ 

/*  The integration buffer location for each baseline is obtained by  
 *  multiplying the baseline_index with:
 *  (words_per_cell = 8) * (num_chan_per_x = 384)
 *  Pol offset = 2* (2*(p0^p1) + p0)
 */
static int baseline_index(const int in0, const int in1, const int n)
{ 
  const int a0 = in0 >> 1; 
  const int a1 = in1 >> 1;
  const int delta = a1-a0;
  const int nant_2 = (n/2) / 2;

  // Three cases: top triangle, middle rectangle, bottom triangle
  const int triangle_size = ((nant_2 + 1) * nant_2)/2;
  const int middle_rect_offset = triangle_size;
  const int last_cell_offset = 4*middle_rect_offset - nant_2 - 1;
  int cell_index;

  if(delta > nant_2) {
    // bottom triangle
    cell_index = last_cell_offset - tri_index(nant_2-2-a0, (n/2)-1-a1);
  } else if (a1 < (n/2)/2) {
    // top triangle
    cell_index = tri_index(a1, a0);
  } else {
    // middle rectangle
    cell_index = middle_rect_offset + (a1-nant_2)*(nant_2+1) + (nant_2-delta);
  }

  return cell_index;
}
  
// Redefine casper ordering to place channels next to each other (since 
// buffers are packetized this way. Thus, the number of complex values 
// per baselines are: 4 stokes * N_CHAN_PER_X 
// Total buffersize: 4 * N_CHAN_PER_X * (N/2 * (N/2 + 1)) / 2
#define N_CASPER_COMPLEX_PER_BASELINE  (N_STOKES * N_CHAN_PER_X)

// Lookup table mapping baseline_idx to regtile_idx
static off_t *regtile_idx_map;

static int init_idx_map()
{
  int a0, a1;
  regtile_idx_map = (off_t *)malloc(N_BASELINES * sizeof(off_t));
  if(!regtile_idx_map) {
    return -1;
  }

  fprintf(stderr,"Number of inputs: %d\n",N_INPUTS);

  for(a0=0; a0<N_INPUTS/2; a0++) {
    for(a1 = a0; a1<N_INPUTS/2; a1++) {
      //fprintf(stderr, "(%d,%d) Baseline idx: %d Regtile index:%lld \n",
      //        a0,a1,baseline_index(2*a0, 2*a1, N_INPUTS), regtile_index(2*a0, 2*a1));
      regtile_idx_map[baseline_index(2*a0, 2*a1, N_INPUTS)] = regtile_index(2*a0, 2*a1);
    }
  }
  return 0;
}

static int init(struct hashpipe_thread_args *args)
{   
    // Get sizing parameters
    xgpuInfo(&xgpu_info);
    bytes_per_dump = xgpu_info.triLength * sizeof(Complex);
    packets_per_dump = bytes_per_dump / OUTPUT_BYTES_PER_PACKET;
    printf("bytes_per_dump = %lu\n", bytes_per_dump);

    if(init_idx_map()) {
      return -1;
    }

    fprintf(stderr,"I have an index map now!\n");

    // Success!
    return 0;
}

static int init_buffer(hera_bda_buf_t *bdabuf)
{   
   // Load baseline averaging params from config file
   FILE *fp;  
   int j, bin, a0, a1, inttime;
   unsigned long ctr[N_BDA_BINS] = {0,0,0,0,0};
   hera_int_bin_buf_t *intbuf;

   // count baselines in each bin   
   if((fp=fopen("bda_config.txt","r")) == NULL){
      printf("Cannot open the configuration file.\n");
      exit(1);
   }
   while(fscanf(fp, "%d %d %d", &a0, &a1, &inttime)!=EOF){
      if((inttime == 0) || !CHECK_PWR2(inttime)){
        printf("(%d,%d): Samples to integrate not power of 2!\n",a0,a1);
        exit(1);
      }
      bdabuf->baselines_per_bin[LOG(inttime)]++;
   }

   // Include autos for the no-integration baselines
   bdabuf->baselines_per_bin[0] += N_ANTS;

   // malloc for storing ant pairs
   for(j=0; j<N_BDA_BINS; j++){
     bdabuf->ant_pair_0[j] = (int *)malloc(bdabuf->baselines_per_bin[j]*sizeof(int));
     bdabuf->ant_pair_1[j] = (int *)malloc(bdabuf->baselines_per_bin[j]*sizeof(int));
   }
   // re-read ant pairs for each integration bin
   rewind(fp);
   while(fscanf(fp, "%d %d %d", &a0, &a1, &inttime)!=EOF){
     bin = LOG(inttime);
     bdabuf->ant_pair_0[bin][ctr[bin]] = a0;
     bdabuf->ant_pair_1[bin][ctr[bin]++] = a1; 
   }
   fclose(fp);

   // Include the autos
   for(a0=0;a0<N_ANTS;a0++){
     bdabuf->ant_pair_0[0][ctr[0]] = a0;
     bdabuf->ant_pair_1[0][ctr[0]++] = a0;
   }

   // Initialise header and malloc for data
   bdabuf->send[0] = 1;
   bdabuf->buf->data = NULL;
   bdabuf->buf->header.datsize = 0;

   for(j=1; j<N_BDA_BINS; j++){
     //fprintf(stderr,"Bin: %d Baselines: %ld\n", j, bdabuf->baselines_per_bin[j]);
     bdabuf->send[j] = 0;
     intbuf = bdabuf->buf+j;
     intbuf->header.sam    = 0;
     intbuf->header.totsam = pow(2,j);
     intbuf->header.datsize = bdabuf->baselines_per_bin[j]*N_COMPLEX_PER_BASELINE*2*sizeof(uint32_t);
     fprintf(stderr,"Init: Size of Buf: %d is %llu\n",j,intbuf->header.datsize);
     intbuf->data = malloc(intbuf->header.datsize);     
     if (!intbuf->data)
       return -1;
     memset(intbuf->data, 0, intbuf->header.datsize);
   }
return 1;
}

#define ELAPSED_NS(start,stop) \
  (((int64_t)stop.tv_sec-start.tv_sec)*1000*1000*1000+(stop.tv_nsec-start.tv_nsec))

static void *run(hashpipe_thread_args_t * args)
{
   // Local aliases to shorten access to args fields
   // Our input buffer happens to be a paper_ouput_databuf
   paper_output_databuf_t *db = (paper_output_databuf_t *)args->ibuf;
   hashpipe_status_t st = args->st;
   const char * status_key = args->thread_desc->skey;

   // Initialize buffers to store averaged baselines
   fprintf(stderr, "Initializing buffers..\n");
   hera_bda_buf_t bdabuf;
   init_buffer(&bdabuf); 

   // Setup socket and message structures
   int sockfd;
   unsigned int xengine_id = 0;
   struct timespec packet_delay = {
     .tv_sec = 0,
     .tv_nsec = PACKET_DELAY_NS
   };

   fprintf(stderr, "Created socket\n");

   hashpipe_status_lock_safe(&st);
   hgetu4(st.buf, "XID", &xengine_id); // No change if not found
   hputu4(st.buf, "XID", xengine_id);
   hputu4(st.buf, "OUTDUMPS", 0);
   hashpipe_status_unlock_safe(&st);

   pkt_t pkt;
   pkt.hdr.xeng_id = XENG_ID(xengine_id);
   pkt.hdr.payload_len = PAYLOAD_LEN(OUTPUT_BYTES_PER_PACKET);

   // TODO Get catcher hostname and port from somewhere

#define stringify2(x) #x
#define stringify(x) stringify2(x)

   // Open socket
   sockfd = open_udp_socket("localhost", stringify(CATCHER_PORT));
   if(sockfd == -1) {
       hashpipe_error(__FUNCTION__, "error opening socket");
       pthread_exit(NULL);
   }

#ifdef TEST_INDEX_CALCS
    int i, j;
    off_t o;
    for(i=0; i<32; i++) {
      for(j=i; j<32; j++) {
        regtile_index(2*i, 2*j);
      }
    }
    for(i=0; i<32; i++) {
      for(j=i; j<32; j++) {
        o = baseline_index(2*i, 2*j, N_INPUTS);
        fprintf(stdout, "%d, %d, %d, %ld\n", i, j, (int) o, 
                                    (long int)regtile_idx_map[o]);
      }
    }
#endif

   /* Main loop */
   int rv;
   int casper_chan, gpu_chan;
   unsigned long bl;
   int j,pol;
   unsigned int dumps = 0;
   int block_idx = 0;
   struct timespec start, stop;
   struct timespec pkt_start, pkt_stop;
   pktdata_t re, im;    //pktdata_t is 32bits
   uint32_t nbytes = 0;
   int offset = 0;
   unsigned long long datoffset = 0;
   int ant0, ant1;
   int catcher_bl_id = 0;
   unsigned long long idx_baseline; 
   off_t idx_regtile;
   hera_int_bin_buf_t *intbuf;

   fprintf(stderr, "Starting main loop...\n");

   /* Test statements */
   fprintf(stderr,"Total number of bins: %d\n",N_BDA_BINS);
   for(j=1; j<N_BDA_BINS; j++){
     fprintf(stderr,"Bin:%d\n",j);
     intbuf = &(bdabuf.buf[j-1]);
     fprintf(stderr,"Integration Bin: %d Datsize: %llu\n",j,intbuf->header.datsize);
   }
   fprintf(stdout, "Size of packet header: %ld\n",sizeof(pkt.hdr));
   fprintf(stdout, "Size of matrix: %lld\n", xgpu_info.matLength);
   fprintf(stdout, "Channels per X-eng: %d\n", N_CHAN_PER_X);

   while (run_threads()) {

     hashpipe_status_lock_safe(&st);
     hputs(st.buf, status_key, "waiting");
     hashpipe_status_unlock_safe(&st);

     // Wait for new block to be filled
     while ((rv=paper_output_databuf_wait_filled(db, block_idx))
           != HASHPIPE_OK) {
       if (rv==HASHPIPE_TIMEOUT) {
           hashpipe_status_lock_safe(&st);
           hputs(st.buf, status_key, "blocked");
           hashpipe_status_unlock_safe(&st);
           continue;
       } else {
           hashpipe_error(__FUNCTION__, "error waiting for filled databuf");
           pthread_exit(NULL);
           break;
       }
     }

     clock_gettime(CLOCK_MONOTONIC, &start);

     // Note processing status, current input block
     hashpipe_status_lock_safe(&st);
     hputs(st.buf, status_key, "processing");
     hputi4(st.buf, "OUTBLKIN", block_idx);
     hashpipe_status_unlock_safe(&st);
     
     // Unpack and convert in packet sized chunks
     // output data in order: baseline x chan x stokes (slowest to fastest varying)
     pktdata_t *pf_re  = db->block[block_idx].data;
     pktdata_t *pf_im  = db->block[block_idx].data + xgpu_info.matLength;
     pktdata_t *p_out = pkt.data;

     /* ------------------------------------------------ */
     /* Loop through baselines and add/send the packets  */
     /* ------------------------------------------------ */

     for(j=0; j<N_BDA_BINS; j++){ //intbuf loop
       intbuf = &(bdabuf.buf[j]);
       intbuf->header.mcnt = db->block[block_idx].header.mcnt;

       if (bdabuf.send[j]){
         fprintf(stderr,"Bin: %d: Send\n",j);
         pkt.hdr.timestamp = TIMESTAMP(intbuf->header.mcnt);
         offset = 0; nbytes = 0;
         p_out = pkt.data; 
         clock_gettime(CLOCK_MONOTONIC, &pkt_start);
 
         for(bl=0; bl<bdabuf.baselines_per_bin[j]; bl++){
           ant0 = bdabuf.ant_pair_0[j][bl];    // ordering set by config file 
           ant1 = bdabuf.ant_pair_1[j][bl]; 
           pkt.hdr.ant1 = ant1;
           pkt.hdr.ant0 = ant0;
           idx_baseline = baseline_index(2*ant0, 2*ant1, N_INPUTS); 
           idx_regtile = regtile_idx_map[idx_baseline];
           //idx_regtile = regtile_index(2*ant0, 2*ant1);
           pkt.hdr.offset = OFFSET(0);
           offset = 0;
           pkt.hdr.baseline_id = BASELINE_ID(baseline_id++);//BASELINE_ID(idx_baseline);
           
           for(casper_chan=0; casper_chan<N_CHAN_PER_X; casper_chan++){
             // This used to de-interleave channels.  De-interleaving is no longer
             // needed, but we choose to continue maintaining the distinction
             // between casper_chan and gpu_chan.
             gpu_chan = casper_chan;

             for(pol=0; pol<N_STOKES; pol++){
               re = pf_re[gpu_chan*REGTILE_CHAN_LENGTH+idx_regtile+pol];
               im = pf_im[gpu_chan*REGTILE_CHAN_LENGTH+idx_regtile+pol];
               if (j==0){  
                 *p_out++ = CONVERT(re);
                 *p_out++ = CONVERT(-im);
               }else{ 
                 datoffset = hera_int_bin_buf_data_idx(bl, gpu_chan, pol);
                 *p_out++ = CONVERT(intbuf->data[datoffset] + re);
                 *p_out++ = CONVERT(intbuf->data[datoffset+1] -im);
               }
               nbytes += 2*sizeof(pktdata_t);
               
               if(nbytes%OUTPUT_BYTES_PER_PACKET == 0){
                 int bytes_sent = send(sockfd, &pkt, 
                                       sizeof(pkt.hdr)+OUTPUT_BYTES_PER_PACKET, 0); 

                 if(bytes_sent == -1){
                   // Send all packets even if catcher is not listening (i.e. we
                   // we get a connection refused error), but abort sending this
                   // dump if we get any other error.
                   if(errno != ECONNREFUSED){
                     perror("send");
                     // Update stats
                     hashpipe_status_lock_safe(&st);
                     hgetu4(st.buf, "OUTDUMPS", &dumps);
                     hputu4(st.buf, "OUTDUMPS", ++dumps);
                     hputr4(st.buf, "OUTSECS", 0.0);
                     hputr4(st.buf, "OUTMBPS", 0.0);
                     hashpipe_status_unlock_safe(&st);
                     // Break out of both for loops
                     goto done_sending;
                   }
                 }else if(bytes_sent != sizeof(pkt.hdr)+OUTPUT_BYTES_PER_PACKET) {
                   printf("only sent %d of %lu bytes!!!\n", bytes_sent, 
                           sizeof(pkt.hdr)+OUTPUT_BYTES_PER_PACKET);
                 }

                 // Delay to prevent overflowing network TX queue
                 clock_gettime(CLOCK_MONOTONIC, &pkt_stop);
                 packet_delay.tv_nsec = PACKET_DELAY_NS - ELAPSED_NS(pkt_start, pkt_stop);
                 if(packet_delay.tv_nsec > 0 && packet_delay.tv_nsec < 1000*1000*1000){
                   nanosleep(&packet_delay, NULL);
                 }
               
                 // Setup for next packet
                 p_out = pkt.data;
                 pkt_start = pkt_stop;
                 offset++;
                 pkt.hdr.offset = OFFSET(offset);
               } // done sending
             } // pol
           } // chan
         } // baseline
         
         if(j!=0) bdabuf.send[j] = 0; 
         intbuf->header.sam = 0;
         memset(intbuf->data, 0, intbuf->header.datsize);
          
       }else{
         fprintf(stderr,"Bin: %d: Copy\n",j);
         for(bl=0; bl<bdabuf.baselines_per_bin[j]; bl++){
           ant0 = bdabuf.ant_pair_0[j][bl]; 
           ant1 = bdabuf.ant_pair_1[j][bl];
           //fprintf(stderr,"Integration Bin: %d: Baseline Num: %ld  Ants:(%d,%d) \n",j,bl,ant0,ant1);
           idx_baseline = baseline_index(2*ant0, 2*ant1, N_INPUTS); 
           idx_regtile = regtile_idx_map[idx_baseline];
           //fprintf(stderr,"Integration Bin: %d: baseline_idx: %lld regtile_idx:%lld\n", 
           //                                  j, idx_baseline, idx_regtile);          
           //fprintf(stderr, "Integration Bin: %d: Data Offset: %lld Dat Size:%lld\n",
           //                                  j,datoffset,intbuf->header.datsize);
 
           for(casper_chan=0; casper_chan<N_CHAN_PER_X; casper_chan++){
             //fprintf(stderr, "Integration: %d: Chan: %d\n", j, casper_chan);
             gpu_chan = casper_chan;
             for(pol=0; pol<N_STOKES; pol++){
               re = pf_re[(gpu_chan*REGTILE_CHAN_LENGTH)+idx_regtile+pol];
               im = pf_im[(gpu_chan*REGTILE_CHAN_LENGTH)+idx_regtile+pol];
               datoffset = hera_int_bin_buf_data_idx(bl, gpu_chan, pol);
               intbuf->data[datoffset] += re;
               intbuf->data[datoffset+1] += -im;
             }
           }
         } fprintf(stderr,"Baselines processed: %ld Total baselines: %ld\n", 
                   bl,bdabuf.baselines_per_bin[j]);
 
        intbuf->header.sam += 1;
        if ((intbuf->header.sam+1) == intbuf->header.totsam)
          bdabuf.send[j] = 1;
       }
     }// end intbuf loop

     clock_gettime(CLOCK_MONOTONIC, &stop);

     hashpipe_status_lock_safe(&st);
     hgetu4(st.buf, "OUTDUMPS", &dumps);
     hputu4(st.buf, "OUTDUMPS", ++dumps);
     hputu4(st.buf, "OUTBYTES", nbytes);
     hputr4(st.buf, "OUTSECS", (float)ELAPSED_NS(start,stop)/1e9);
     hputr4(st.buf, "OUTMBPS", (1e3*8*nbytes)/ELAPSED_NS(start,stop));
     hashpipe_status_unlock_safe(&st);

done_sending:

     // Mark block as free
     paper_output_databuf_set_free(db, block_idx);

     // Setup for next block
     block_idx = (block_idx + 1) % db->header.n_block;

     /* Will exit if thread has been cancelled */
     pthread_testcancel();
  }

    // Thread success!
    return NULL;
}

static hashpipe_thread_desc_t gpu_output_thread = {
    name: "hera_gpu_output_thread",
    skey: "OUTSTAT",
    init: init,
    run:  run,
    ibuf_desc: {paper_output_databuf_create},
    obuf_desc: {NULL}
};

static __attribute__((constructor)) void ctor()
{
  register_hashpipe_thread(&gpu_output_thread);
}

//     /* ---------------------------------------------------- */
//     /* Loop through 1-sample baselines and send the packets */
//     /* ---------------------------------------------------- */
//       
//     fprintf(stderr, "Processing 1-sample baselines\n");
//
//     // Update header's timestamp for this dump.
//     pkt.hdr.timestamp = TIMESTAMP(db->block[block_idx].header.mcnt);
//     pkt.hdr.offset = OFFSET(0);
//     clock_gettime(CLOCK_MONOTONIC, &pkt_start);
//
//     // All stokes are sent in adjacent words
//     for(bl=0; bl < bdabuf.baselines_per_bin[0]; bl++){
//       ant0 = bdabuf.ant_pair_0[0][bl];
//       ant1 = bdabuf.ant_pair_1[0][bl];
//       idx_baseline = baseline_index(2*ant0, 2*ant1, N_INPUTS);
//       idx_regtile = regtile_idx_map[idx_baseline];
//       //idx_regtile = regtile_index(2*ant0, 2*ant1);
//       pkt.hdr.baseline_id = BASELINE_ID(idx_baseline);
//       offset = 0;
//       pkt.hdr.offset = OFFSET(0);
//
//       for(casper_chan=0; casper_chan<N_CHAN_PER_X; casper_chan++){
//         // This used to de-interleave channels.  De-interleaving is no longer
//         // needed, but we choose to continue maintaining the distinction
//         // between casper_chan and gpu_chan.
//         gpu_chan = casper_chan;
//
//         for(pol=0; pol<N_STOKES; pol++){
//           re = pf_re[(gpu_chan*REGTILE_CHAN_LENGTH)+idx_regtile+pol];
//           im = pf_im[(gpu_chan*REGTILE_CHAN_LENGTH)+idx_regtile+pol];
//           *p_out++ = CONVERT(re);
//           *p_out++ = CONVERT(-im); // Conjugate to match downstream expectations.
//           nbytes += 2*sizeof(pktdata_t);
//
//           if(nbytes % OUTPUT_BYTES_PER_PACKET == 0) {
//             int bytes_sent = send(sockfd, &pkt, 
//                                   sizeof(pkt.hdr)+OUTPUT_BYTES_PER_PACKET, 0); 
//
//             if(bytes_sent == -1){
//               // Send all packets even if catcher is not listening (i.e. we
//               // we get a connection refused error), but abort sending this
//               // dump if we get any other error.
//               if(errno != ECONNREFUSED){
//                 perror("send");
//                 // Update stats
//                 hashpipe_status_lock_safe(&st);
//                 hgetu4(st.buf, "OUTDUMPS", &dumps);
//                 hputu4(st.buf, "OUTDUMPS", ++dumps);
//                 hputr4(st.buf, "OUTSECS", 0.0);
//                 hputr4(st.buf, "OUTMBPS", 0.0);
//                 hashpipe_status_unlock_safe(&st);
//                 // Break out of both for loops
//                 goto done_sending;
//               }
//             }else if(bytes_sent != sizeof(pkt.hdr)+OUTPUT_BYTES_PER_PACKET) {
//               printf("only sent %d of %lu bytes!!!\n", bytes_sent, 
//                       sizeof(pkt.hdr)+OUTPUT_BYTES_PER_PACKET);
//             }
//
//             // Delay to prevent overflowing network TX queue
//             clock_gettime(CLOCK_MONOTONIC, &pkt_stop);
//             packet_delay.tv_nsec = PACKET_DELAY_NS - ELAPSED_NS(pkt_start, pkt_stop);
//             if(packet_delay.tv_nsec > 0 && packet_delay.tv_nsec < 1000*1000*1000){
//               nanosleep(&packet_delay, NULL);
//             }
//
//             // Setup for next packet
//             p_out = pkt.data;
//             pkt_start = pkt_stop;
//             offset++;
//             pkt.hdr.offset = OFFSET(offset);
//           } // done sending
//         } // pol loop
//       } // chan loop
//     } // baseline loop


