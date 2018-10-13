#ifndef _PAPER_DATABUF_H
#define _PAPER_DATABUF_H

#include <stdint.h>
#include "hashpipe_databuf.h"
#include "config.h"

// Determined by F engine
#define N_CHAN_TOTAL_GENERATED (8192)
#define FENG_SAMPLE_RATE (500000000)
#define N_CHAN_TOTAL 6144
#define N_FENGINES   192
#define N_CHAN_PER_F N_CHAN_TOTAL

// Number of separate X-engines which deal with
// alternate time chunks
#define TIME_DEMUX 2

// Determined by F engine packetizer
#define N_INPUTS_PER_PACKET  6
#define N_CHAN_PER_PACKET    384
#define N_TIME_PER_PACKET    2
// N_BYTES_PER_PACKET excludes header!
#define N_BYTES_PER_PACKET  (N_INPUTS_PER_PACKET*N_CHAN_PER_PACKET*N_TIME_PER_PACKET)

// X engine sizing (from xGPU)
#define N_ANTS               XGPU_NSTATION
#define N_INPUTS          (2*XGPU_NSTATION)
#define N_TIME_PER_BLOCK     XGPU_NTIME
#define N_CHAN_PER_X         XGPU_NFREQUENCY

// Derived from above quantities
//#define N_SUB_BLOCKS_PER_INPUT_BLOCK (N_TIME_PER_BLOCK / N_TIME_PER_PACKET)
//#define N_SUB_BLOCKS_PER_INPUT_BLOCK N_TIME_PER_BLOCK
//#define N_SUB_BLOCKS_PER_INPUT_BLOCK (N_TIME_PER_BLOCK / 2048)
#define N_BYTES_PER_BLOCK            (N_TIME_PER_BLOCK * N_CHAN_PER_X * N_INPUTS)
#define N_PACKETS_PER_BLOCK          (N_BYTES_PER_BLOCK / N_BYTES_PER_PACKET)
#define N_PACKETS_PER_BLOCK_PER_F    (N_PACKETS_PER_BLOCK * N_INPUTS_PER_PACKET / 2 / N_FENGINES)

// Validate packet dimensions
#if    N_BYTES_PER_PACKET != (N_TIME_PER_PACKET*N_CHAN_PER_PACKET*N_INPUTS_PER_PACKET)
#error N_BYTES_PER_PACKET != (N_TIME_PER_PACKET*N_CHAN_PER_PACKET*N_INPUTS_PER_PACKET)
#endif

#define N_FLUFFED_BYTES_PER_BLOCK  ((N_PACKETS_PER_BLOCK * N_BYTES_PER_PACKET) * 2)
#define N_FLUFFED_WORDS_PER_BLOCK (N_FLUFFED_BYTES_PER_BLOCK / 8) 

// Number of floats in xGPU's "register tile order" output matrix.
#define N_OUTPUT_MATRIX (2 * N_CHAN_PER_X * (N_INPUTS/2 + 2) * N_INPUTS)

#define PAGE_SIZE (4096)
#define CACHE_ALIGNMENT (128)

// Correlator Output parameters
#define CATCHER_PORT 10000
#define OUTPUT_BYTES_PER_PACKET (4096)
#define CATCHER_N_BLOCKS 2
#define CATCHER_CHAN_SUM 4
#define VIS_MATRIX_ENTRIES (N_CHAN_TOTAL * (N_INPUTS * ((N_INPUTS>>1) + 1)))
#define VIS_MATRIX_ENTRIES_PER_CHAN (N_INPUTS * ((N_INPUTS>>1) + 1))
#define PACKETS_PER_VIS_MATRIX ((8L*TIME_DEMUX*VIS_MATRIX_ENTRIES) / OUTPUT_BYTES_PER_PACKET)
#define N_STOKES 4


// The HERA correlator is based largely on the PAPER correlator.  The main 
// difference will be in the F engines.  The ROACH2 based F engines are being
// replaced by SNAP based F engines.  Because SNAPs have a different number of
// inputs, the F engine packet format must change. 
// Even though the X engines are the same, the network thread and the
// fluffing threads for HERA are somewhat different.  The data buffer format
// for the GPU thread's input buffer is exactly the same.  The fluffing
// thread's input buffer is therefore sized the same (i.e. half the size of the
// GPU input buffer), but the format is different since the packets stored by
// the network thread will be different size/format.  Also changed is the
// number of time samples per packet i(and therefore per mcount).

/*
 * INPUT BUFFER STRUCTURES
 */

#define N_INPUT_BLOCKS 4
#ifndef N_DEBUG_INPUT_BLOCKS
#define N_DEBUG_INPUT_BLOCKS 0
#endif

typedef struct paper_input_header {
  int64_t good_data; // functions as a boolean, 64 bit to maintain word alignment
  uint64_t mcnt;     // mcount of first packet
} paper_input_header_t;

typedef uint8_t paper_input_header_cache_alignment[
  CACHE_ALIGNMENT - (sizeof(paper_input_header_t)%CACHE_ALIGNMENT)
];

// For HERA:
//
// * HERA net thread output
// * HERA fluffer thread input
// * Ordered sequence of packets in "data" field:
//
//   +-- mcount
//   |       +-- Antenna
//   |       |       |<-packets->|
//   V       V       |           |
//   m       a       |c        t |
//   ==      ==      |==       ==|
//   m0 }--> a0 }--> |c0 }---> t0|\.
//   m1      a1      |c1       t1| \.
//   m2      a2      |c2       t2|  \ chan 0
//   :       :       |:        : |  / packet
//   :       :       |:        : | /.
//   :       :       |:        : |/.
//   :       :       +---- ------+
//   :       :       |Nc/2 }-> t0|\.
//   :       :       |:        t1| \.
//   :       :       |:        t2|  \ chan Nc/2
//   :       :       |:        : |  / packet
//   :       :       |:        : | /.
//   :       :       |:        : |/.
//   ==      ==       ===      ==
//   Nm      Na       Nc       Nt
//
//   Each time sample is 2 bytes comprising two 4b+4b complex values.
//   Note that the bits may not all be contiguous to facilitate fluffing.
//
// m = packet's mcount - block's first mcount
// a = packet's antenna number
// c = channel within within X engine
// t = time sample within channel (0 to Nt-1)

#define Nm (N_TIME_PER_BLOCK/N_TIME_PER_PACKET)
#define Na N_ANTS
#define Nc N_CHAN_PER_X
#define Nt N_TIME_PER_PACKET
// Number of pols = 2
#define Np 2
// Number of complex elements = 2
#define Nx 2

// HERA mcount index = m * Na * Nc * Nt * HERA_N_INPUTS_PER_PACKET / sizeof(uint64_t)
// HERA ant index = a * Nc * Nt * HERA_N_INPUTS_PER_PACKET / sizeof(uint64_t)
// HERA chan index = c * Nt * HERA_N_INPUTS_PER_PACKET / sizeof(uint64_t)
//
// Computes paper_input_databuf_t data word (uint64_t) offset for complex data
// word corresponding to the given parameters for HERA F engine packets.

//#define paper_input_databuf_data_idx(m,a,c,t) \ //foo
//  (((((m * Na + a) * (Nc) + c)*Nt + t) * 2) / sizeof(uint64_t))
#define paper_input_databuf_data_idx(m,a,c,t) \
  ((((m) * Na*Nc*Nt*Np) + ((a) * Nc*Nt*Np) + ((c)*Nt*Np) + ((t)*Np)) / sizeof(uint64_t))

#define paper_input_databuf_data_idx8(m,a,p,c,t) \
  (((m) * Na*Nc*Nt*Np) + ((a) * Nc*Nt*Np) + ((c)*Nt*Np) + ((t)*Np) + (p))

#define paper_input_databuf_data_idx256(m,a,c,t) \
  ((((((m) * Na) + (a))*Nc + (c))*Nt + (t))*Np / sizeof(__m256i))
  //((((m) * Na*Nc*Nt*Np) + ((a) * Nc*Nt*Np) + ((c)*Nt*Np) + ((t)*Np)) / sizeof(__m256i))
//TODO  (((((m * Na + a) * (Nc) + c)*Nt + t) * N_INPUTS_PER_PACKET) / sizeof(uint64_t))

typedef struct paper_input_block {
  paper_input_header_t header;
  paper_input_header_cache_alignment padding; // Maintain cache alignment
  uint64_t data[N_BYTES_PER_BLOCK/sizeof(uint64_t)];
} paper_input_block_t;

// Used to pad after hashpipe_databuf_t to maintain cache alignment
typedef uint8_t hashpipe_databuf_cache_alignment[
  CACHE_ALIGNMENT - (sizeof(hashpipe_databuf_t)%CACHE_ALIGNMENT)
];

typedef struct paper_input_databuf {
  hashpipe_databuf_t header;
  hashpipe_databuf_cache_alignment padding; // Maintain cache alignment
  paper_input_block_t block[N_INPUT_BLOCKS+N_DEBUG_INPUT_BLOCKS];
} paper_input_databuf_t;

/*
 * GPU INPUT BUFFER STRUCTURES
 */

#define N_GPU_INPUT_BLOCKS 4

// == paper_gpu_input_databuf_t ==
//
// * Fluffer thread output
// * GPU thread input
// GPU input buffer is input[time/4][channel][station][pol][complexity][time%4]
// * Multidimensional array in "data" field:
//
//   +--time   +--chan +--fid  +--complexity +--time%4
//   |         |       |       |             |
//   V         V       V       V             V
//   m         c       f       r             i
//   ==        ==      ==
//   m0 }----> c0 }--> f0 }--> real }------> i0
//   m1        c1      f1      imag          i1
//   m2        c2      f2                    i2
//   :         :       :                     i3
//   ==        ==      ==      ==            ==
//   Nm*Nt/4   Nc      Na      2             t%4
//
//   Each group of eight bytes contains four timesamples for
//   the real and imaginary part of a single ant-pol-chan

// Returns word (uint64_t) offset for real input data word
// corresponding to the given parameters.  Corresponding imaginary data
// word is 4 bytes later. Only valid for t a multiple of 4

//#define paper_gpu_input_databuf_data_idx(m,a,t,c) \ //foo
//  (((m)*Nt*Nc*Na) + (((t)/4)*Nc*Na*4) + ((c)*Na*4) + (a)*4)

#define paper_gpu_input_databuf_data_idx(t, c, a, p, x) \
  (((((t)>>2)*4*Nx*Np*Na*Nc) + ((c)*4*Nx*Np*Na) + ((a)*4*Np*Nx) + ((p)*4*Nx) + ((x)*4) + ((t)%4)) / sizeof(uint64_t))

#define paper_gpu_input_databuf_data_idx8(t, c, a, p, x) \
  ((((t)>>2)*4*Nx*Np*Na*Nc) + ((c)*4*Nx*Np*Na) + ((a)*4*Np*Nx) + ((p)*4*Nx) + ((x)*4) + ((t)%4))

#define paper_gpu_input_databuf_data_idx256(t, c, a, p, x) \
  (((((t)>>2)*4*Nx*Np*Na*Nc) + ((c)*4*Nx*Np*Na) + ((a)*4*Np*Nx) + ((p)*4*Nx) + ((x)*4) + ((t)%4)) / sizeof(__m256i))

#define paper_gpu_input_databuf_data_idx_tca_256(t, c, a) \
  (((((((t)*Nc) + 4*(c))*Na + 4*(a))*Nx*Np) + (t)%4) / sizeof(__m256i))
  //(((((((t)>>2)*Nc) + (c))*Na + (a))*Nx*4*Np + ((t)%4)) / sizeof(__m256i))

//#define paper_gpu_input_databuf_data_idx256(m,a,t,c) \ //foo
//  (4*((m*Nt*Nc*Na) + (t*Nc*Na) + (c*Na) + a) / sizeof(__m256i))

typedef struct paper_gpu_input_block {
  paper_input_header_t header;
  paper_input_header_cache_alignment padding; // Maintain cache alignment
  uint64_t data[(2*N_BYTES_PER_BLOCK/sizeof(uint64_t))];
} paper_gpu_input_block_t;

typedef struct paper_gpu_input_databuf {
  hashpipe_databuf_t header;
  hashpipe_databuf_cache_alignment padding; // Maintain cache alignment
  paper_gpu_input_block_t block[N_GPU_INPUT_BLOCKS];
} paper_gpu_input_databuf_t;

/*
 * OUTPUT BUFFER STRUCTURES
 */

#define N_OUTPUT_BLOCKS 2

typedef struct paper_output_header {
  uint64_t mcnt;
  uint64_t flags[(N_CHAN_PER_X+63)/64];
} paper_output_header_t;

typedef uint8_t paper_output_header_cache_alignment[
  CACHE_ALIGNMENT - (sizeof(paper_output_header_t)%CACHE_ALIGNMENT)
];

typedef struct paper_output_block {
  paper_output_header_t header;
  paper_output_header_cache_alignment padding; // Maintain cache alignment
  int32_t data[N_OUTPUT_MATRIX];
} paper_output_block_t;

typedef struct paper_output_databuf {
  hashpipe_databuf_t header;
  hashpipe_databuf_cache_alignment padding; // Maintain cache alignment
  paper_output_block_t block[N_OUTPUT_BLOCKS];
} paper_output_databuf_t;

/*
 * CATCHER BUFFER STRUCTURES
 */

typedef struct hera_catcher_input_header{
  uint64_t mcnt;
  uint64_t good_data;
  //uint8_t flags[VIS_MATRIX_ENTRIES*TIME_DEMUX];
} hera_catcher_input_header_t;

typedef uint8_t hera_catcher_input_header_cache_alignment[
  CACHE_ALIGNMENT - (sizeof(hera_catcher_input_header_t)%CACHE_ALIGNMENT)
];

typedef struct hera_catcher_input_block {
  hera_catcher_input_header_t header;
  hera_catcher_input_header_cache_alignment padding; // Maintain cache alignment
  uint64_t data[VIS_MATRIX_ENTRIES*TIME_DEMUX];
} hera_catcher_input_block_t;

typedef struct hera_catcher_input_databuf {
  hashpipe_databuf_t header;
  hashpipe_databuf_cache_alignment padding; // Maintain cache alignment
  hera_catcher_input_block_t block[CATCHER_N_BLOCKS];
} hera_catcher_input_databuf_t;

#define hera_catcher_input_databuf_idx32(x, o) \
  (2L*(VIS_MATRIX_ENTRIES_PER_CHAN * N_CHAN_PER_X*x) + (o>>2))

/*
 * INPUT BUFFER FUNCTIONS
 */

hashpipe_databuf_t *paper_input_databuf_create(int instance_id, int databuf_id);

static inline paper_input_databuf_t *paper_input_databuf_attach(int instance_id, int databuf_id)
{
    return (paper_input_databuf_t *)hashpipe_databuf_attach(instance_id, databuf_id);
}

static inline int paper_input_databuf_detach(paper_input_databuf_t *d)
{
    return hashpipe_databuf_detach((hashpipe_databuf_t *)d);
}

static inline void paper_input_databuf_clear(paper_input_databuf_t *d)
{
    hashpipe_databuf_clear((hashpipe_databuf_t *)d);
}

static inline int paper_input_databuf_block_status(paper_input_databuf_t *d, int block_id)
{
    return hashpipe_databuf_block_status((hashpipe_databuf_t *)d, block_id);
}

static inline int paper_input_databuf_total_status(paper_input_databuf_t *d)
{
    return hashpipe_databuf_total_status((hashpipe_databuf_t *)d);
}


int paper_input_databuf_wait_free(paper_input_databuf_t *d, int block_id);

int paper_input_databuf_busywait_free(paper_input_databuf_t *d, int block_id);

int paper_input_databuf_wait_filled(paper_input_databuf_t *d, int block_id);

int paper_input_databuf_busywait_filled(paper_input_databuf_t *d, int block_id);

int paper_input_databuf_set_free(paper_input_databuf_t *d, int block_id);

int paper_input_databuf_set_filled(paper_input_databuf_t *d, int block_id);

/*
 * CATCHER BUFFER FUNCTIONS
 */

hashpipe_databuf_t *hera_catcher_input_databuf_create(int instance_id, int databuf_id);

static inline hera_catcher_input_databuf_t *hera_catcher_input_databuf_attach(int instance_id, int databuf_id)
{
    return (hera_catcher_input_databuf_t *)hashpipe_databuf_attach(instance_id, databuf_id);
}

static inline int hera_catcher_input_databuf_detach(hera_catcher_input_databuf_t *d)
{
    return hashpipe_databuf_detach((hashpipe_databuf_t *)d);
}

static inline void hera_catcher_input_databuf_clear(hera_catcher_input_databuf_t *d)
{
    hashpipe_databuf_clear((hashpipe_databuf_t *)d);
}

static inline int hera_catcher_input_databuf_block_status(hera_catcher_input_databuf_t *d, int block_id)
{
    return hashpipe_databuf_block_status((hashpipe_databuf_t *)d, block_id);
}

static inline int hera_catcher_input_databuf_total_status(hera_catcher_input_databuf_t *d)
{
    return hashpipe_databuf_total_status((hashpipe_databuf_t *)d);
}

int hera_catcher_input_databuf_wait_free(hera_catcher_input_databuf_t *d, int block_id);

int hera_catcher_input_databuf_busywait_free(hera_catcher_input_databuf_t *d, int block_id);

int hera_catcher_input_databuf_wait_filled(hera_catcher_input_databuf_t *d, int block_id);

int hera_catcher_input_databuf_busywait_filled(hera_catcher_input_databuf_t *d, int block_id);

int hera_catcher_input_databuf_set_free(hera_catcher_input_databuf_t *d, int block_id);

int hera_catcher_input_databuf_set_filled(hera_catcher_input_databuf_t *d, int block_id);


/*
 * GPU INPUT BUFFER FUNCTIONS
 */

hashpipe_databuf_t *paper_gpu_input_databuf_create(int instance_id, int databuf_id);

static inline void paper_gpu_input_databuf_clear(paper_gpu_input_databuf_t *d)
{
    hashpipe_databuf_clear((hashpipe_databuf_t *)d);
}

static inline paper_gpu_input_databuf_t *paper_gpu_input_databuf_attach(int instance_id, int databuf_id)
{
    return (paper_gpu_input_databuf_t *)hashpipe_databuf_attach(instance_id, databuf_id);
}

static inline int paper_gpu_input_databuf_detach(paper_gpu_input_databuf_t *d)
{
    return hashpipe_databuf_detach((hashpipe_databuf_t *)d);
}

static inline int paper_gpu_input_databuf_block_status(paper_gpu_input_databuf_t *d, int block_id)
{
    return hashpipe_databuf_block_status((hashpipe_databuf_t *)d, block_id);
}

static inline int paper_gpu_input_databuf_total_status(paper_gpu_input_databuf_t *d)
{
    return hashpipe_databuf_total_status((hashpipe_databuf_t *)d);
}

static inline int paper_gpu_input_databuf_wait_free(paper_gpu_input_databuf_t *d, int block_id)
{
    return hashpipe_databuf_wait_free((hashpipe_databuf_t *)d, block_id);
}

static inline int paper_gpu_input_databuf_busywait_free(paper_gpu_input_databuf_t *d, int block_id)
{
    return hashpipe_databuf_busywait_free((hashpipe_databuf_t *)d, block_id);
}

static inline int paper_gpu_input_databuf_wait_filled(paper_gpu_input_databuf_t *d, int block_id)
{
    return hashpipe_databuf_wait_filled((hashpipe_databuf_t *)d, block_id);
}

static inline int paper_gpu_input_databuf_busywait_filled(paper_gpu_input_databuf_t *d, int block_id)
{
    return hashpipe_databuf_busywait_filled((hashpipe_databuf_t *)d, block_id);
}

static inline int paper_gpu_input_databuf_set_free(paper_gpu_input_databuf_t *d, int block_id)
{
    return hashpipe_databuf_set_free((hashpipe_databuf_t *)d, block_id);
}

static inline int paper_gpu_input_databuf_set_filled(paper_gpu_input_databuf_t *d, int block_id)
{
    return hashpipe_databuf_set_filled((hashpipe_databuf_t *)d, block_id);
}

/*
 * OUTPUT BUFFER FUNCTIONS
 */

hashpipe_databuf_t *paper_output_databuf_create(int instance_id, int databuf_id);

static inline void paper_output_databuf_clear(paper_output_databuf_t *d)
{
    hashpipe_databuf_clear((hashpipe_databuf_t *)d);
}

static inline paper_output_databuf_t *paper_output_databuf_attach(int instance_id, int databuf_id)
{
    return (paper_output_databuf_t *)hashpipe_databuf_attach(instance_id, databuf_id);
}

static inline int paper_output_databuf_detach(paper_output_databuf_t *d)
{
    return hashpipe_databuf_detach((hashpipe_databuf_t *)d);
}

static inline int paper_output_databuf_block_status(paper_output_databuf_t *d, int block_id)
{
    return hashpipe_databuf_block_status((hashpipe_databuf_t *)d, block_id);
}

static inline int paper_output_databuf_total_status(paper_output_databuf_t *d)
{
    return hashpipe_databuf_total_status((hashpipe_databuf_t *)d);
}

static inline int paper_output_databuf_wait_free(paper_output_databuf_t *d, int block_id)
{
    return hashpipe_databuf_wait_free((hashpipe_databuf_t *)d, block_id);
}

static inline int paper_output_databuf_busywait_free(paper_output_databuf_t *d, int block_id)
{
    return hashpipe_databuf_busywait_free((hashpipe_databuf_t *)d, block_id);
}

static inline int paper_output_databuf_wait_filled(paper_output_databuf_t *d, int block_id)
{
    return hashpipe_databuf_wait_filled((hashpipe_databuf_t *)d, block_id);
}

static inline int paper_output_databuf_busywait_filled(paper_output_databuf_t *d, int block_id)
{
    return hashpipe_databuf_busywait_filled((hashpipe_databuf_t *)d, block_id);
}

static inline int paper_output_databuf_set_free(paper_output_databuf_t *d, int block_id)
{
    return hashpipe_databuf_set_free((hashpipe_databuf_t *)d, block_id);
}

static inline int paper_output_databuf_set_filled(paper_output_databuf_t *d, int block_id)
{
    return hashpipe_databuf_set_filled((hashpipe_databuf_t *)d, block_id);
}

#endif // _PAPER_DATABUF_H
