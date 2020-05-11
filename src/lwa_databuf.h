#ifndef _LWA_DATABUF_H
#define _LWA_DATABUF_H

#include <stdint.h>
#include "hashpipe_databuf.h"
#include "config.h"

// Determined by F engine
#define N_CHAN_TOTAL_GENERATED (4096)
#define FENG_SAMPLE_RATE (200000000)
#define N_CHAN_TOTAL 3072
#define N_CHAN_PER_F N_CHAN_TOTAL

// Number of F-engine source ports used.
// Using more ports means more IBV receive buffers,
// which potentially increases performance
#define N_SRC_PORTS 8

// Number of separate X-engines which deal with
// alternate time chunks
#define TIME_DEMUX 1

// Determined by F engine packetizer
#define N_INPUTS_PER_PACKET  64
#define N_CHAN_PER_PACKET    96
#define N_TIME_PER_PACKET    1
// N_BYTES_PER_PACKET excludes header!
#define N_BYTES_PER_PACKET  (N_INPUTS_PER_PACKET*N_CHAN_PER_PACKET*N_TIME_PER_PACKET)

// X engine sizing (from xGPU)
#define N_ANTS               XGPU_NSTATION
#define N_INPUTS          (2*XGPU_NSTATION)
#define N_TIME_PER_BLOCK     XGPU_NTIME
#define N_CHAN_PER_X         XGPU_NFREQUENCY

// Derived from above quantities
#define N_FENGINES		((N_INPUTS / N_INPUTS_PER_PACKET) * N_INPUTS_PER_PACKET)
#define N_BYTES_PER_BLOCK            (N_TIME_PER_BLOCK * N_CHAN_PER_X * N_INPUTS)
// Number of bytes per block if we assume an integral number of SNAP boards,
// with all other antennas being fake padding required for XGPU sizing compliance.
#define N_BYTES_PER_BLOCK_NET        (N_TIME_PER_BLOCK * N_CHAN_PER_X * N_FENGINES)
#define N_PACKETS_PER_BLOCK          (N_BYTES_PER_BLOCK_NET / N_BYTES_PER_PACKET)
#define N_PACKETS_PER_BLOCK_PER_F    (N_INPUTS_PER_PACKET * N_PACKETS_PER_BLOCK / N_FENGINES)
// Number of X-engines per time slice. E.g. for LWA: 16.
#define N_XENGINES_PER_TIME (N_CHAN_TOTAL / N_CHAN_PER_X)

// Validate packet dimensions
#if    N_BYTES_PER_PACKET != (N_TIME_PER_PACKET*N_CHAN_PER_PACKET*N_INPUTS_PER_PACKET)
#error N_BYTES_PER_PACKET != (N_TIME_PER_PACKET*N_CHAN_PER_PACKET*N_INPUTS_PER_PACKET)
#endif

#define N_FLUFFED_BYTES_PER_BLOCK  (N_BYTES_PER_BLOCK * 2)
#define N_FLUFFED_WORDS_PER_BLOCK (N_FLUFFED_BYTES_PER_BLOCK / 8) 

// Number of floats in xGPU's "register tile order" output matrix.
#define N_OUTPUT_MATRIX (2 * N_CHAN_PER_X * (N_INPUTS/2 + 2) * N_INPUTS)

#define PAGE_SIZE (4096)
#define CACHE_ALIGNMENT (128)

// This LWA X-Engine is based on the HERA X-Enginee, which in turn
// is based on the LWA X-Engine.
// All these X-engines have the same GPU input buffer format, but different
// network threads, owing to the different packet formats of the systems.
// The difference come about because each telescope has a different number
// of antennas per F-Engine hardware platform.

/*
 * INPUT BUFFER STRUCTURES
 */

#define N_INPUT_BLOCKS 16
#ifndef N_DEBUG_INPUT_BLOCKS
#define N_DEBUG_INPUT_BLOCKS 0
#endif

typedef struct lwa_input_header {
  int64_t good_data; // functions as a boolean, 64 bit to maintain word alignment
  uint64_t mcnt;     // mcount of first packet
} lwa_input_header_t;

typedef uint8_t lwa_input_header_cache_alignment[
  CACHE_ALIGNMENT - (sizeof(lwa_input_header_t)%CACHE_ALIGNMENT)
];

// For LWA:
//
// * LWA net thread output
// * LWA fluffer thread input
// * Ordered sequence of packets in LWA "data" field:
//
// *** Note that for LWA the t is a length-1 array [t0] ***
// *** I.e. N_TIME_PER_PACKET = 1 ***
// *** Unlike for HERA, we explicitly index by pol ID, not ANT,POL
//
//   +-- mcount
//   |       +-- Chan
//   |       |<----- packets ----->|
//   V       V         Pol         |
//   m       |c        |p        t |
//   ==      |==       |==       ==|
//   m0 }--> |c0  }--> |p0 }---> t0|\.
//   m1      |c1       |p1       t1| \.
//   m2      |c2       |p2       t2|  \ chan 0
//   :       |:        |:        : |  / packet
//   :       |:        |:        : | /.
//   :       |:        |:        : |/.
//   :       |:        |:  }---> t0|\.
//   :       |:        |:        t1| \.
//   :       |:        |:        t2|  \ chan Nc/2
//   :       |:        |:        : |  / packet
//   :       |:        |:        : | /.
//   :       |:        |Np-1     : |/.
//   ==       ==        ==       ==
//   Nm       Nc        Np       Nt(=1 for LWA)
//
//   Each time sample is 2 bytes comprising two 4b+4b complex values.
//   Note that the bits may not all be contiguous to facilitate fluffing.
//
// m = packet's mcount - block's first mcount
// a = packet's antenna number
// c = channel within within X engine
// t = time sample within channel (0 to Nt-1)

#define Nm (N_TIME_PER_BLOCK/N_TIME_PER_PACKET)
#define Np N_INPUTS
#define Nc N_CHAN_PER_X
#define Nt N_TIME_PER_PACKET
// Number of complex elements = 2
#define Nx 2

// LWA mcount index = m * Nc * Np * Nt(=1) * LWA_N_INPUTS_PER_PACKET / sizeof(uint64_t)
// LWA pol index = p * Nc * Nt * LWA_N_INPUTS_PER_PACKET / sizeof(uint64_t)
// LWA chan index = c * Nt * LWA_N_INPUTS_PER_PACKET / sizeof(uint64_t)
//
// Computes net_input_databuf_t data word (uint64_t) offset for complex data
// word corresponding to the given parameters for LWA F engine packets.

#define lwa_input_databuf_data_idx(m,p,c,t) \
  ((((m) * Np*Nc*Nt) + ((c) * Np*Nt) + ((p)*Nt) + ((t))) / sizeof(uint64_t))

#define lwa_input_databuf_data_idx8(m,p,c,t) \
  (((m) * Np*Nc*Nt) + ((c) * Np*Nt) + ((p)*Nt) + ((t)))

#define lwa_input_databuf_data_idx256(m,p,c,t) \
  ((((m) * Np*Nc*Nt) + ((c) * Np*Nt) + ((p)*Nt) + ((t))) / sizeof(__m256i))

typedef struct lwa_input_block {
  lwa_input_header_t header;
  lwa_input_header_cache_alignment padding; // Maintain cache alignment
  uint64_t data[N_BYTES_PER_BLOCK/sizeof(uint64_t)];
} lwa_input_block_t;

// Used to pad after hashpipe_databuf_t to maintain cache alignment
typedef uint8_t hashpipe_databuf_cache_alignment[
  CACHE_ALIGNMENT - (sizeof(hashpipe_databuf_t)%CACHE_ALIGNMENT)
];

typedef struct lwa_input_databuf {
  hashpipe_databuf_t header;
  hashpipe_databuf_cache_alignment padding; // Maintain cache alignment
  lwa_input_block_t block[N_INPUT_BLOCKS+N_DEBUG_INPUT_BLOCKS];
} lwa_input_databuf_t;

/*
 * GPU INPUT BUFFER STRUCTURES
 */

#define N_GPU_INPUT_BLOCKS 4

// == lwa_gpu_input_databuf_t ==
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
//   Nm*Nt/4   Nc      Np      2             t%4
//
//   Each group of eight bytes contains four timesamples for
//   the real and imaginary part of a single ant-pol-chan

// Returns word (uint64_t) offset for real input data word
// corresponding to the given parameters.  Corresponding imaginary data
// word is 4 bytes later. Only valid for t a multiple of 4

#define lwa_gpu_input_databuf_data_idx(t, c, p, x) \
  (((((t)>>2)*4*Nx*Np*Nc) + ((c)*4*Nx*Np) + ((p)*4*Nx) + ((x)*4) + ((t)%4)) / sizeof(uint64_t))

#define lwa_gpu_input_databuf_data_idx8(t, c, p, x) \
  ((((t)>>2)*4*Nx*Np*Nc) + ((c)*4*Nx*Np) + ((p)*4*Nx) + ((x)*4) + ((t)%4))

#define lwa_gpu_input_databuf_data_idx256(t, c, p, x) \
  (((((t)>>2)*4*Nx*Np*Nc) + ((c)*4*Nx*Np) + ((p)*4*Nx) + ((x)*4) + ((t)%4)) / sizeof(__m256i))

#define lwa_gpu_input_databuf_data_idx_tca_256(t, c, p) \
  (((((((t)*Nc) + 4*(c))*Np + 4*(p))*Nx*Np) + (t)%4) / sizeof(__m256i))

typedef struct lwa_gpu_input_block {
  lwa_input_header_t header;
  lwa_input_header_cache_alignment padding; // Maintain cache alignment
  uint64_t data[(2*N_BYTES_PER_BLOCK/sizeof(uint64_t))];
} lwa_gpu_input_block_t;

typedef struct lwa_gpu_input_databuf {
  hashpipe_databuf_t header;
  hashpipe_databuf_cache_alignment padding; // Maintain cache alignment
  lwa_gpu_input_block_t block[N_GPU_INPUT_BLOCKS];
} lwa_gpu_input_databuf_t;

/*
 * OUTPUT BUFFER STRUCTURES
 */

#define N_OUTPUT_BLOCKS 2

typedef struct lwa_output_header {
  uint64_t mcnt;
  uint64_t flags[(N_CHAN_PER_X+63) /64];
} lwa_output_header_t;

typedef uint8_t lwa_output_header_cache_alignment[
  CACHE_ALIGNMENT - (sizeof(lwa_output_header_t)%CACHE_ALIGNMENT)
];

typedef struct lwa_output_block {
  lwa_output_header_t header;
  lwa_output_header_cache_alignment padding; // Maintain cache alignment
  int32_t data[N_OUTPUT_MATRIX];
} lwa_output_block_t;

typedef struct lwa_output_databuf {
  hashpipe_databuf_t header;
  hashpipe_databuf_cache_alignment padding; // Maintain cache alignment
  lwa_output_block_t block[N_OUTPUT_BLOCKS];
} lwa_output_databuf_t;


/*
 * INPUT BUFFER FUNCTIONS
 */

hashpipe_databuf_t *lwa_input_databuf_create(int instance_id, int databuf_id);

static inline lwa_input_databuf_t *lwa_input_databuf_attach(int instance_id, int databuf_id)
{
    return (lwa_input_databuf_t *)hashpipe_databuf_attach(instance_id, databuf_id);
}

static inline int lwa_input_databuf_detach(lwa_input_databuf_t *d)
{
    return hashpipe_databuf_detach((hashpipe_databuf_t *)d);
}

static inline void lwa_input_databuf_clear(lwa_input_databuf_t *d)
{
    hashpipe_databuf_clear((hashpipe_databuf_t *)d);
}

static inline int lwa_input_databuf_block_status(lwa_input_databuf_t *d, int block_id)
{
    return hashpipe_databuf_block_status((hashpipe_databuf_t *)d, block_id);
}

static inline int lwa_input_databuf_total_status(lwa_input_databuf_t *d)
{
    return hashpipe_databuf_total_status((hashpipe_databuf_t *)d);
}


int lwa_input_databuf_wait_free(lwa_input_databuf_t *d, int block_id);

int lwa_input_databuf_busywait_free(lwa_input_databuf_t *d, int block_id);

int lwa_input_databuf_wait_filled(lwa_input_databuf_t *d, int block_id);

int lwa_input_databuf_busywait_filled(lwa_input_databuf_t *d, int block_id);

int lwa_input_databuf_set_free(lwa_input_databuf_t *d, int block_id);

int lwa_input_databuf_set_filled(lwa_input_databuf_t *d, int block_id);

/*
 * GPU INPUT BUFFER FUNCTIONS
 */

hashpipe_databuf_t *lwa_gpu_input_databuf_create(int instance_id, int databuf_id);

static inline void lwa_gpu_input_databuf_clear(lwa_gpu_input_databuf_t *d)
{
    hashpipe_databuf_clear((hashpipe_databuf_t *)d);
}

static inline lwa_gpu_input_databuf_t *lwa_gpu_input_databuf_attach(int instance_id, int databuf_id)
{
    return (lwa_gpu_input_databuf_t *)hashpipe_databuf_attach(instance_id, databuf_id);
}

static inline int lwa_gpu_input_databuf_detach(lwa_gpu_input_databuf_t *d)
{
    return hashpipe_databuf_detach((hashpipe_databuf_t *)d);
}

static inline int lwa_gpu_input_databuf_block_status(lwa_gpu_input_databuf_t *d, int block_id)
{
    return hashpipe_databuf_block_status((hashpipe_databuf_t *)d, block_id);
}

static inline int lwa_gpu_input_databuf_total_status(lwa_gpu_input_databuf_t *d)
{
    return hashpipe_databuf_total_status((hashpipe_databuf_t *)d);
}

static inline int lwa_gpu_input_databuf_wait_free(lwa_gpu_input_databuf_t *d, int block_id)
{
    return hashpipe_databuf_wait_free((hashpipe_databuf_t *)d, block_id);
}

static inline int lwa_gpu_input_databuf_busywait_free(lwa_gpu_input_databuf_t *d, int block_id)
{
    return hashpipe_databuf_busywait_free((hashpipe_databuf_t *)d, block_id);
}

static inline int lwa_gpu_input_databuf_wait_filled(lwa_gpu_input_databuf_t *d, int block_id)
{
    return hashpipe_databuf_wait_filled((hashpipe_databuf_t *)d, block_id);
}

static inline int lwa_gpu_input_databuf_busywait_filled(lwa_gpu_input_databuf_t *d, int block_id)
{
    return hashpipe_databuf_busywait_filled((hashpipe_databuf_t *)d, block_id);
}

static inline int lwa_gpu_input_databuf_set_free(lwa_gpu_input_databuf_t *d, int block_id)
{
    return hashpipe_databuf_set_free((hashpipe_databuf_t *)d, block_id);
}

static inline int lwa_gpu_input_databuf_set_filled(lwa_gpu_input_databuf_t *d, int block_id)
{
    return hashpipe_databuf_set_filled((hashpipe_databuf_t *)d, block_id);
}

/*
 * OUTPUT BUFFER FUNCTIONS
 */

hashpipe_databuf_t *lwa_output_databuf_create(int instance_id, int databuf_id);

static inline void lwa_output_databuf_clear(lwa_output_databuf_t *d)
{
    hashpipe_databuf_clear((hashpipe_databuf_t *)d);
}

static inline lwa_output_databuf_t *lwa_output_databuf_attach(int instance_id, int databuf_id)
{
    return (lwa_output_databuf_t *)hashpipe_databuf_attach(instance_id, databuf_id);
}

static inline int lwa_output_databuf_detach(lwa_output_databuf_t *d)
{
    return hashpipe_databuf_detach((hashpipe_databuf_t *)d);
}

static inline int lwa_output_databuf_block_status(lwa_output_databuf_t *d, int block_id)
{
    return hashpipe_databuf_block_status((hashpipe_databuf_t *)d, block_id);
}

static inline int lwa_output_databuf_total_status(lwa_output_databuf_t *d)
{
    return hashpipe_databuf_total_status((hashpipe_databuf_t *)d);
}

static inline int lwa_output_databuf_wait_free(lwa_output_databuf_t *d, int block_id)
{
    return hashpipe_databuf_wait_free((hashpipe_databuf_t *)d, block_id);
}

static inline int lwa_output_databuf_busywait_free(lwa_output_databuf_t *d, int block_id)
{
    return hashpipe_databuf_busywait_free((hashpipe_databuf_t *)d, block_id);
}

static inline int lwa_output_databuf_wait_filled(lwa_output_databuf_t *d, int block_id)
{
    return hashpipe_databuf_wait_filled((hashpipe_databuf_t *)d, block_id);
}

static inline int lwa_output_databuf_busywait_filled(lwa_output_databuf_t *d, int block_id)
{
    return hashpipe_databuf_busywait_filled((hashpipe_databuf_t *)d, block_id);
}

static inline int lwa_output_databuf_set_free(lwa_output_databuf_t *d, int block_id)
{
    return hashpipe_databuf_set_free((hashpipe_databuf_t *)d, block_id);
}

static inline int lwa_output_databuf_set_filled(lwa_output_databuf_t *d, int block_id)
{
    return hashpipe_databuf_set_filled((hashpipe_databuf_t *)d, block_id);
}

#endif // _LWA_DATABUF_H
