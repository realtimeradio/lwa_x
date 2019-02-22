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

typedef struct BDA_params {
  int intbin_cntr[5];  // No. of baselines in each bin
  int *ant_pair_0[5];  // Ant0 of pair
  int *ant_pair_1[5];  // Ant1 of pair
} BDA_params_t; 

BDA_params_t bdaparams;
static XGPUInfo xgpu_info;

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
static off_t casper_index(const int in0, const int in1, const int n)
{
  const int a0 = in0 >> 1;
  const int a1 = in1 >> 1;
  const int p0 = in0 & 1;
  const int p1 = in1 & 1;
  const int delta = a1-a0;
  const int num_words_per_cell = 8;
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
  //printf("%s: a0=%d, a1=%d, delta=%d, cell_index=%d\n", __FUNCTION__, a0, a1, delta, cell_index);
  // Pol offset
  const int pol_offset = 2*(2*(p0^p1) + p0);
  // Word index (in units of words (i.e. floats) of real component
  const int index = (cell_index * num_words_per_cell) + pol_offset;
  return index;
}

// For each channel, a casper ordered buffer contains four complex values for
// each pair of input pairs.  Thus, the number of complex values in a casper
// ordered buffer are: 4 * (N/2 * (N/2 + 1)) / 2 = N * (N/2 + 1)
#define N_CASPER_COMPLEX_PER_CHAN (N_INPUTS * (N_INPUTS/2 + 1))

// Lookup table mapping casper_idx to regtile_idx
static off_t *idx_map;

static int init_idx_map()
{
  int a0, a1, p0, p1, i, j;
  idx_map = malloc(N_CASPER_COMPLEX_PER_CHAN * sizeof(off_t));
  if(!idx_map) {
    return -1;
  }

  for(a1=0; a1<N_INPUTS/2; a1++) {
    for(a0=0; a0<=a1; a0++) {
      for(p0=0; p0<2; p0++) {
        for(p1=0; p1<2; p1++) {
          i = 2*a0 + p0;
          j = 2*a1 + p1;
          idx_map[casper_index(i,j,N_INPUTS)/2] = regtile_index(i,j);
        }
      }
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

    // Load baseline averaging params from config file
    FILE *fp;  
    int i,j, bin, a0, a1, inttime;
    int cntr[5] = {0,0,0,0,0};
    
    if((fp=fopen("bda_config.txt","r")) == NULL){
       printf("Cannot open the configuration file.\n");
       exit(1);
    }
    while(fscanf(fp, "%d %d %d", &a0, &a1, &inttime)!=EOF){
       if((inttime == 0) || !CHECK_PWR2(inttime)){
         printf("(%d,%d): Samples to integrate not power of 2!\n",a0,a1);
         exit(1);
       }
       bda_params.intbin_cntr[LOG(inttime)]++;
    }

    // malloc for storing ant pairs
    for(j=0; j<5; j++){
      bda_params.ant_pair_0[j] = (int *)malloc(bda_params.intbin_cntr[j]*sizeof(int));
      bda_params.ant_pair_1[j] = (int *)malloc(bda_params.intbin_cntr[j]*sizeof(int));
    }
    rewind(fp);
    while(fscanf(fp, "%d %d %d", &a0, &a1, &inttime)!=EOF){
      bin = LOG(inttime);
      bda_params.ant_pair_0[bin][ctr[bin]] = a0;
      bda_params.ant_pair_1[bin][ctr[bin]++] = a1; 
    }
    fclose(fp);

    // Malloc for integration buffers
     

    // Success!
    return 0;
}

#define ELAPSED_NS(start,stop) \
  (((int64_t)stop.tv_sec-start.tv_sec)*1000*1000*1000+(stop.tv_nsec-start.tv_nsec))

