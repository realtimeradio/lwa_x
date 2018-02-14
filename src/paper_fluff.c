#include <stdio.h>
#include <pthread.h>
//#include <stdlib.h>
//#include <unistd.h>
//#include <stdint.h>
//#include <time.h>

#include <smmintrin.h>
#include <immintrin.h>

#include "paper_fluff.h"

#define N_WORD128_PER_PACKET  (N_BYTES_PER_PACKET/sizeof(__m128i))
// OUTPUT_STRIDE in is units of __m256i
//#define OUTPUT_STRIDE ((2*Nf*N_INPUTS_PER_PACKET)/sizeof(__m256i))
//#define OUTPUT_STRIDE ((2*Na)/sizeof(__m256i))
#define OUTPUT_STRIDE 1

#define paper_input_databuf_data_idx256(m,a,c,t) \
  ((((m * Na + a) * Nc + c) * t) * N_INPUTS_PER_PACKET / sizeof(__m256i))
#define paper_gpu_input_databuf_data_idx256(m,a,t,c) \
  (4*((m*Nt*Nc*Na) + (t*Nc*Na) + (c*Na) + a) / sizeof(__m256i))

#if 0
typedef union {
  __m128i m128;
  __m64   m64[2];
} vec_t;
#define m128(v) (v.m128)
#define m64(v,i) (v.m64[i])
#endif // 0

int paper_fluff_diag(const uint64_t const * const in, uint64_t * out)
{
  const __m128i mask = _mm_set_epi64((__m64)0xf0f0f0f0f0f0f0f0ULL, (__m64)0xf0f0f0f0f0f0f0f0ULL);

  __m128i * p128_in  = (__m128i *)in;
  __m256i * p256_out = (__m256i *)out;

  // Load 128 bits (16 bytes) with _mm_stream_load_si128
  __m128i lo = _mm_stream_load_si128(p128_in);

  // Shift lo left by 4 bits to make hi
  __m128i hi = _mm_slli_epi64(lo, 4);

  // Mask off unwanted bits
  lo = _mm_and_si128(lo, mask);
  hi = _mm_and_si128(hi, mask);

  // Transfer lo to lower half of __m256i variable 'out256'
  // Use pragmas to avoid "uninitialized use of out256" warning/error
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wuninitialized"
  __m256i out256 = _mm256_insertf128_si256(out256, lo, 0);
#pragma GCC diagnostic pop

  // Transfer hi to upper half of __m256i variable 'out256'
  out256 = _mm256_insertf128_si256(out256, hi, 1);

  // Store 256 bits (32 bytes) with _mm256_stream_si256
  _mm256_stream_si256(p256_out, out256);

  return 0;
}


int paper_fluff(const uint64_t const * const in, uint64_t * out)
{
  //uint64_t v0, v1, v2, v3;
  //vec_t in0, in1, in2, in3;
  //vec_t v0, v1, v2, v3;
  int m, c, t, a, i, j, k, l, n, p;

  __m256i int0, int1, outlo, int2, int3, outhi;

  __m256i *p_in = (__m256i *)in;
  __m256i *p_out;
  const __m256i maskhi = _mm256_set_epi64x(0xff00ff00ff00ff00ULL, 0xff00ff00ff00ff00ULL,
                                           0xff00ff00ff00ff00ULL, 0xff00ff00ff00ff00ULL);
  const __m256i masklo = _mm256_set_epi64x(0x00ff00ff00ff00ffULL, 0x00ff00ff00ff00ffULL,
                                           0x00ff00ff00ff00ffULL, 0x00ff00ff00ff00ffULL);
  //__m256i row[16], temp0[16], temp1[16], temp2[16];
  __m256i row[16], temp[16];

  p_out = (__m256i *)out;
  // Loop over different packets
  for(m=0; m<Nm; m++) {
    // Loop over X-engine channels
    for(c=0; c<Nc; c++) {
      // Loop over 16-sample time chunks
      for(t=0; t<Nt/16; t++) {
        // Loop over 16-antenna chunks
        for(a=0; a<Na/16; a++) {
          // Load a 16-sample by 16-antenna matrix into 16 x 128 bit registers
          for(i=0; i<16; i++) {
            //row[i]  = _mm256_stream_load_si256(p_in + paper_input_databuf_data_idx256(m,16*a+i,c,t));
            row[i]  = _mm256_load_si256(p_in + paper_input_databuf_data_idx256(m,16*a+i,c,t));
            //row[i]  = _mm256_stream_load_si256(p_in + paper_input_databuf_data_idx256(m,16*a+i,c,t));
          }
          /* TRANSPOSE */
          for(j=0; j<8; j++) {
            temp[j]   = _mm256_unpacklo_epi16(row[2*j], row[2*j+1]);
            temp[j+8] = _mm256_unpackhi_epi16(row[2*j], row[2*j+1]);
          }
          for(k=0; k<8; k++) {
            row[k]   = _mm256_unpacklo_epi32(temp[2*k], temp[2*k+1]);
            row[k+8] = _mm256_unpackhi_epi32(temp[2*k], temp[2*k+1]);
          }
          for(l=0; l<8; l++) {
            temp[l]   = _mm256_unpacklo_epi64(row[2*l], row[2*l+1]);
            temp[l+8] = _mm256_unpackhi_epi64(row[2*l], row[2*l+1]);
          }
          /* END TRANSPOSE -- actually we have transposed 4 8x8 blocks of the matrix.
             For the full transpose we'd need a unpackX_epi128, which doesn't exist.*/
          /* Fluff each 256 bit row of the transposed matrix into
             two 256-bit fluffed words. The loop here is over rows, which are times.
             The matrix here should be:
             ANT[ 0: 7]T0, ANT[ 0: 7]T8
             ANT[ 0: 7]T1, ANT[ 0: 7]T9
             ANT[ 0: 7]T2, ANT[ 0: 7]T10
             ANT[ 0: 7]T3, ANT[ 0: 7]T11
             ANT[ 0: 7]T4, ANT[ 0: 7]T12
             ANT[ 0: 7]T5, ANT[ 0: 7]T13
             ANT[ 0: 7]T6, ANT[ 0: 7]T14
             ANT[ 0: 7]T7, ANT[ 0: 7]T15
             ANT[ 8:15]T0, ANT[ 8:15]T8
             ANT[ 8:15]T1, ANT[ 8:15]T9
             ANT[ 8:15]T2, ANT[ 8:15]T10
             ANT[ 8:15]T3, ANT[ 8:15]T11
             ANT[ 8:15]T4, ANT[ 8:15]T12
             ANT[ 8:15]T5, ANT[ 8:15]T13
             ANT[ 8:15]T6, ANT[ 8:15]T14
             ANT[ 8:15]T7, ANT[ 8:15]T15 */
          for(p=0; p<2; p++) {
            for(n=0; n<8; n++) {
              // fluff and load
              // shift so that we have four bitvectors, where the top 4 bits
              // of each 16-bit word are consecutive values
              int0 = temp[8*p + n];
              int1 = _mm256_slli_epi16(temp[8*p + n], 4);
              int2 = _mm256_slli_epi16(temp[8*p + n], 8);
              int3 = _mm256_slli_epi16(temp[8*p + n], 12);

              // shift down with sign extension. AKA, "fluff"
              // This makes the top 8 bits of the int0 and int2 values valid 8-bit signed numbers
              // The shift by 12 makes the lower 8-bits of the int1 and in3 values 8-bit signed numbers
              int0 = _mm256_srai_epi16(int0, 4);
              int1 = _mm256_srai_epi16(int1, 12);
              int2 = _mm256_srai_epi16(int2, 4);
              int3 = _mm256_srai_epi16(int3, 12);

              // Mask bits, ready to AND together
              int0 = _mm256_and_si256(int0, maskhi);
              int1 = _mm256_and_si256(int1, masklo);
              int2 = _mm256_and_si256(int2, maskhi);
              int3 = _mm256_and_si256(int3, masklo);

              // AND together pairs of vectors
              outlo = _mm256_and_si256(int0, int1);
              outhi = _mm256_and_si256(int2, int3);

              // Store two 256-bit words with _mm256_stream_si256
              //fprintf(stdout, "m:%d, a:%d, c:%d, t:%d, word1: %d, word2: %d\n", m,a,c,t,paper_gpu_input_databuf_data_idx(m,16*a+8*p,t+n,c), paper_gpu_input_databuf_data_idx(m,16*a+8*p,t+n+8,c));
              _mm256_store_si256(p_out + paper_gpu_input_databuf_data_idx256(m,16*a+8*p,t+n,c), outlo);
              _mm256_store_si256(p_out + paper_gpu_input_databuf_data_idx256(m,16*a+8*p,t+n+8,c), outhi);
              //_mm256_stream_si256(p_out + paper_gpu_input_databuf_data_idx256(m,16*a+8*p,t+n,c), outlo);
              //_mm256_stream_si256(p_out + paper_gpu_input_databuf_data_idx256(m,16*a+8*p,t+n+8,c), outhi);
              //_mm256_stream_si256(foo256, outlo);
              //_mm256_stream_si256(bar256, outhi);
            
            }
          }
        }
      }
    }
  }
  // Return number of 64 bit words fluffed
  //return Nm*Na*2*N_WORD128_PER_PACKET;
  return N_BYTES_PER_BLOCK >> 3;
}


int paper_fluff_lut(const uint64_t const * const in, uint64_t * out)
{
  uint16_t * in16  = (uint16_t *) in;
  uint32_t * out32 = (uint32_t *) out;
  int i;

  for(i=0; i<N_BYTES_PER_BLOCK>>1; i++) {
    out32[i] = lut[in16[i]];
  }

  //uint64_t * in64  = (uint64_t *) in;
  //uint64_t * out64 = (uint64_t *) out;
  //int i;

  //for(i=0; i<N_BYTES_PER_BLOCK>>3; i++) {
  //  out64[2*i]     = (lut[(in64[i]>>48)] << 32) +        (lut[(in64[i]>>32) & 0xff]);
  //  out64[2*i + 1] = (lut[(in64[i]>>16) & 0xff] << 32) + (lut[(in64[i] & 0xff)]);
  //}

  // Return number of 64 bit words fluffed
  return N_BYTES_PER_BLOCK >> 3;
}
