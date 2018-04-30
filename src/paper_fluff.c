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

/*
 * The xGPU input buffer is time x channel x antenna x pol. Each pol input
 * is a 1Byte + 1Byte complex value.
 * So a 256 bit (32 Byte) chunk of memory should contain a single time/channel,
 * and 16 inputs (8 antennas x dual-pol).
 *
 * The netthread output buffer has a shape dictated by the structure of a single
 * packet. It has dimensions time x antenna x channel x times-per-packet x pol.
 * Each pol is a 4bit + 4bit complex value.
 * Thus, the fluffing thread has to do some gymnastics to write the GPU input buffer
 * more than one antenna at a time.
*/
int paper_fluff(const uint64_t const * const in, uint64_t * out)
{
  //uint64_t v0, v1, v2, v3;
  //vec_t in0, in1, in2, in3;
  //vec_t v0, v1, v2, v3;
  int m, ct, c, t, a, i, j, k, l, n, p, wc, wa, wt;

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
  // Loop over different time chunks
  for(m=0; m<Nm; m++) {
    // Loop over X-engine channel-times
    for(ct=0; ct<(Nc*Nt); ct=ct+16) {
      c = ct / Nt;
      t = 0; //Only valid if Nt < 16
      for(a=0; a<Na; a=a+16) {
        // Load a 16-chan-sample by 16-antenna matrix into 16 x 256 bit registers
        for(i=0; i<16; i++) {
          //row[i]  = _mm256_stream_load_si256(p_in + paper_input_databuf_data_idx256(m,16*a+i,c,t));
          row[i]  = _mm256_load_si256(p_in + paper_input_databuf_data_idx256(m,a+i,c,t));
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
           For the full transpose we'd need a unpackX_epi128, which isn't available.*/
        /* Fluff each 256 bit row of the transposed matrix into
           two 256-bit fluffed words. The loop here is over rows, which are times.
           The matrix here should be:
           ANT[ 0: 7]C0T0, ANT[ 0: 7]C4T0
           ANT[ 0: 7]C0T1, ANT[ 0: 7]C4T1
           ANT[ 0: 7]C1T0, ANT[ 0: 7]C5T0
           ANT[ 0: 7]C1T1, ANT[ 0: 7]C5T1
           ANT[ 0: 7]C2T0, ANT[ 0: 7]C6T0
           ANT[ 0: 7]C2T1, ANT[ 0: 7]C6T1
           ANT[ 0: 7]C3T0, ANT[ 0: 7]C7T0
           ANT[ 0: 7]C3T1, ANT[ 0: 7]C7T1
           ANT[ 8:15]C0T0, ANT[ 8:15]C4T0
           ANT[ 8:15]C0T1, ANT[ 8:15]C4T1
           ANT[ 8:15]C1T0, ANT[ 8:15]C5T0
           ANT[ 8:15]C1T1, ANT[ 8:15]C5T1
           ANT[ 8:15]C2T0, ANT[ 8:15]C6T0
           ANT[ 8:15]C2T1, ANT[ 8:15]C6T1
           ANT[ 8:15]C3T0, ANT[ 8:15]C7T0
           ANT[ 8:15]C3T1, ANT[ 8:15]C7T1 */

        // Convenience loops: Separate into the top and bottom
        // part of this matrix
        for(p=0; p<2; p++) {
          wa = p * 8; // antenna to write
          for(n=0; n<8; n++) {
            wc = n / Nt; // channel to write (the first one appearing in the row)
            wt = n % Nt; // time sample to write
            // fluff and load
            // shift so that we have four bitvectors, where the top 4 bits
            // of each 16-bit word are consecutive values
            // eg, for n=0, p=0; int0 is: A0C0T0xr, A0C0T0xi, A0C0T0yr, A0C0T0yi 
            int0 = temp[8*p + n];
            // eg, for n=0, p=0; int1 is: A0C0T0xi, A0C0T0yr, A0C0T0yi, 00000000
            int1 = _mm256_slli_epi16(temp[8*p + n], 4);
            // eg, for n=0, p=0; int2 is: A0C0T0yr, A0C0T0yi, 00000000, 00000000
            int2 = _mm256_slli_epi16(temp[8*p + n], 8);
            // eg, for n=0, p=0; int3 is: A0C0T0yi, 00000000, 00000000, 00000000
            int3 = _mm256_slli_epi16(temp[8*p + n], 12);

            // shift down with sign extension. AKA, "fluff"
            // This makes the top 8 bits of the int0 and int2 values valid 8-bit signed numbers
            // The shift by 12 makes the lower 8-bits of the int1 and in3 values 8-bit signed numbers
            // eg, for n=0, p=0; int0 is: A0C0T0xr-8b, A0C0T0xi, A0C0T0yr
            int0 = _mm256_srai_epi16(int0, 4);
            // eg, for n=0, p=0; int1 is: A0C0T0xi-16b
            int1 = _mm256_srai_epi16(int1, 12);
            // eg, for n=0, p=0; int2 is: A0C0T0yr-8b, A0C0T0yi, 00000000
            int2 = _mm256_srai_epi16(int2, 4);
            // eg, for n=0, p=0; int3 is: A0C0T0yi-16b
            int3 = _mm256_srai_epi16(int3, 12);

            // Mask bits, ready to AND together
            // eg, for n=0, p=0; int0 is: A0C0T0xr-8b, 8'b0
            int0 = _mm256_and_si256(int0, maskhi);
            // eg, for n=0, p=0; int1 is: 8'b0, A0C0T0xi-8b
            int1 = _mm256_and_si256(int1, masklo);
            // eg, for n=0, p=0; int2 is: A0C0T0yr-8b, 8'b0
            int2 = _mm256_and_si256(int2, maskhi);
            // eg, for n=0, p=0; int3 is: 8'b0, A0C0T0yi-8b
            int3 = _mm256_and_si256(int3, masklo);

            // AND together pairs of vectors
            // eg, for n=0, p=0; outlo is: A0C0T0xr-8b, A0C0T0xi-8b
            outlo = _mm256_and_si256(int0, int1);
            // eg, for n=0, p=0; outlo is: A0C0T0yr-8b, A0C0T0yi-8b
            outhi = _mm256_and_si256(int2, int3);

            // For a full 256-bit word (for n=0, p=0):
            // outlo: A0C0T0x(r/i), A1C0T0x(r/i), ..., A7C0T0x(r/i), A0C4T0x(r/i), A7C4T0x(r/i)
            // outhi: A0C0T0y(r/i), A1C0T0y(r/i), ..., A7C0T0y(r/i), A0C4T0y(r/i), A7C4T0y(r/i)

            // Interleave these to gather X and Y pols
            // eg, for n=0, p=0; int0 = A0C0T0x(r/i), A0C0T0y(r/i), ..., A7C0T0x(r/i),  A7C0T0y(r/i)
            int0 = _mm256_unpacklo_epi16(outlo, outhi);
            // eg, for n=0, p=0; int1 = A0C4T0x(r/i), A0C4T0y(r/i), ..., A7C4T0x(r/i),  A7C4T0y(r/i)
            int1 = _mm256_unpackhi_epi16(outlo, outhi);

            // Store these two 256-bit words with _mm256_stream_si256
            //fprintf(stdout, "m:%d, a:%d, c:%d, t:%d, word1: %d, word2: %d\n", m,a,c,t,paper_gpu_input_databuf_data_idx(m,16*a+8*p,t+n,c), paper_gpu_input_databuf_data_idx(m,16*a+8*p,t+n+8,c));
            _mm256_store_si256(p_out + paper_gpu_input_databuf_data_idx256(m, a+wa, wt, wc), int0);
            _mm256_store_si256(p_out + paper_gpu_input_databuf_data_idx256(m, a+wa, wt, wc + (8/Nt)), int1);
            //_mm256_stream_si256(p_out + paper_gpu_input_databuf_data_idx256(m,16*a+8*p,t+n,c), outlo);
            //_mm256_stream_si256(p_out + paper_gpu_input_databuf_data_idx256(m,16*a+8*p,t+n+8,c), outhi);
            //_mm256_stream_si256(foo256, outlo);
            //_mm256_stream_si256(bar256, outhi);
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
