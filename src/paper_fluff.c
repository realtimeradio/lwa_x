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

/*
 * The xGPU input buffer is time x channel x antenna x pol x complexity x 4[samples]. Each real/imag
 * input component is a 1Byte + 1Byte complex value.
 *
 * The netthread output buffer has a shape dictated by the structure of a single
 * packet. It has dimensions time x antenna x channel x times-per-packet x pol.
 * Each pol is a 4bit + 4bit complex value.
 * Thus, the fluffing thread has to do some gymnastics to write the GPU input buffer
 * correctly. You'll probably want this to hand -- https://software.intel.com/sites/landingpage/IntrinsicsGuide/ (!)
 * 
 * This code isn't generic anymore, it assumes that Nt = 2. Sometimes the variable Nt is used to provide
 * hints about what might need changing for different Nt values.
*/
int paper_fluff(const uint64_t const * const in, uint64_t * out)
{
  //uint64_t v0, v1, v2, v3;
  //vec_t in0, in1, in2, in3;
  //vec_t v0, v1, v2, v3;
  int m, c, t, a, i, j, k, n, p, wc, wa, wt;

  __m256i int0, int1, int2, int3, outx, outy;

  __m256i *p_in = (__m256i *)in;
  __m256i *p_out;
  const __m256i maskhi = _mm256_set_epi64x(0xff00ff00ff00ff00ULL, 0xff00ff00ff00ff00ULL,
                                           0xff00ff00ff00ff00ULL, 0xff00ff00ff00ff00ULL);
  const __m256i masklo = _mm256_set_epi64x(0x00ff00ff00ff00ffULL, 0x00ff00ff00ff00ffULL,
                                           0x00ff00ff00ff00ffULL, 0x00ff00ff00ff00ffULL);
  const __m256i shuffle_map = _mm256_set_epi64x(0x0f0b07030e0a0602ULL, 0x0d0905010c080400ULL,
                                                0x0f0b07030e0a0602ULL, 0x0d0905010c080400ULL);
  //__m256i row[16], temp0[16], temp1[16], temp2[16];
  __m256i row[16], temp[16];

  p_out = (__m256i *)out;
  // Loop over different time chunks
  for(m=0; m<Nm*Nt/4; m++) { //dealing with 4 time samples
    // Loop over channels
    for(c=0; c<Nc; c=c+8) { // dealing with 8 channels every step
      for(t=0; t<(4/Nt); t=t+Nt) { // Need to load blocks of 4 times to satisfy XGPU. Can only load Nt per row read
        for(a=0; a<Na; a=a+4) {  // dealing with 4 anetnnas per transpose
          // Load a 16-chan-sample by 4 antenna matrix into 8 x 256 bit registers
          for(i=0; i<8; i++) {
            //row[i]  = _mm256_stream_load_si256(p_in + paper_input_databuf_data_idx256(m,16*a+i,c,t));
            row[i]  = _mm256_load_si256(p_in + paper_input_databuf_data_idx256(m,a+i/Nt,c,t+(Nt*(i%Nt))));
            //row[i]  = _mm256_stream_load_si256(p_in + paper_input_databuf_data_idx256(m,16*a+i,c,t));
          }
          /* TRANSPOSE in blocks of 32 bits*/
          for(j=0; j<4; j++) {
            temp[j]   = _mm256_unpacklo_epi32(row[2*j], row[2*j+1]);
            temp[j+4] = _mm256_unpackhi_epi32(row[2*j], row[2*j+1]);
          }
          for(k=0; k<4; k++) {
            row[k]   = _mm256_unpacklo_epi64(temp[2*k], temp[2*k+1]);
            row[k+4] = _mm256_unpackhi_epi64(temp[2*k], temp[2*k+1]);
          }
          /* END TRANSPOSE -- actually we have transposed 4 4x4 blocks of the matrix.
             For the full transpose we'd need an unpackX_epi128, which isn't available.*/

          /* 
             Fluff each 256 bit row of the transposed matrix into
             two 256-bit fluffed words. The loop here is over rows, each of which is 4 times
             for a pair of antennas..
             The matrix here should be (for Nt=2):

             ANT[ 0: 1]C0T[0:3], ANT[ 0: 1]C4T[0:3]
             ANT[ 0: 1]C1T[0:3], ANT[ 0: 1]C5T[0:3]
             ANT[ 0: 1]C2T[0:3], ANT[ 0: 1]C6T[0:3]
             ANT[ 0: 1]C3T[0:3], ANT[ 0: 1]C7T[0:3]

             ANT[ 2: 3]C0T[0:3], ANT[ 2: 3]C4T[0:3]
             ANT[ 2: 3]C1T[0:3], ANT[ 2: 3]C5T[0:3]
             ANT[ 2: 3]C2T[0:3], ANT[ 2: 3]C6T[0:3]
             ANT[ 2: 3]C3T[0:3], ANT[ 2: 3]C7T[0:3]
          */

          /*
             We still need to swizzle the blocks of 32 bits.
             They are currently time x pol x complexity
             We require pol x complexity x time
             And we must fluff each real/imag component to 8 bits
          */

          // Convenience loops: Separate into the top and bottom
          // part of this matrix
          for(p=0; p<2; p++) {
            wa = p * 2; // antenna offset. I.e., antenna being processed is a+wa
            // Loop over the 4 rows in each of the top/bottom halves
            for(n=0; n<4; n++) {
              wc = n; // channel to write (the first one appearing in the row)
              wt = 0; // time sample to write
              // fluff and load
              // shift so that we have four bitvectors, where the top 4 bits
              // of each 16-bit word are consecutive values
              // eg, for n=0, p=0; int0 is:
              //    A0C0T0xr, A0C0T0xi, A0C0T0yr, A0C0T0yi // 16-bits
              //    A0C0T1xr, A0C0T1xi, A0C0T1yr, A0C0T1yi 
              //    A0C0T2xr, A0C0T2xi, A0C0T2yr, A0C0T2yi 
              //    A0C0T3xr, A0C0T3xi, A0C0T3yr, A0C0T3yi 
              //    A1C0T0xr, A1C0T0xi, A1C0T0yr, A1C0T0yi // 16-bits
              //    A1C0T1xr, A1C0T1xi, A1C0T1yr, A1C0T1yi 
              //    A1C0T2xr, A1C0T2xi, A1C0T2yr, A1C0T2yi 
              //    A1C0T3xr, A1C0T3xi, A1C0T3yr, A1C0T3yi 
              //    A0C4T0xr, A0C4T0xi, A0C4T0yr, A0C4T0yi // 16-bits
              //    A0C4T1xr, A0C4T1xi, A0C4T1yr, A0C4T1yi 
              //    A0C4T2xr, A0C4T2xi, A0C4T2yr, A0C4T2yi 
              //    A0C4T3xr, A0C4T3xi, A0C4T3yr, A0C4T3yi 
              //    A1C4T0xr, A1C4T0xi, A1C4T0yr, A1C4T0yi // 16-bits
              //    A1C4T1xr, A1C4T1xi, A1C4T1yr, A1C4T1yi 
              //    A1C4T2xr, A1C4T2xi, A1C4T2yr, A1C4T2yi 
              //    A1C4T3xr, A1C4T3xi, A1C4T3yr, A1C4T3yi 
              int0 = row[4*p + n];
              // eg, for n=0, p=0; int1 is:
              //    A0C0T0xi, A0C0T0yr, A0C0T0yi, 00000000 // 16-bits
              //    A0C0T1xi, A0C0T1yr, A0C0T1yi, 00000000 
              //    A0C0T2xi, A0C0T2yr, A0C0T2yi, 00000000 
              //    A0C0T3xi, A0C0T3yr, A0C0T3yi, 00000000 
              //    A1C0T0xi, A1C0T0yr, A1C0T0yi, 00000000 // 16-bits
              //    A1C0T1xi, A1C0T1yr, A1C0T1yi, 00000000 
              //    A1C0T2xi, A1C0T2yr, A1C0T2yi, 00000000 
              //    A1C0T3xi, A1C0T3yr, A1C0T3yi, 00000000 
              //    A0C4T0xi, A0C4T0yr, A0C4T0yi, 00000000 // 16-bits
              //    A0C4T1xi, A0C4T1yr, A0C4T1yi, 00000000 
              //    A0C4T2xi, A0C4T2yr, A0C4T2yi, 00000000 
              //    A0C4T3xi, A0C4T3yr, A0C4T3yi, 00000000 
              //    A1C4T0xi, A1C4T0yr, A1C4T0yi, 00000000 // 16-bits
              //    A1C4T1xi, A1C4T1yr, A1C4T1yi, 00000000 
              //    A1C4T2xi, A1C4T2yr, A1C4T2yi, 00000000 
              //    A1C4T3xi, A1C4T3yr, A1C4T3yi, 00000000 
              int1 = _mm256_slli_epi16(row[4*p + n], 4);
              // eg, for n=0, p=0; int2 is:
              //    A0C0T0yr, A0C0T0yi, 00000000, 00000000 // 16-bits
              //    A0C0T1yr, A0C0T1yi, 00000000, 00000000 
              //    A0C0T2yr, A0C0T2yi, 00000000, 00000000 
              //    A0C0T3yr, A0C0T3yi, 00000000, 00000000 
              //    A1C0T0yr, A1C0T0yi, 00000000, 00000000 // 16-bits
              //    A1C0T1yr, A1C0T1yi, 00000000, 00000000 
              //    A1C0T2yr, A1C0T2yi, 00000000, 00000000 
              //    A1C0T3yr, A1C0T3yi, 00000000, 00000000 
              //    A0C4T0yr, A0C4T0yi, 00000000, 00000000 // 16-bits
              //    A0C4T1yr, A0C4T1yi, 00000000, 00000000 
              //    A0C4T2yr, A0C4T2yi, 00000000, 00000000 
              //    A0C4T3yr, A0C4T3yi, 00000000, 00000000 
              //    A1C4T0yr, A1C4T0yi, 00000000, 00000000 // 16-bits
              //    A1C4T1yr, A1C4T1yi, 00000000, 00000000 
              //    A1C4T2yr, A1C4T2yi, 00000000, 00000000 
              //    A1C4T3yr, A1C4T3yi, 00000000, 00000000 
              int2 = _mm256_slli_epi16(row[4*p + n], 8);
              // eg, for n=0, p=0; int3 is:
              //    A0C0T0yi, 00000000, 00000000, 00000000 // 16-bits
              //    A0C0T1yi, 00000000, 00000000, 00000000 
              //    A0C0T2yi, 00000000, 00000000, 00000000 
              //    A0C0T3yi, 00000000, 00000000, 00000000 
              //    A1C0T0yi, 00000000, 00000000, 00000000 // 16-bits
              //    A1C0T1yi, 00000000, 00000000, 00000000 
              //    A1C0T2yi, 00000000, 00000000, 00000000 
              //    A1C0T3yi, 00000000, 00000000, 00000000 
              //    A0C4T0yi, 00000000, 00000000, 00000000 // 16-bits
              //    A0C4T1yi, 00000000, 00000000, 00000000 
              //    A0C4T2yi, 00000000, 00000000, 00000000 
              //    A0C4T3yi, 00000000, 00000000, 00000000 
              //    A1C4T0yi, 00000000, 00000000, 00000000 // 16-bits
              //    A1C4T1yi, 00000000, 00000000, 00000000 
              //    A1C4T2yi, 00000000, 00000000, 00000000 
              //    A1C4T3yi, 00000000, 00000000, 00000000 
              int3 = _mm256_slli_epi16(row[4*p + n], 12);

              // shift down with sign extension. AKA, "fluff"
              // This makes the top 8 bits of the int0 and int2 values valid 8-bit signed numbers
              // The shift by 12 makes the lower 8-bits of the int1 and in3 values 8-bit signed numbers
              // eg, for n=0, p=0; int0 is: A0C0T0xr-8b, A0C0T0xi, A0C0T0yr ...
              int0 = _mm256_srai_epi16(int0, 4);
              // eg, for n=0, p=0; int1 is: A0C0T0xi-16b ...
              int1 = _mm256_srai_epi16(int1, 12);
              // eg, for n=0, p=0; int2 is: A0C0T0yr-8b, A0C0T0yi, 00000000 ...
              int2 = _mm256_srai_epi16(int2, 4);
              // eg, for n=0, p=0; int3 is: A0C0T0yi-16b ...
              int3 = _mm256_srai_epi16(int3, 12);

              // Mask bits, ready to AND together
              // eg, for n=0, p=0; int0 is: A0C0T0xr-8b, 8'b0 ...
              int0 = _mm256_and_si256(int0, maskhi);
              // eg, for n=0, p=0; int1 is: 8'b0, A0C0T0xi-8b ...
              int1 = _mm256_and_si256(int1, masklo);
              // eg, for n=0, p=0; int2 is: A0C0T0yr-8b, 8'b0 ...
              int2 = _mm256_and_si256(int2, maskhi);
              // eg, for n=0, p=0; int3 is: 8'b0, A0C0T0yi-8b ...
              int3 = _mm256_and_si256(int3, masklo);

              // AND together pairs of vectors
              // eg, for n=0, p=0; outlo is: A0C0T0xr-8b, A0C0T0xi-8b ...
              //    A0C0T0xr, A0C0T0xi // 16-bits
              //    A0C0T1xr, A0C0T1xi 
              //    A0C0T2xr, A0C0T2xi 
              //    A0C0T3xr, A0C0T3xi 
              //    A1C0T0xr, A1C0T0xi // 16-bits
              //    A1C0T1xr, A1C0T1xi 
              //    A1C0T2xr, A1C0T2xi 
              //    A1C0T3xr, A1C0T3xi 
              //    A0C4T0xr, A0C4T0xi // 16-bits
              //    A0C4T1xr, A0C4T1xi 
              //    A0C4T2xr, A0C4T2xi 
              //    A0C4T3xr, A0C4T3xi 
              //    A1C4T0xr, A1C4T0xi // 16-bits
              //    A1C4T1xr, A1C4T1xi 
              //    A1C4T2xr, A1C4T2xi 
              //    A1C4T3xr, A1C4T3xi 
              outx = _mm256_and_si256(int0, int1);
              // eg, for n=0, p=0; outlo is: A0C0T0yr-8b, A0C0T0yi-8b
              //    A0C0T0yr, A0C0T0yi // 16-bits
              //    A0C0T1yr, A0C0T1yi 
              //    A0C0T2yr, A0C0T2yi 
              //    A0C0T3yr, A0C0T3yi 
              //    A1C0T0yr, A1C0T0yi // 16-bits
              //    A1C0T1yr, A1C0T1yi 
              //    A1C0T2yr, A1C0T2yi 
              //    A1C0T3yr, A1C0T3yi 
              //    A0C4T0yr, A0C4T0yi // 16-bits
              //    A0C4T1yr, A0C4T1yi 
              //    A0C4T2yr, A0C4T2yi 
              //    A0C4T3yr, A0C4T3yi 
              //    A1C4T0yr, A1C4T0yi // 16-bits
              //    A1C4T1yr, A1C4T1yi 
              //    A1C4T2yr, A1C4T2yi 
              //    A1C4T3yr, A1C4T3yi 
              outy = _mm256_and_si256(int2, int3);

              // For a full 256-bit word (for n=0, p=0):
              // Interleave these to gather X and Y pols
              // eg, for n=0, p=0; int0 = 
              //    A0C0T0 [x/y][r/i] // 32-bits
              //    A0C0T1 [x/y][r/i]
              //    A0C0T2 [x/y][r/i]
              //    A0C0T3 [x/y][r/i]
              //    A1C0T0 [x/y][r/i]
              //    A1C0T1 [x/y][r/i]
              //    A1C0T2 [x/y][r/i]
              //    A1C0T3 [x/y][r/i]
              int0 = _mm256_unpacklo_epi16(outx, outy);
              // eg, for n=0, p=0; int1 =
              //    A0C4T0 [x/y][r/i] // 32-bits
              //    A0C4T1 [x/y][r/i]
              //    A0C4T2 [x/y][r/i]
              //    A0C4T3 [x/y][r/i]
              //    A1C4T0 [x/y][r/i]
              //    A1C4T1 [x/y][r/i]
              //    A1C4T2 [x/y][r/i]
              //    A1C4T3 [x/y][r/i]
              int1 = _mm256_unpackhi_epi16(outx, outy);

              // Finally perform a bytewise transpose using a shuffle
              // Outputs are now ordered pol x complexity x time ([ants 0..1][real/imag][time 0..3])
              int0 = _mm256_shuffle_epi8(int0, shuffle_map); // Channel 0
              int1 = _mm256_shuffle_epi8(int1, shuffle_map); // Channel 4

              // Store these two 256-bit words with _mm256_stream_si256
              //fprintf(stdout, "m:%d, a:%d, c:%d, t:%d, word1: %d, word2: %d\n", m,a,c,t,paper_gpu_input_databuf_data_idx(m,16*a+8*p,t+n,c), paper_gpu_input_databuf_data_idx(m,16*a+8*p,t+n+8,c));
              _mm256_store_si256(p_out + paper_gpu_input_databuf_data_idx256(m, a+wa, wt, c + wc), int0);
              _mm256_store_si256(p_out + paper_gpu_input_databuf_data_idx256(m, a+wa, wt, c + wc + 4), int1);
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
