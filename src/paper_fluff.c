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
#define OUTPUT_STRIDE ((2*Na)/sizeof(__m256i))

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
  int m, a, w;

  __m128i lo128, hi128;
  __m256i out256 = _mm256_setzero_si256();

  __m128i *p_in = (__m128i *)in;
  __m256i *p_out;
  const __m128i mask = _mm_set_epi64((__m64)0xf0f0f0f0f0f0f0f0ULL, (__m64)0xf0f0f0f0f0f0f0f0ULL);

  for(m=0; m<Nm; m++) {
    //for(f=0; f<Nf; f++) {
    for(a=0; a<Na; a++) {
        // Each F engine provides 32 samples per channel per time, which fluffs
        // to 64 bytes, which equals two 256 bit words.
        // p_out = out + mcnt offset + feng offset
        p_out = (__m256i *)out + m*Nt*Nc*Na*2 + a*2;
        for(w=0; w<N_WORD128_PER_PACKET/2; w++) {

          // Load 128 bits (16 bytes) with _mm_stream_load_si128
          lo128 = _mm_stream_load_si128(p_in++);
          // Shift lo left by 4 bits to make hi
          hi128 = _mm_slli_epi64(lo128, 4);
          // Mask off unwanted bits
          lo128 = _mm_and_si128(lo128, mask);
          hi128 = _mm_and_si128(hi128, mask);
          // Transfer lo128 to lower half of __m256i variable 'out256'
          // Use pragmas to avoid "uninitialized use of out256" warning/error
          out256 = _mm256_insertf128_si256(out256, lo128, 0);
          // Transfer hi128 to upper half of __m256i variable 'out256'
          out256 = _mm256_insertf128_si256(out256, hi128, 1);
          // Store 256 bits (32 bytes) with _mm256_stream_si256
          _mm256_stream_si256(p_out++, out256);

          // Now do it again...

          // Load 128 bits (16 bytes) with _mm_stream_load_si128
          lo128 = _mm_stream_load_si128(p_in++);
          // Shift lo left by 4 bits to make hi
          hi128 = _mm_slli_epi64(lo128, 4);
          // Mask off unwanted bits
          lo128 = _mm_and_si128(lo128, mask);
          hi128 = _mm_and_si128(hi128, mask);
          // Transfer lo128 to lower half of __m256i variable 'out256'
          // Use pragmas to avoid "uninitialized use of out256" warning/error
          out256 = _mm256_insertf128_si256(out256, lo128, 0);
          // Transfer hi128 to upper half of __m256i variable 'out256'
          out256 = _mm256_insertf128_si256(out256, hi128, 1);
          // Store 256 bits (32 bytes) with _mm256_stream_si256
          _mm256_stream_si256(p_out, out256);

          p_out += OUTPUT_STRIDE-1;

        }
    }
  }

  // Return number of 64 bit words fluffed
  return Nm*Na*2*N_WORD128_PER_PACKET;
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
