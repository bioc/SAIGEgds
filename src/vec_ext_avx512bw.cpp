// ===============================================================
//
// vec_ext_avx512bw.cpp: optimization with AVX512F+AVX512BW vectorization
//
// Copyright (C) 2022   Xiuwen Zheng
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include "vec_ext.h"
#include "vectorization.h"

#ifdef VEC_CPU_ARCH_X86_AVX512BW
extern const bool VEC_ALGORITHM_AVX512BW = true;
#else
extern const bool VEC_ALGORITHM_AVX512BW = false;
#endif


#define SIMD_NAME(NAME)  NAME ## _avx512bw
#define THROW_ERROR      throw "No AVX512BW support!"


#ifdef VEC_CPU_ARCH_X86_AVX512BW

#ifdef __ICC
	#pragma intel optimization_parameter target_arch=CORE-AVX512
#   define TARGET_AVX512    __attribute__((target("avx512f")))
#else
#   if !defined(__AVX512F__) && !defined(__clang__)
		#pragma GCC target("avx512f")
#   endif
#   if !defined(__AVX512BW__) && !defined(__clang__)
		#pragma GCC target("avx512bw")
#   endif
#   define TARGET_AVX512    __attribute__((target("avx512f,avx512bw")))
#endif

#include <xmmintrin.h>  // SSE
#include <emmintrin.h>  // SSE2
#include <immintrin.h>  // AVX, AVX2, AVX512F


#undef SIMD_NAME
#define SIMD_NAME(NAME)    MATH_OFAST TARGET_AVX512 NAME ## _avx512bw


extern "C" void SIMD_NAME(grm_calc_update_f32)(
	const uint8_t p_i[], const uint8_t p_j[], size_t n, const float p_G[],
	int &miss_n, double &sum)
{
	// constants
	const __m512i shr_n4 = _mm512_set_epi32(
		30, 28, 26, 24,  22, 20, 18, 16,
		14, 12, 10, 8,   6, 4, 2, 0);
	const __m512i shr_n8 = _mm512_set_epi64(14, 12, 10, 8, 6, 4, 2, 0);
	const __m512i X03 = _mm512_set1_epi8(0x03);
	const __m512i XFF = _mm512_set1_epi32(0xFF);
	const __m512i T16 = _mm512_set_epi64(
		// look-up table instead of bitwise in vec_ext_def.cpp
		0x0707070707060504LL, 0x0705030107040100LL,
		0x0707070707060504LL, 0x0705030107040100LL,
		0x0707070707060504LL, 0x0705030107040100LL,
		0x0707070707060504LL, 0x0705030107040100LL);
	const __m512i BASE8 = _mm512_set_epi32(
		0x00F800B8, 0x00780038,  0x00F000B0, 0x00700030,
		0x00E800A8, 0x00680028,  0x00E000A0, 0x00600020,
		0x00D80098, 0x00580018,  0x00D00090, 0x00500010,
		0x00C80088, 0x00480008,  0x00C00080, 0x00400000 );
	const __m512i BASE4 = _mm512_set_epi32(
		0x78, 0x70,  0x68, 0x60,  0x58, 0x50,  0x48, 0x40,
		0x38, 0x30,  0x28, 0x20,  0x18, 0x10,  0x08, 0x00 );
	// float[16] together
	__m512 sum16 = _mm512_setzero_ps();

	// 8 bytes together
	for (; n >= 8; n-=8, p_G+=256)
	{
		// packed genotypes
		const uint64_t b_i = *((const unaligned_uint64*)p_i);
		const uint64_t b_j = *((const unaligned_uint64*)p_j);
		p_i += 8; p_j += 8;
		// missing genotype if b_i & 0x03 == 0x03
		const uint64_t missing =
			(((b_i << 1) & b_i) | ((b_j << 1) & b_j)) & 0xAAAAAAAAAAAAAAAAULL;
		miss_n += __builtin_popcountll(missing);
		// bitwise
		__m512i b1 = _mm512_srlv_epi64(_mm512_set1_epi64(b_i), shr_n8) & X03;
		__m512i b2 = _mm512_srlv_epi64(_mm512_set1_epi64(b_j), shr_n8) & X03;
		// product of two normalized genotypes in the look-up table
		__m512i g = _mm512_shuffle_epi8(T16, (b2 << 2) | b1) | BASE8;
		// update sum
		__m512i i1 = g & XFF;
		__m512i i2 = _mm512_srli_epi32(g, 16) & XFF;
		__m512 s1 = _mm512_i32gather_ps(i1, p_G, 4);
		__m512 s2 = _mm512_i32gather_ps(i2, p_G, 4);
		// finally update sum
		sum16 += s1 + s2;
	}
	// 4 bytes together
	for (; n >= 4; n-=4, p_G+=128)
	{
		// packed genotypes
		const uint32_t b_i = *((const unaligned_uint32*)p_i);
		const uint32_t b_j = *((const unaligned_uint32*)p_j);
		p_i += 4; p_j += 4;
		// missing genotype if b_i & 0x03 == 0x03
		const uint32_t missing =
			(((b_i << 1) & b_i) | ((b_j << 1) & b_j)) & 0xAAAAAAAA;
		miss_n += __builtin_popcount(missing);
		// bitwise
		__m512i b1 = _mm512_srlv_epi32(_mm512_set1_epi32(b_i), shr_n4) & X03;
		__m512i b2 = _mm512_srlv_epi32(_mm512_set1_epi32(b_j), shr_n4) & X03;
		// product of two normalized genotypes in the look-up table
		__m512i ii = _mm512_shuffle_epi8(T16, (b2 << 2) | b1) & XFF;
		// update sum
		sum16 += _mm512_i32gather_ps(ii | BASE4, p_G, 4);
	}
	sum += _mm512_reduce_add_ps(sum16);
	// the remaining bytes
	for (; n > 0; n--, p_G+=32)
	{
		// packed genotypes
		const uint8_t b_i = *p_i++;
		const uint8_t b_j = *p_j++;
		// missing genotype if b_i & 0x03 == 0x03
		const uint8_t missing = (((b_i << 1) & b_i) | ((b_j << 1) & b_j)) & 0xAA;
		miss_n += __builtin_popcount(missing);
		// bitwise
		const uint8_t a_or_b = b_i | b_j;
		const uint8_t a_and_b = b_i & b_j;
		const uint8_t bit12 = (a_or_b & 0x55) |
			(((a_and_b << 1) | a_and_b | missing) & 0xAA);
		const uint8_t bit3 = a_or_b;
		// product of two normalized genotypes in the look-up table
		// geno 1
		sum += p_G[(bit12 & 0x03) | ((bit3 & 0x02) << 1)];
		// geno 2
		sum += p_G[((bit12 & 0x0C) >> 2) | ((bit3 & 0x08) >> 1) | 0x08];
		// geno 3
		sum += p_G[((bit12 & 0x30) >> 4) | ((bit3 & 0x20) >> 3) | 0x10];
		// geno 4
		sum += p_G[((bit12 & 0xC0) >> 6) | ((bit3 & 0x80) >> 5) | 0x18];
	}
}


extern "C" void SIMD_NAME(grm_calc_update_f64)(
	const uint8_t p_i[], const uint8_t p_j[], size_t n, const double p_G[],
	int &miss_n, double &sum)
{
	// constants
	const __m512i shr_n = _mm512_set_epi64(14, 12, 10, 8, 6, 4, 2, 0);
	const __m512i X03 = _mm512_set1_epi8(0x03);
	const __m512i XFF = _mm512_set1_epi64(0xFF);
	const __m512i T16 = _mm512_set_epi64(
		// look-up table instead of bitwise in vec_ext_def.cpp
		0x0707070707060504LL, 0x0705030107040100LL,
		0x0707070707060504LL, 0x0705030107040100LL,
		0x0707070707060504LL, 0x0705030107040100LL,
		0x0707070707060504LL, 0x0705030107040100LL);
	const __m512i BASE = _mm512_set_epi64(
		0x00F800B800780038LL, 0x00F000B000700030LL,
		0x00E800A800680028LL, 0x00E000A000600020LL,
		0x00D8009800580018LL, 0x00D0009000500010LL,
		0x00C8008800480008LL, 0x00C0008000400000LL);
	// double[8] together
	__m512d sum8 = _mm512_setzero_pd();

#ifdef VEC_CPU_LP64
	// 8 bytes together
	for (; n >= 8; n-=8, p_G+=256)
	{
		// packed genotypes
		const uint64_t b_i = *((const unaligned_uint64*)p_i);
		const uint64_t b_j = *((const unaligned_uint64*)p_j);
		p_i += 8; p_j += 8;
		// missing genotype if b_i & 0x03 == 0x03
		const uint64_t missing =
			(((b_i << 1) & b_i) | ((b_j << 1) & b_j)) & 0xAAAAAAAAAAAAAAAAULL;
		miss_n += __builtin_popcountll(missing);
		// bitwise
		__m512i b1 = _mm512_srlv_epi64(_mm512_set1_epi64(b_i), shr_n) & X03;
		__m512i b2 = _mm512_srlv_epi64(_mm512_set1_epi64(b_j), shr_n) & X03;
		// product of two normalized genotypes in the look-up table
		// indices for geno 1 .. 4
		__m512i g  = _mm512_shuffle_epi8(T16, (b2 << 2) | b1) | BASE;
		// update sum
		__m512i i1 = g & XFF;
		__m512i i2 = _mm512_srli_epi64(g, 16) & XFF;
		__m512i i3 = _mm512_srli_epi64(g, 32) & XFF;
		__m512i i4 = _mm512_srli_epi64(g, 48) & XFF;
		__m512d s1 = _mm512_i64gather_pd(i1, p_G, 8);
		__m512d s2 = _mm512_i64gather_pd(i2, p_G, 8);
		__m512d s3 = _mm512_i64gather_pd(i3, p_G, 8);
		__m512d s4 = _mm512_i64gather_pd(i4, p_G, 8);
		// finally update sum
		sum8 += s1 + s2 + s3 + s4;
	}
#endif
	// 4 bytes together
	for (; n >= 4; n-=4, p_G+=128)
	{
		// packed genotypes
		const uint32_t b_i = *((const unaligned_uint32*)p_i);
		const uint32_t b_j = *((const unaligned_uint32*)p_j);
		p_i += 4; p_j += 4;
		// missing genotype if b_i & 0x03 == 0x03
		const uint32_t missing =
			(((b_i << 1) & b_i) | ((b_j << 1) & b_j)) & 0xAAAAAAAA;
		miss_n += __builtin_popcount(missing);
		// bitwise
		__m512i b1 = _mm512_srlv_epi64(_mm512_set1_epi64(b_i), shr_n) & X03;
		__m512i b2 = _mm512_srlv_epi64(_mm512_set1_epi64(b_j), shr_n) & X03;
		// product of two normalized genotypes in the look-up table
		// indices for geno 1 .. 4
		__m512i g  = _mm512_shuffle_epi8(T16, (b2 << 2) | b1) | BASE;
		// update sum
		__m512i i1 = g & XFF;
		__m512i i2 = _mm512_srli_epi64(g, 16) & XFF;
		__m512d s1 = _mm512_i64gather_pd(i1, p_G, 8);
		__m512d s2 = _mm512_i64gather_pd(i2, p_G, 8);
		// finally update sum
		sum8 += s1 + s2;
	}
	sum += _mm512_reduce_add_pd(sum8);
	// the remaining bytes
	for (; n > 0; n--, p_G+=32)
	{
		// packed genotypes
		const uint8_t b_i = *p_i++;
		const uint8_t b_j = *p_j++;
		// missing genotype if b_i & 0x03 == 0x03
		const uint8_t missing = (((b_i << 1) & b_i) | ((b_j << 1) & b_j)) & 0xAA;
		miss_n += __builtin_popcount(missing);
		// bitwise
		const uint8_t a_or_b = b_i | b_j;
		const uint8_t a_and_b = b_i & b_j;
		const uint8_t bit12 = (a_or_b & 0x55) |
			(((a_and_b << 1) | a_and_b | missing) & 0xAA);
		const uint8_t bit3 = a_or_b;
		// product of two normalized genotypes in the look-up table
		// geno 1
		sum += p_G[(bit12 & 0x03) | ((bit3 & 0x02) << 1)];
		// geno 2
		sum += p_G[((bit12 & 0x0C) >> 2) | ((bit3 & 0x08) >> 1) | 0x08];
		// geno 3
		sum += p_G[((bit12 & 0x30) >> 4) | ((bit3 & 0x20) >> 3) | 0x10];
		// geno 4
		sum += p_G[((bit12 & 0xC0) >> 6) | ((bit3 & 0x80) >> 5) | 0x18];
	}
}


extern "C" double SIMD_NAME(get_dot_sp_b)(const double p_std_geno_d[],
	const double p_b[], const uint8_t p_sp_g[])
{
	const int *s256 = (const int *)p_sp_g;  // no worry about the unalign
	const uint8_t *s = p_sp_g + (*s256++);
	double dot = 0;
	// g1 * b, g2 * b, g3 * b
	for (int k=1; k <= 3; k++)
	{
		__m512d sum8 = _mm512_setzero_pd();
		for (int n256 = *s256++; n256 > 0; n256--)
		{
			const double *p = p_b + (*s256++);  // the start position
			size_t n = (*s++) + 1;
			// for-each nonzero
			for (; n >= 16; n-=16, s+=16)
			{
				__m128i z = _mm_loadu_si128((__m128i const*)s);
				__m512i i1 = _mm512_cvtepu8_epi64(z);
				__m512i i2 = _mm512_cvtepu8_epi64(_mm_shuffle_epi32(z, 0x0E));
				sum8 += _mm512_i64gather_pd(i1, p, 8) +
					_mm512_i64gather_pd(i2, p, 8);
			}
			if (n >= 8)
			{
				__m128i z; z[0] = *((unaligned_uint64*)s);
				__m512i ii = _mm512_cvtepu8_epi64(z);
				sum8 += _mm512_i64gather_pd(ii, p, 8);
				n -= 8; s += 8;
			}
			/*
			if (n >= 4)
			{
				__m128i z; z[0] = *((unaligned_uint32*)s);
				__m256i i4 = _mm256_cvtepu8_epi64(z);
				sum8 += _mm512_zextpd256_pd512(_mm256_i64gather_pd(p, i4, 8));
				n -= 4; s += 4;
			}
			*/
			for (; n > 0; n--) sum8[0] += p[*s++];
		}
		dot += _mm512_reduce_add_pd(sum8) * p_std_geno_d[k];
	}
	// output
	return dot;
}


extern "C" void SIMD_NAME(set_dot_sp_b)(double p_add[], double dot,
	const double p_std_geno_d[], const uint8_t p_sp_g[])
{
	const int *s256 = (const int *)p_sp_g;  // no worry about the unalign
	const uint8_t *s = p_sp_g + (*s256++);
	// g1 * dot, g2 * dot, g3 * dot
	for (int k=1; k <= 3; k++)
	{
		const __m512d v = _mm512_set1_pd(dot * p_std_geno_d[k]);
		for (int n256 = *s256++; n256 > 0; n256--)
		{
			double *p = p_add + (*s256++);  // the start position
			size_t n = (*s++) + 1;
			// for-each nonzero
			for (; n >= 8; n-=8, s+=8)
			{
				__m128i z; z[0] = *((unaligned_uint64*)s);
				__m512i ii = _mm512_cvtepu8_epi64(z);
				_mm512_i64scatter_pd(p, ii, v+_mm512_i64gather_pd(ii, p, 8), 8);
			}
			for (; n > 0; n--) p[*s++] += v[0];
		}
	}
}


#else

extern "C" void SIMD_NAME(grm_calc_update_f32)(
	const uint8_t p_i[], const uint8_t p_j[], size_t n, const float p_G[],
	int &miss_n, double &sum)
{
	THROW_ERROR;
}

extern "C" void SIMD_NAME(grm_calc_update_f64)(
	const uint8_t p_i[], const uint8_t p_j[], size_t n, const double p_G[],
	int &miss_n, double &sum)
{
	THROW_ERROR;
}

extern "C" double SIMD_NAME(get_dot_sp_b)(const double p_std_geno_d[],
	const double p_b[], const uint8_t p_sp_g[])
{
	THROW_ERROR;
}

extern "C" void SIMD_NAME(set_dot_sp_b)(double p_add[], double dot,
	const double p_std_geno_d[], const uint8_t p_sp_g[])
{
	THROW_ERROR;
}

#endif
