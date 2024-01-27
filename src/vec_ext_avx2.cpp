// ===============================================================
//
// vec_ext_avx2.cpp: optimization with AVX2 vectorization
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

#include "vectorization.h"
#include "vec_ext.h"

// need a patch for gcc_v4.8
#if defined(VEC_CPU_ARCH_X86) && defined(__GNUC__) && (__GNUC__==4) && (__GNUC_MINOR__==8)
#   pragma GCC target("avx2")
#   define __AVX2__      1
#   define __AVX__       1
#   define __SSE4_1__    1
#   define __SSE4_2__    1
#   define __SSE3__      1
#   define __SSSE3__     1
#   define __POPCNT__    1
#   include <xmmintrin.h>  // SSE
#   include <emmintrin.h>  // SSE2
#   include <immintrin.h>  // AVX
#   define TARGET_AVX2
#endif


#ifdef VEC_CPU_ARCH_X86_AVX2
extern const bool VEC_ALGORITHM_AVX2 = true;
#else
extern const bool VEC_ALGORITHM_AVX2 = false;
#endif


#define SIMD_NAME(NAME)  NAME ## _avx2
#define THROW_ERROR      throw "No AVX2 support!"


#ifdef VEC_CPU_ARCH_X86_AVX2

#ifdef __ICC
	#pragma intel optimization_parameter target_arch=CORE-AVX2
#elif !defined(__AVX2__) && !defined(__clang__)
	#pragma GCC target("avx2")
#endif

#include <xmmintrin.h>  // SSE
#include <emmintrin.h>  // SSE2
#include <immintrin.h>  // AVX, AVX2

#ifndef TARGET_AVX2
#   define TARGET_AVX2    __attribute__((target("avx2")))
#endif

#undef SIMD_NAME
#define SIMD_NAME(NAME)    MATH_O3 TARGET_AVX2 NAME ## _avx2


extern "C" void SIMD_NAME(grm_calc_update_f32)(
	const uint8_t p_i[], const uint8_t p_j[], size_t n, const float p_G[],
	int &miss_n, double &sum)
{
	// constants
	const __m256i shr_n4 = _mm256_set_epi32(14, 12, 10, 8, 6, 4, 2, 0);
	const __m256i shr_n8 = _mm256_set_epi64x(6, 4, 2, 0);
	const __m256i X03 = _mm256_set1_epi8(0x03);
	// const __m256i XAA = _mm256_set1_epi8(0xAA);
	// const __m256i X08 = _mm256_set1_epi8(0x08);
	const __m256i XFF = _mm256_set1_epi32(0xFF);
	// const __m256i BASE32 = _mm256_set1_epi64x(0xE0C0A08060402000LL);
	// const __m256i BASE32_2 = _mm256_set_epi32(
	//	0x300, 0x300, 0x200, 0x200, 0x100, 0x100, 0x000, 0x000);
	const __m256i BASE8 = _mm256_set_epi64x(
		0xF8D8B89878583818LL, 0xF0D0B09070503010LL,
		0xE8C8A88868482808LL, 0xE0C0A08060402000LL);
	const __m256i BASE4 = _mm256_set_epi32(
		0x00780038, 0x00700030,  0x00680028, 0x00600020,
		0x00580018, 0x00500010,  0x00480008, 0x00400000 );
	const __m256i T16 = _mm256_set_epi64x(
		// look-up table instead of bitwise in vec_ext_def.cpp
		0x0707070707060504LL, 0x0705030107040100LL,
		0x0707070707060504LL, 0x0705030107040100LL);
	// double[8] together
	__m256 sum8 = _mm256_setzero_ps();

/*
	// 32 bytes together
	for (; n >= 32; n-=32, p_G+=1024)
	{
		// packed genotypes
		const __m256i b1 = _mm256_loadu_si256((const __m256i*)p_i);
		const __m256i b2 = _mm256_loadu_si256((const __m256i*)p_j);
		p_i += 32; p_j += 32;
		// missing genotype if b_i & 0x03 == 0x03
		const __m256i missing =
			(((b1 << 1) & b1) | ((b2 << 1) & b2)) & XAA;
		miss_n += __builtin_popcountll(missing[0]) +
			__builtin_popcountll(missing[1]) +
			__builtin_popcountll(missing[2]) +
			__builtin_popcountll(missing[3]);
		// product of two normalized genotypes in the look-up table
		// indices for geno 1 .. 4
		__m256i g1 = ((b2 & X03) << 2) | (b1 & X03);
		__m256i g2 = (((b2 >> 2) & X03) << 2) | ((b1 >> 2) & X03);
		__m256i g3 = (((b2 >> 4) & X03) << 2) | ((b1 >> 4) & X03);
		__m256i g4 = (((b2 >> 6) & X03) << 2) | ((b1 >> 6) & X03);
		g1 = _mm256_shuffle_epi8(T16, g1) | BASE32;
		g2 = _mm256_shuffle_epi8(T16, g2) | BASE32 | X08;
		g3 = _mm256_shuffle_epi8(T16, g3) | BASE32 | (X08 << 1);
		g4 = _mm256_shuffle_epi8(T16, g4) | BASE32 | X08 | (X08 << 1);
		// shuffle
		...
		// update sum
		__m256i i1, i2, i3, i4;
		// byte 1
		i1 = (g1 & XFF) | BASE32_2;
		i2 = (g2 & XFF) | BASE32_2;
		i3 = (g3 & XFF) | BASE32_2;
		i4 = (g4 & XFF) | BASE32_2;
		__m256 s1 = _mm256_i32gather_ps(p_G, i1, 4) + _mm256_i32gather_ps(p_G, i2, 4) +
			_mm256_i32gather_ps(p_G, i3, 4) + _mm256_i32gather_ps(p_G, i4, 4);
		// byte 2
		i1 = (_mm256_srli_epi32(g1, 8) & XFF) | BASE32_2;
		i2 = (_mm256_srli_epi32(g2, 8) & XFF) | BASE32_2;
		i3 = (_mm256_srli_epi32(g3, 8) & XFF) | BASE32_2;
		i4 = (_mm256_srli_epi32(g4, 8) & XFF) | BASE32_2;
		__m256 s2 = _mm256_i32gather_ps(p_G, i1, 4) + _mm256_i32gather_ps(p_G, i2, 4) +
			_mm256_i32gather_ps(p_G, i3, 4) + _mm256_i32gather_ps(p_G, i4, 4);
		// byte 3
		i1 = (_mm256_srli_epi32(g1, 16) & XFF) | BASE32_2;
		i2 = (_mm256_srli_epi32(g2, 16) & XFF) | BASE32_2;
		i3 = (_mm256_srli_epi32(g3, 16) & XFF) | BASE32_2;
		i4 = (_mm256_srli_epi32(g4, 16) & XFF) | BASE32_2;
		__m256 s3 = _mm256_i32gather_ps(p_G, i1, 4) + _mm256_i32gather_ps(p_G, i2, 4) +
			_mm256_i32gather_ps(p_G, i3, 4) + _mm256_i32gather_ps(p_G, i4, 4);
		// byte 4
		i1 = (_mm256_srli_epi32(g1, 24) & XFF) | BASE32_2;
		i2 = (_mm256_srli_epi32(g2, 24) & XFF) | BASE32_2;
		i3 = (_mm256_srli_epi32(g3, 24) & XFF) | BASE32_2;
		i4 = (_mm256_srli_epi32(g4, 24) & XFF) | BASE32_2;
		__m256 s4 = _mm256_i32gather_ps(p_G, i1, 4) + _mm256_i32gather_ps(p_G, i2, 4) +
			_mm256_i32gather_ps(p_G, i3, 4) + _mm256_i32gather_ps(p_G, i4, 4);
		// finally update sum
		sum8 += s1 + s2 + s3 + s4;
	}
*/

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
		__m256i b1 = _mm256_srlv_epi64(_mm256_set1_epi64x(b_i), shr_n8) & X03;
		__m256i b2 = _mm256_srlv_epi64(_mm256_set1_epi64x(b_j), shr_n8) & X03;
		// product of two normalized genotypes in the look-up table
		// indices for geno 1 .. 4
		__m256i g  = _mm256_shuffle_epi8(T16, (b2 << 2) | b1) | BASE8;
		// update sum
		__m256i i1 = g & XFF;
		__m256i i2 = _mm256_srli_epi32(g, 8) & XFF;
		__m256i i3 = _mm256_srli_epi32(g, 16) & XFF;
		__m256i i4 = _mm256_srli_epi32(g, 24) & XFF;
		__m256 s1 = _mm256_i32gather_ps(p_G, i1, 4);
		__m256 s2 = _mm256_i32gather_ps(p_G, i2, 4);
		__m256 s3 = _mm256_i32gather_ps(p_G, i3, 4);
		__m256 s4 = _mm256_i32gather_ps(p_G, i4, 4);
		// finally update sum
		sum8 += s1 + s2 + s3 + s4;
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
		__m256i b1 = _mm256_srlv_epi32(_mm256_set1_epi32(b_i), shr_n4) & X03;
		__m256i b2 = _mm256_srlv_epi32(_mm256_set1_epi32(b_j), shr_n4) & X03;
		// product of two normalized genotypes in the look-up table
		// indices for geno 1 .. 4
		__m256i g  = _mm256_shuffle_epi8(T16, (b2 << 2) | b1) | BASE4;
		// update sum
		__m256i i1 = g & XFF;
		__m256i i2 = _mm256_srli_epi32(g, 16) & XFF;
		__m256 s1 = _mm256_i32gather_ps(p_G, i1, 4);
		__m256 s2 = _mm256_i32gather_ps(p_G, i2, 4);
		// finally update sum
		sum8 += s1 + s2;
	}
	sum += sum8[0] + sum8[1] + sum8[2] + sum8[3] +
		sum8[4] + sum8[5] + sum8[6] + sum8[7];
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
	const __m256i shr_n = _mm256_set_epi64x(6, 4, 2, 0);
	const __m256i X03 = _mm256_set1_epi8(0x03);
	const __m256i XFF = _mm256_set1_epi64x(0xFF);
	const __m256i BASE = _mm256_set_epi64x(
		0xF8D8B89878583818LL, 0xF0D0B09070503010LL,
		0xE8C8A88868482808LL, 0xE0C0A08060402000LL);
	const __m256i T16 = _mm256_set_epi64x(
		// look-up table instead of bitwise in vec_ext_def.cpp
		0x0707070707060504LL, 0x0705030107040100LL,
		0x0707070707060504LL, 0x0705030107040100LL);
	// double[4] together
	__m256d sum4 = _mm256_setzero_pd();

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
		__m256i b1 = _mm256_srlv_epi64(_mm256_set1_epi64x(b_i), shr_n) & X03;
		__m256i b2 = _mm256_srlv_epi64(_mm256_set1_epi64x(b_j), shr_n) & X03;
		// product of two normalized genotypes in the look-up table
		// indices for geno 1 .. 4
		__m256i g  = _mm256_shuffle_epi8(T16, (b2 << 2) | b1) | BASE;
		// update sum
		__m256i i1 = g & XFF;
		__m256i i2 = _mm256_srli_epi64(g, 8) & XFF;
		__m256i i3 = _mm256_srli_epi64(g, 16) & XFF;
		__m256i i4 = _mm256_srli_epi64(g, 24) & XFF;
		__m256i i5 = _mm256_srli_epi64(g, 32) & XFF;
		__m256i i6 = _mm256_srli_epi64(g, 40) & XFF;
		__m256i i7 = _mm256_srli_epi64(g, 48) & XFF;
		__m256i i8 = _mm256_srli_epi64(g, 56) & XFF;
		__m256d s1 = _mm256_i64gather_pd(p_G, i1, 8);
		__m256d s2 = _mm256_i64gather_pd(p_G, i2, 8);
		__m256d s3 = _mm256_i64gather_pd(p_G, i3, 8);
		__m256d s4 = _mm256_i64gather_pd(p_G, i4, 8);
		__m256d s5 = _mm256_i64gather_pd(p_G, i5, 8);
		__m256d s6 = _mm256_i64gather_pd(p_G, i6, 8);
		__m256d s7 = _mm256_i64gather_pd(p_G, i7, 8);
		__m256d s8 = _mm256_i64gather_pd(p_G, i8, 8);
		// finally update sum
		sum4 += s1 + s2 + s3 + s4 + s5 + s6 + s7 + s8;
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
		__m256i b1 = _mm256_srlv_epi64(_mm256_set1_epi64x(b_i), shr_n) & X03;
		__m256i b2 = _mm256_srlv_epi64(_mm256_set1_epi64x(b_j), shr_n) & X03;
		// product of two normalized genotypes in the look-up table
		// indices for geno 1 .. 4
		__m256i g  = _mm256_shuffle_epi8(T16, (b2 << 2) | b1) | BASE;
		// update sum
		__m256i i1 = g & XFF;
		__m256i i2 = _mm256_srli_epi64(g, 8) & XFF;
		__m256i i3 = _mm256_srli_epi64(g, 16) & XFF;
		__m256i i4 = _mm256_srli_epi64(g, 24) & XFF;
		__m256d s1 = _mm256_i64gather_pd(p_G, i1, 8);
		__m256d s2 = _mm256_i64gather_pd(p_G, i2, 8);
		__m256d s3 = _mm256_i64gather_pd(p_G, i3, 8);
		__m256d s4 = _mm256_i64gather_pd(p_G, i4, 8);
		// finally update sum
		sum4 += s1 + s2 + s3 + s4;
	}
	sum += sum4[0] + sum4[1] + sum4[2] + sum4[3];
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
		__m256d sum4 = _mm256_setzero_pd();
		for (int n256 = *s256++; n256 > 0; n256--)
		{
			const double *p = p_b + (*s256++);  // the start position
			size_t n = (*s++) + 1;
			// for-each nonzero
			for (; n >= 16; n-=16, s+=16)
			{
				__m128i z = _mm_loadu_si128((__m128i const*)s);
				__m256i i1 = _mm256_cvtepu8_epi64(z);
				__m256i i2 = _mm256_cvtepu8_epi64(_mm_shuffle_epi32(z, 1));
				__m256i i3 = _mm256_cvtepu8_epi64(_mm_shuffle_epi32(z, 2));
				__m256i i4 = _mm256_cvtepu8_epi64(_mm_shuffle_epi32(z, 3));
				sum4 += _mm256_i64gather_pd(p, i1, 8) +
					_mm256_i64gather_pd(p, i2, 8) +
					_mm256_i64gather_pd(p, i3, 8) +
					_mm256_i64gather_pd(p, i4, 8);
			}
			for (; n >= 4; n-=4, s+=4)
			{
				__m128i z; z[0] = *((unaligned_uint32*)s);
				__m256i i4 = _mm256_cvtepu8_epi64(z);
				sum4 += _mm256_i64gather_pd(p, i4, 8);
			}
			for (; n > 0; n--) sum4[0] += p[*s++];
		}
		dot += (sum4[0]+sum4[1]+sum4[2]+sum4[3]) * p_std_geno_d[k];
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
		const __m256d v4 = _mm256_set1_pd(dot * p_std_geno_d[k]);
		for (int n256 = *s256++; n256 > 0; n256--)
		{
			double *p = p_add + (*s256++);  // the start position
			size_t n = (*s++) + 1;
			// for-each nonzero
			for (; n >= 4; n-=4, s+=4)
			{
				__m128i z; z[0] = *((unaligned_uint32*)s);
				__m256i i4 = _mm256_cvtepu8_epi64(z);
				__m256d a = _mm256_i64gather_pd(p, i4, 8) + v4;
				p[i4[0]] = a[0]; p[i4[1]] = a[1];
				p[i4[2]] = a[2]; p[i4[3]] = a[3];
			}
			for (; n > 0; n--) p[*s++] += v4[0];
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
