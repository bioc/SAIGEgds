// ===============================================================
//
// vec_ext_def.cpp: optimization with the default setting
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


extern "C"
{
	bool fc_use_avx512f = true;
	bool fc_use_avx2 = true;
}

/// Initialize the target-specific functions
extern "C" void vec_init_function()
{
	extern const bool VEC_ALGORITHM_AVX2;
	extern const bool VEC_ALGORITHM_AVX512BW;

	// initialized the defaults
	fc_grm_calc_update_f32 = &grm_calc_update_f32_def;
	fc_grm_calc_update_f64 = &grm_calc_update_f64_def;
	fc_get_dot_sp_b = &get_dot_sp_b_def;
	fc_set_dot_sp_b = &set_dot_sp_b_def;

#ifdef VEC_CPU_ARCH_X86

#ifdef VEC_BUILTIN_CPU
	__builtin_cpu_init();
#endif

// ========  Intel AVX512BW  ========
#if defined(__AVX512F__) && defined(__AVX512BW__)
	const bool has_avx512bw = true;
#elif defined(VEC_BUILTIN_CPU_AVX512BW)
#if defined(__ICC) && defined(_FEATURE_AVX512F) && defined(_FEATURE_AVX512BW)
	// since __builtin_cpu_supports("avx512bw") in ICC always return 0
	const bool has_avx512bw = _may_i_use_cpu_feature(
		_FEATURE_AVX512F | _FEATURE_AVX512BW) && VEC_ALGORITHM_AVX512BW;
#else
	const bool has_avx512bw = __builtin_cpu_supports("avx512f") &&
		__builtin_cpu_supports("avx512bw") && VEC_ALGORITHM_AVX512BW;
#endif
#else
	const bool has_avx512bw = false;
#endif

/*
// ========  Intel AVX512F  ========
#if defined(__AVX512F__)
	const bool has_avx512f = true;
#elif defined(VEC_BUILTIN_CPU_AVX512F)
	const bool has_avx512f = __builtin_cpu_supports("avx512f") &&
		VEC_ALGORITHM_AVX512F;
#else
	const bool has_avx512f = false;
#endif
*/

// ========  Intel AVX2  ========
#if defined(__AVX2__)
	const bool has_avx2 = true;
#elif defined(VEC_BUILTIN_CPU_AVX2)
	const bool has_avx2 = __builtin_cpu_supports("avx2") && VEC_ALGORITHM_AVX2;
#else
	const bool has_avx2 = false;
#endif

	// set the functions
	if (has_avx512bw && fc_use_avx512f)
	{
		fc_grm_calc_update_f32 = &grm_calc_update_f32_avx512bw;
		fc_grm_calc_update_f64 = &grm_calc_update_f64_avx512bw;
		fc_get_dot_sp_b = &get_dot_sp_b_avx512bw;
		fc_set_dot_sp_b = &set_dot_sp_b_avx512bw;
	} else if (has_avx2 && fc_use_avx2)
	{
		fc_grm_calc_update_f32 = &grm_calc_update_f32_avx2;
		fc_grm_calc_update_f64 = &grm_calc_update_f64_avx2;
		fc_get_dot_sp_b = &get_dot_sp_b_avx2;
		fc_set_dot_sp_b = &set_dot_sp_b_avx2;
	}
#endif
}


extern "C"
{
	// function pointers for GRM calculation
	type_grm_calc_update_f32 fc_grm_calc_update_f32 = NULL;
	type_grm_calc_update_f64 fc_grm_calc_update_f64 = NULL;

	// function pointers for the crossproduct of GRM and a vector b
	type_get_dot_sp_b fc_get_dot_sp_b = NULL;
	type_set_dot_sp_b fc_set_dot_sp_b = NULL;
}

/// Update GRM values using a look-up table
extern "C" MATH_OFAST
void grm_calc_update_f32_def(const uint8_t p_i[], const uint8_t p_j[],
	size_t n, const float p_G[], int &miss_n, double &sum)
{
#ifdef VEC_CPU_LP64
	// 8 bytes together
	for (; n >= 8; n-=8, p_G+=256)
	{
		// packed genotypes
		const uint64_t b_i = *((const unaligned_uint64*)p_i);
		const uint64_t b_j = *((const unaligned_uint64*)p_j);
		p_i += 8; p_j += 8;
		// using bitwise operations
		// missing genotype if b_i & 0x03 == 0x03
		const uint64_t missing =
			(((b_i << 1) & b_i) | ((b_j << 1) & b_j)) & 0xAAAAAAAAAAAAAAAAULL;
		miss_n += __builtin_popcountll(missing);
		// bitwise
		const uint64_t a_or_b = b_i | b_j;
		const uint64_t a_and_b = b_i & b_j;
		const uint64_t bit12 = (a_or_b & 0x5555555555555555ULL) |
			(((a_and_b << 1) | a_and_b | missing) & 0xAAAAAAAAAAAAAAAAULL);
		const uint64_t bit3 = a_or_b;
		// product of two normalized genotypes in the look-up table
		// indices for geno 1 .. 4
		const uint64_t
			g1 = (bit12 & 0x0303030303030303ULL) |
				((bit3 & 0x0202020202020202ULL) << 1) | 0xE0C0A08060402000ULL,
			g2 = ((bit12 & 0x0C0C0C0C0C0C0C0CULL) >> 2) |
				((bit3 & 0x0808080808080808ULL) >> 1) | 0xE8C8A88868482808ULL,
			g3 = ((bit12 & 0x3030303030303030ULL) >> 4) |
				((bit3 & 0x2020202020202020ULL) >> 3) | 0xF0D0B09070503010ULL,
			g4 = ((bit12 & 0xC0C0C0C0C0C0C0C0ULL) >> 6) |
				((bit3 & 0x8080808080808080ULL) >> 5) | 0xF8D8B89878583818ULL;
		// update sum
		float s1, s2, s3, s4;
		// byte 1
		s1 = p_G[g1 & 0xFF];  s2 = p_G[g2 & 0xFF];
		s3 = p_G[g3 & 0xFF];  s4 = p_G[g4 & 0xFF];
		// byte 2
		s1 += p_G[(g1 >> 8) & 0xFF];  s2 += p_G[(g2 >> 8) & 0xFF];
		s3 += p_G[(g3 >> 8) & 0xFF];  s4 += p_G[(g4 >> 8) & 0xFF];
		// byte 3
		s1 += p_G[(g1 >> 16) & 0xFF];  s2 += p_G[(g2 >> 16) & 0xFF];
		s3 += p_G[(g3 >> 16) & 0xFF];  s4 += p_G[(g4 >> 16) & 0xFF];
		// byte 4
		s1 += p_G[(g1 >> 24) & 0xFF];  s2 += p_G[(g2 >> 24) & 0xFF];
		s3 += p_G[(g3 >> 24) & 0xFF];  s4 += p_G[(g4 >> 24) & 0xFF];
		// byte 5
		s1 += p_G[(g1 >> 32) & 0xFF];  s2 += p_G[(g2 >> 32) & 0xFF];
		s3 += p_G[(g3 >> 32) & 0xFF];  s4 += p_G[(g4 >> 32) & 0xFF];
		// byte 6
		s1 += p_G[(g1 >> 40) & 0xFF];  s2 += p_G[(g2 >> 40) & 0xFF];
		s3 += p_G[(g3 >> 40) & 0xFF];  s4 += p_G[(g4 >> 40) & 0xFF];
		// byte 7
		s1 += p_G[(g1 >> 48) & 0xFF];  s2 += p_G[(g2 >> 48) & 0xFF];
		s3 += p_G[(g3 >> 48) & 0xFF];  s4 += p_G[(g4 >> 48) & 0xFF];
		// byte 8
		s1 += p_G[g1 >> 56];  s2 += p_G[g2 >> 56];
		s3 += p_G[g3 >> 56];  s4 += p_G[g4 >> 56];
		// finally update sum
		sum += s1 + s2 + s3 + s4;
	}
#endif
	// 4 bytes together
	for (; n >= 4; n-=4, p_G+=128)
	{
		// packed genotypes
		const uint32_t b_i = *((const unaligned_uint32*)p_i);
		const uint32_t b_j = *((const unaligned_uint32*)p_j);
		p_i += 4; p_j += 4;
		// using bitwise operations
		// missing genotype if b_i & 0x03 == 0x03
		const uint32_t missing =
			(((b_i << 1) & b_i) | ((b_j << 1) & b_j)) & 0xAAAAAAAA;
		miss_n += __builtin_popcount(missing);
		// bitwise
		const uint32_t a_or_b = b_i | b_j;
		const uint32_t a_and_b = b_i & b_j;
		const uint32_t bit12 = (a_or_b & 0x55555555) |
			(((a_and_b << 1) | a_and_b | missing) & 0xAAAAAAAA);
		const uint32_t bit3 = a_or_b;
		// product of two normalized genotypes in the look-up table
		// indices for geno 1 .. 4
		const uint32_t
			g1 = (bit12 & 0x03030303) | ((bit3 & 0x02020202) << 1) | 0x60402000,
			g2 = ((bit12 & 0x0C0C0C0C) >> 2) | ((bit3 & 0x08080808) >> 1) | 0x68482808,
			g3 = ((bit12 & 0x30303030) >> 4) | ((bit3 & 0x20202020) >> 3) | 0x70503010,
			g4 = ((bit12 & 0xC0C0C0C0) >> 6) | ((bit3 & 0x80808080) >> 5) | 0x78583818;
		// update sum
		float s1, s2, s3, s4;
		// byte 1
		s1 = p_G[g1 & 0xFF];  s2 = p_G[g2 & 0xFF];
		s3 = p_G[g3 & 0xFF];  s4 = p_G[g4 & 0xFF];
		// byte 2
		s1 += p_G[(g1 >> 8) & 0xFF];  s2 += p_G[(g2 >> 8) & 0xFF];
		s3 += p_G[(g3 >> 8) & 0xFF];  s4 += p_G[(g4 >> 8) & 0xFF];
		// byte 3
		s1 += p_G[(g1 >> 16) & 0xFF];  s2 += p_G[(g2 >> 16) & 0xFF];
		s3 += p_G[(g3 >> 16) & 0xFF];  s4 += p_G[(g4 >> 16) & 0xFF];
		// byte 4
		s1 += p_G[g1 >> 24];  s2 += p_G[g2 >> 24];
		s3 += p_G[g3 >> 24];  s4 += p_G[g4 >> 24];
		// finally update sum
		sum += s1 + s2 + s3 + s4;
	}
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


/// Update GRM values using a look-up table
extern "C" MATH_OFAST
void grm_calc_update_f64_def(const uint8_t p_i[], const uint8_t p_j[],
	size_t n, const double p_G[], int &miss_n, double &sum)
{
#ifdef VEC_CPU_LP64
	// 8 bytes together
	for (; n >= 8; n-=8, p_G+=256)
	{
		// packed genotypes
		const uint64_t b_i = *((const unaligned_uint64*)p_i);
		const uint64_t b_j = *((const unaligned_uint64*)p_j);
		p_i += 8; p_j += 8;
		// using bitwise operations
		// missing genotype if b_i & 0x03 == 0x03
		const uint64_t missing =
			(((b_i << 1) & b_i) | ((b_j << 1) & b_j)) & 0xAAAAAAAAAAAAAAAAULL;
		miss_n += __builtin_popcountll(missing);
		// bitwise
		const uint64_t a_or_b = b_i | b_j;
		const uint64_t a_and_b = b_i & b_j;
		const uint64_t bit12 = (a_or_b & 0x5555555555555555ULL) |
			(((a_and_b << 1) | a_and_b | missing) & 0xAAAAAAAAAAAAAAAAULL);
		const uint64_t bit3 = a_or_b;
		// product of two normalized genotypes in the look-up table
		// indices for geno 1 .. 4
		const uint64_t
			g1 = (bit12 & 0x0303030303030303ULL) |
				((bit3 & 0x0202020202020202ULL) << 1) | 0xE0C0A08060402000ULL,
			g2 = ((bit12 & 0x0C0C0C0C0C0C0C0CULL) >> 2) |
				((bit3 & 0x0808080808080808ULL) >> 1) | 0xE8C8A88868482808ULL,
			g3 = ((bit12 & 0x3030303030303030ULL) >> 4) |
				((bit3 & 0x2020202020202020ULL) >> 3) | 0xF0D0B09070503010ULL,
			g4 = ((bit12 & 0xC0C0C0C0C0C0C0C0ULL) >> 6) |
				((bit3 & 0x8080808080808080ULL) >> 5) | 0xF8D8B89878583818ULL;
		// update sum
		double s1, s2, s3, s4;
		// byte 1
		s1 = p_G[g1 & 0xFF];  s2 = p_G[g2 & 0xFF];
		s3 = p_G[g3 & 0xFF];  s4 = p_G[g4 & 0xFF];
		// byte 2
		s1 += p_G[(g1 >> 8) & 0xFF];  s2 += p_G[(g2 >> 8) & 0xFF];
		s3 += p_G[(g3 >> 8) & 0xFF];  s4 += p_G[(g4 >> 8) & 0xFF];
		// byte 3
		s1 += p_G[(g1 >> 16) & 0xFF];  s2 += p_G[(g2 >> 16) & 0xFF];
		s3 += p_G[(g3 >> 16) & 0xFF];  s4 += p_G[(g4 >> 16) & 0xFF];
		// byte 4
		s1 += p_G[(g1 >> 24) & 0xFF];  s2 += p_G[(g2 >> 24) & 0xFF];
		s3 += p_G[(g3 >> 24) & 0xFF];  s4 += p_G[(g4 >> 24) & 0xFF];
		// byte 5
		s1 += p_G[(g1 >> 32) & 0xFF];  s2 += p_G[(g2 >> 32) & 0xFF];
		s3 += p_G[(g3 >> 32) & 0xFF];  s4 += p_G[(g4 >> 32) & 0xFF];
		// byte 6
		s1 += p_G[(g1 >> 40) & 0xFF];  s2 += p_G[(g2 >> 40) & 0xFF];
		s3 += p_G[(g3 >> 40) & 0xFF];  s4 += p_G[(g4 >> 40) & 0xFF];
		// byte 7
		s1 += p_G[(g1 >> 48) & 0xFF];  s2 += p_G[(g2 >> 48) & 0xFF];
		s3 += p_G[(g3 >> 48) & 0xFF];  s4 += p_G[(g4 >> 48) & 0xFF];
		// byte 8
		s1 += p_G[g1 >> 56];  s2 += p_G[g2 >> 56];
		s3 += p_G[g3 >> 56];  s4 += p_G[g4 >> 56];
		// finally update sum
		sum += s1 + s2 + s3 + s4;
	}
#endif
	// 4 bytes together
	for (; n >= 4; n-=4, p_G+=128)
	{
		// packed genotypes
		const uint32_t b_i = *((const unaligned_uint32*)p_i);
		const uint32_t b_j = *((const unaligned_uint32*)p_j);
		p_i += 4; p_j += 4;
		// using bitwise operations
		// missing genotype if b_i & 0x03 == 0x03
		const uint32_t missing =
			(((b_i << 1) & b_i) | ((b_j << 1) & b_j)) & 0xAAAAAAAA;
		miss_n += __builtin_popcount(missing);
		// bitwise
		const uint32_t a_or_b = b_i | b_j;
		const uint32_t a_and_b = b_i & b_j;
		const uint32_t bit12 = (a_or_b & 0x55555555) |
			(((a_and_b << 1) | a_and_b | missing) & 0xAAAAAAAA);
		const uint32_t bit3 = a_or_b;
		// product of two normalized genotypes in the look-up table
		// indices for geno 1 .. 4
		const uint32_t
			g1 = (bit12 & 0x03030303) | ((bit3 & 0x02020202) << 1) | 0x60402000,
			g2 = ((bit12 & 0x0C0C0C0C) >> 2) | ((bit3 & 0x08080808) >> 1) | 0x68482808,
			g3 = ((bit12 & 0x30303030) >> 4) | ((bit3 & 0x20202020) >> 3) | 0x70503010,
			g4 = ((bit12 & 0xC0C0C0C0) >> 6) | ((bit3 & 0x80808080) >> 5) | 0x78583818;
		// update sum
		double s1, s2, s3, s4;
		// byte 1
		s1 = p_G[g1 & 0xFF];  s2 = p_G[g2 & 0xFF];
		s3 = p_G[g3 & 0xFF];  s4 = p_G[g4 & 0xFF];
		// byte 2
		s1 += p_G[(g1 >> 8) & 0xFF];  s2 += p_G[(g2 >> 8) & 0xFF];
		s3 += p_G[(g3 >> 8) & 0xFF];  s4 += p_G[(g4 >> 8) & 0xFF];
		// byte 3
		s1 += p_G[(g1 >> 16) & 0xFF];  s2 += p_G[(g2 >> 16) & 0xFF];
		s3 += p_G[(g3 >> 16) & 0xFF];  s4 += p_G[(g4 >> 16) & 0xFF];
		// byte 4
		s1 += p_G[g1 >> 24];  s2 += p_G[g2 >> 24];
		s3 += p_G[g3 >> 24];  s4 += p_G[g4 >> 24];
		// finally update sum
		sum += s1 + s2 + s3 + s4;
	}
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


// =========================================================================
// Get and set dot using sparse genotypes in get_crossprod_b_grm()

extern "C" MATH_OFAST
double get_dot_sp_b_def(const double p_std_geno_d[], const double p_b[],
	const uint8_t p_sp_g[])
{
	const int *s256 = (const int *)p_sp_g;  // no worry about the unalign
	const uint8_t *s = p_sp_g + (*s256++);
	double dot = 0;
	// g1 * b, g2 * b, g3 * b
	for (int k=1; k <= 3; k++)
	{
		double sum = 0;
		for (int n256 = *s256++; n256 > 0; n256--)
		{
			const double *p = p_b + (*s256++);  // the start position
			const size_t n = (*s++) + 1;
			for (size_t j=0; j < n; j++)  // for-each nonzero
				sum += p[s[j]];
			s += n;
		}
		dot += sum * p_std_geno_d[k];
	}
	// output
	return dot;
}


extern "C" MATH_OFAST
void set_dot_sp_b_def(double p_add[], double dot, const double p_std_geno_d[],
	const uint8_t p_sp_g[])
{
	const int *s256 = (const int *)p_sp_g;  // no worry about the unalign
	const uint8_t *s = p_sp_g + (*s256++);
	// g1 * dot, g2 * dot, g3 * dot
	for (int k=1; k <= 3; k++)
	{
		const double v = dot * p_std_geno_d[k];
		for (int n256 = *s256++; n256 > 0; n256--)
		{
			double *p = p_add + (*s256++);  // the start position
			const size_t n = (*s++) + 1;
			for (size_t j=0; j < n; j++)  // for-each nonzero
				p[s[j]] += v;
			s += n;
		}
	}
}

