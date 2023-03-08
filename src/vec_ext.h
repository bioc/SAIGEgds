// ===============================================================
//
// vec_ext.h: optimization with vectorization
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

#ifndef VEC_EXT_H_
#define VEC_EXT_H_


// Detect whether x86 microprocessor architecture or not
#if defined(__i386__) || defined(__X86__) || defined(_M_IX86) || defined(__I86__) || defined(__INTEL__) || defined(__amd64__) || defined(__x86_64__) || defined(_M_AMD64)
#   define VEC_CPU_ARCH_X86
#endif


// 32-bit or 64-bit registers
#ifdef __LP64__
#   define VEC_CPU_LP64
#else
#   ifdef VEC_CPU_LP64
#      undef VEC_CPU_LP64
#   endif
#endif


// remove VEC_CPU_ARCH_X86 if define VEC_NO_X86_SIMD
#if defined(VEC_CPU_ARCH_X86) && defined(VEC_NO_X86_SIMD)
#   undef VEC_CPU_ARCH_X86
#endif

// whether has __builtin_cpu_supports or not
#if defined(__GNUC__) && ((__GNUC__>4) || (__GNUC__==4 && __GNUC_MINOR__>=8))
#   define VEC_BUILTIN_CPU
#elif defined(__clang__) && defined(__clang_major__) && defined(__clang_minor__) && ((__clang_major__>3) || (__clang_major__==3 && __clang_minor__>=7))
#   define VEC_BUILTIN_CPU
#endif



// SSE2
#ifdef VEC_CPU_ARCH_X86
#   if defined(__GNUC__) && ((__GNUC__>3) || (__GNUC__==3 && __GNUC_MINOR__>=1))
#       define VEC_CPU_ARCH_X86_SSE2
#   elif defined(__clang_major__) && defined(__clang_minor__) && ((__clang_major__>3) || (__clang_major__==3 && __clang_minor__>=3))
#       define VEC_CPU_ARCH_X86_SSE2
#   endif
#   ifdef VEC_BUILTIN_CPU
#       define VEC_BUILTIN_CPU_SSE2
#   endif
#endif
#if defined(__SSE2__) && !defined(VEC_CPU_ARCH_X86_SSE2)
#   define VEC_CPU_ARCH_X86_SSE2
#endif


// SSE4.2
#ifdef VEC_CPU_ARCH_X86
#   if defined(__GNUC__) && ((__GNUC__>4) || (__GNUC__==4 && __GNUC_MINOR__>=3))
#       define VEC_CPU_ARCH_X86_SSE4_2
#   elif defined(__clang_major__) && defined(__clang_minor__) && ((__clang_major__>3) || (__clang_major__==3 && __clang_minor__>=3))
#       define VEC_CPU_ARCH_X86_SSE4_2
#   endif
#   ifdef VEC_BUILTIN_CPU
#       define VEC_BUILTIN_CPU_SSE4_2
#   endif
#endif
#if defined(__SSE4_2__) && !defined(VEC_CPU_ARCH_X86_SSE4_2)
#   define VEC_CPU_ARCH_X86_SSE4_2
#endif


// AVX
#ifdef VEC_CPU_ARCH_X86
#   if defined(__GNUC__) && ((__GNUC__>4) || (__GNUC__==4 && __GNUC_MINOR__>=4))
#       define VEC_CPU_ARCH_X86_AVX
#   elif defined(__clang_major__) && defined(__clang_minor__) && ((__clang_major__>3) || (__clang_major__==3 && __clang_minor__>=3))
#       define VEC_CPU_ARCH_X86_AVX
#   endif
#   ifdef VEC_BUILTIN_CPU
#       define VEC_BUILTIN_CPU_AVX
#   endif
#endif
#if defined(__AVX__) && !defined(VEC_CPU_ARCH_X86_AVX)
#   define VEC_CPU_ARCH_X86_AVX
#endif


// AVX2
#ifdef VEC_CPU_ARCH_X86
#   if defined(__GNUC__) && ((__GNUC__>4) || (__GNUC__==4 && __GNUC_MINOR__>=7))
#       define VEC_CPU_ARCH_X86_AVX2
#   elif defined(__clang_major__) && defined(__clang_minor__) && ((__clang_major__>3) || (__clang_major__==3 && __clang_minor__>=7))
#       define VEC_CPU_ARCH_X86_AVX2
#   endif
#   ifdef VEC_BUILTIN_CPU
#       define VEC_BUILTIN_CPU_AVX2
#   endif
#endif
#if defined(__AVX2__) && !defined(VEC_CPU_ARCH_X86_AVX2)
#   define VEC_CPU_ARCH_X86_AVX2
#endif


// AVX512F
#ifdef VEC_CPU_ARCH_X86
#   if defined(__GNUC__) && (__GNUC__>=5)
#       define VEC_CPU_ARCH_X86_AVX512F
#   elif defined(__clang_major__) && defined(__clang_minor__) && ((__clang_major__>3) || (__clang_major__==3 && __clang_minor__>=9))
#       define VEC_CPU_ARCH_X86_AVX512F
#   endif
#   if defined(__GNUC__) && (__GNUC__>=6)
#       define VEC_BUILTIN_CPU_AVX512F
#   elif defined(__clang_major__) && defined(__clang_minor__) && ((__clang_major__>3) || (__clang_major__==3 && __clang_minor__>=9))
#       define VEC_BUILTIN_CPU_AVX512F
#   elif defined(__ICC)
#       define VEC_BUILTIN_CPU_AVX512F
#   endif
#endif
#if defined(__AVX512F__) && !defined(VEC_CPU_ARCH_X86_AVX512F)
#   define VEC_CPU_ARCH_X86_AVX512F
#endif


// AVX512BW
#ifdef VEC_CPU_ARCH_X86
#   if defined(__GNUC__) && (__GNUC__>=5)
#       define VEC_CPU_ARCH_X86_AVX512BW
#   elif defined(__clang_major__) && defined(__clang_minor__) && ((__clang_major__>3) || (__clang_major__==3 && __clang_minor__>=9))
#       define VEC_CPU_ARCH_X86_AVX512BW
#   endif
#   if defined(__GNUC__) && (__GNUC__>=6)
#       define VEC_BUILTIN_CPU_AVX512BW
#   elif defined(__clang_major__) && defined(__clang_minor__) && ((__clang_major__>3) || (__clang_major__==3 && __clang_minor__>=9))
#       define VEC_BUILTIN_CPU_AVX512BW
#   elif defined(__ICC)
#       define VEC_BUILTIN_CPU_AVX512BW
#   endif
#endif
#if defined(__AVX512BW__) && !defined(VEC_CPU_ARCH_X86_AVX512BW)
#   define VEC_CPU_ARCH_X86_AVX512BW
#endif


#include <stdint.h>
#include <stdlib.h>


// allows unaligned access
typedef uint32_t unaligned_uint32 __attribute__((aligned(1)));
typedef uint64_t unaligned_uint64 __attribute__((aligned(1)));


#ifdef __cplusplus
extern "C" {
#endif

/// Initialize the target-specific functions
void vec_init_function();


// =========================================================================
// GRM entry calculated from a pair of 2-bit packed genotypes

typedef void (*type_grm_calc_update_f32)(
	const uint8_t p_i[], const uint8_t p_j[], size_t n, const float p_G[],
	int &miss_n, double &sum);
typedef void (*type_grm_calc_update_f64)(
	const uint8_t p_i[], const uint8_t p_j[], size_t n, const double p_G[],
	int &miss_n, double &sum);

extern type_grm_calc_update_f32 fc_grm_calc_update_f32;
extern type_grm_calc_update_f64 fc_grm_calc_update_f64;

void grm_calc_update_f32_def(const uint8_t p_i[], const uint8_t p_j[],
	size_t n, const float p_G[], int &miss_n, double &sum);
void grm_calc_update_f32_avx2(const uint8_t p_i[], const uint8_t p_j[],
	size_t n, const float p_G[], int &miss_n, double &sum);
void grm_calc_update_f32_avx512bw(const uint8_t p_i[], const uint8_t p_j[],
	size_t n, const float p_G[], int &miss_n, double &sum);

void grm_calc_update_f64_def(const uint8_t p_i[], const uint8_t p_j[],
	size_t n, const double p_G[], int &miss_n, double &sum);
void grm_calc_update_f64_avx2(const uint8_t p_i[], const uint8_t p_j[],
	size_t n, const double p_G[], int &miss_n, double &sum);
void grm_calc_update_f64_avx512bw(const uint8_t p_i[], const uint8_t p_j[],
	size_t n, const double p_G[], int &miss_n, double &sum);


// =========================================================================
// Get and set dot using sparse genotypes in get_crossprod_b_grm()

typedef double (*type_get_dot_sp_b)(const double p_std_geno_d[],
	const double p_b[], const uint8_t p_sp_g[]);
typedef void (*type_set_dot_sp_b)(double p_add[], double dot,
	const double p_std_geno_d[], const uint8_t p_sp_g[]);

extern type_get_dot_sp_b fc_get_dot_sp_b;
extern type_set_dot_sp_b fc_set_dot_sp_b;


double get_dot_sp_b_def(const double p_std_geno_d[], const double p_b[],
	const uint8_t p_sp_g[]);
double get_dot_sp_b_avx2(const double p_std_geno_d[], const double p_b[],
	const uint8_t p_sp_g[]);
double get_dot_sp_b_avx512bw(const double p_std_geno_d[], const double p_b[],
	const uint8_t p_sp_g[]);

void set_dot_sp_b_def(double p_add[], double dot, const double p_std_geno_d[],
	const uint8_t p_sp_g[]);
void set_dot_sp_b_avx2(double p_add[], double dot, const double p_std_geno_d[],
	const uint8_t p_sp_g[]);
void set_dot_sp_b_avx512bw(double p_add[], double dot, const double p_std_geno_d[],
	const uint8_t p_sp_g[]);



#ifdef __cplusplus
}
#endif


#endif /* VEC_EXT_H_ */
