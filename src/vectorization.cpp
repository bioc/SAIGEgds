// ===========================================================
//
// vectorization.cpp: optimization with vectorization
//
// Copyright (C) 2019    Xiuwen Zheng
//
// This file is part of SAIGEgds.
//
// SAIGEgds is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License Version 3 as
// published by the Free Software Foundation.
//
// SAIGEgds is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public
// License along with SAIGEgds.
// If not, see <http://www.gnu.org/licenses/>.

#include <cstring>
#include <Rdefines.h>
#include <R.h>


// Function multiversioning
#if (defined(__GNUC__) && ((__GNUC__ > 4) || (__GNUC__==4 && __GNUC_MINOR__>=8)))
#   define COREARRAY_HAVE_TARGET
#   define COREARRAY_TARGET(opt)    __attribute__((__target__(opt)))
#   if (__GNUC__ >= 6)
#       define COREARRAY_HAVE_TARGET_CLONES
#       define COREARRAY_TARGET_CLONES(opt)    __attribute__((target_clones(opt)))
#   endif
//#elif defined(__clang__)  // not support
//#   define COREARRAY_HAVE_TARGET
//#   define COREARRAY_TARGET(opt)           __attribute__((target(opt)))
//#   define COREARRAY_TARGET_CLONES(opt)    __attribute__((target_clones(opt)))
#else
#   define COREARRAY_TARGET(opt)
#   define COREARRAY_TARGET_CLONES(opt)
#endif

#ifdef COREARRAY_HAVE_TARGET
#   define COREARRAY_TARGET_DEFAULT    COREARRAY_TARGET("default")
#   define COREARRAY_TARGET_SSE2       COREARRAY_TARGET("sse2")
#   define COREARRAY_TARGET_AVX        COREARRAY_TARGET("avx")
#else
#   if defined(__AVX__)
#       define COREARRAY_TARGET_AVX
#   elif defined(__SSE2__)
#       define COREARRAY_TARGET_SSE2
#   else
#       define COREARRAY_TARGET_DEFAULT
#   endif
#endif


#ifdef COREARRAY_HAVE_TARGET
#   include <xmmintrin.h>  // SSE
#   include <emmintrin.h>  // SSE2
#   include <immintrin.h>  // AVX, AVX2
#else
#   if (defined(__SSE__) && defined(__SSE2__))
#       include <xmmintrin.h>  // SSE
#       include <emmintrin.h>  // SSE2
#   endif
#   if defined(__AVX__) || defined(__AVX2__)
#       include <immintrin.h>  // AVX, AVX2
#   endif
#endif



using namespace std;

// ========================================================================= //

#ifdef COREARRAY_TARGET_DEFAULT
static const char *simd_version() COREARRAY_TARGET_DEFAULT
	{ return "generic"; }
#endif

#ifdef COREARRAY_TARGET_SSE2
static const char *simd_version() COREARRAY_TARGET_SSE2
	{ return "SSE2"; }
#endif

#ifdef COREARRAY_TARGET_AVX
static const char *simd_version() COREARRAY_TARGET_AVX
	{ return "AVX"; }
#endif

/// SIMD version
extern "C" SEXP saige_simd_version()
{
	return mkString(simd_version());
}


// ========================================================================= //
// sum_i x[i]*y[i]

#ifdef COREARRAY_TARGET_DEFAULT
inline static double d_dot(size_t n, const double *x, const double *y)
	COREARRAY_TARGET_DEFAULT
{
	double sum = 0;
	for (; n > 0; n--) sum += (*x++) * (*y++);
	return sum;
}
#endif

#ifdef COREARRAY_TARGET_SSE2
inline static double d_dot(size_t n, const double *x, const double *y)
	COREARRAY_TARGET_SSE2
{
	__m128d sum2 = _mm_setzero_pd();
	for (; n >= 2; n-=2)
	{
		__m128d xx = _mm_loadu_pd(x); x += 2;
		__m128d yy = _mm_loadu_pd(y); y += 2;
		sum2 = _mm_add_pd(sum2, _mm_mul_pd(xx, yy));
	}
	sum2 = _mm_add_pd(sum2, _mm_shuffle_pd(sum2, sum2, 0x01));
	double sum = _mm_cvtsd_f64(sum2);
	if (n > 0) sum += (*x) * (*y);
	return sum;
}
#endif

#ifdef COREARRAY_TARGET_AVX
inline static double d_dot(size_t n, const double *x, const double *y)
	COREARRAY_TARGET_AVX
{
	// AVX
	__m256d sum4 = _mm256_setzero_pd();
	for (; n >= 4; n-=4)
	{
		__m256d xx = _mm256_loadu_pd(x); x += 4;
		__m256d yy = _mm256_loadu_pd(y); y += 4;
		sum4 = _mm256_add_pd(sum4, _mm256_mul_pd(xx, yy));
	}
	sum4 = _mm256_add_pd(sum4, _mm256_permute2f128_pd(sum4, sum4, 0x01));
	// SSE
	__m128d sum2 = _mm256_castpd256_pd128(sum4);
	if (n >= 2)
	{
		__m128d xx = _mm_loadu_pd(x); x += 2;
		__m128d yy = _mm_loadu_pd(y); y += 2;
		sum2 = _mm_add_pd(sum2, _mm_mul_pd(xx, yy));
		n -= 2;
	}
	sum2 = _mm_add_pd(sum2, _mm_shuffle_pd(sum2, sum2, 0x01));
	// remainder
	double sum = _mm_cvtsd_f64(sum2);
	if (n > 0) sum += (*x) * (*y);
	return sum;
}
#endif

/// sum_i x[i]*y[i]
extern "C" double f64_dot(size_t n, const double *x, const double *y)
{
	return d_dot(n, x, y);
}


// ========================================================================= //
// sum_i x[i]*y[i]*y[i]

#ifdef COREARRAY_TARGET_DEFAULT
inline static double d_dot_sp(size_t n, const double *x, const double *y)
	COREARRAY_TARGET_DEFAULT
{
	double sum = 0;
	for (; n > 0; n--) sum += (*x++) * (*y++);
	return sum;
}
#endif

#ifdef COREARRAY_TARGET_SSE2
inline static double d_dot_sp(size_t n, const double *x, const double *y)
	COREARRAY_TARGET_SSE2
{
	__m128d sum2 = _mm_setzero_pd();
	for (; n >= 2; n-=2)
	{
		__m128d xx = _mm_loadu_pd(x); x += 2;
		__m128d yy = _mm_loadu_pd(y); y += 2;
		sum2 = _mm_add_pd(sum2, _mm_mul_pd(xx, _mm_mul_pd(yy, yy)));
	}
	sum2 = _mm_add_pd(sum2, _mm_shuffle_pd(sum2, sum2, 0x01));
	double sum = _mm_cvtsd_f64(sum2);
	if (n > 0) sum += (*x) * (*y) * (*y);
	return sum;
}
#endif

#ifdef COREARRAY_TARGET_AVX
inline static double d_dot_sp(size_t n, const double *x, const double *y)
	COREARRAY_TARGET_AVX
{
	// AVX
	__m256d sum4 = _mm256_setzero_pd();
	for (; n >= 4; n-=4)
	{
		__m256d xx = _mm256_loadu_pd(x); x += 4;
		__m256d yy = _mm256_loadu_pd(y); y += 4;
		sum4 = _mm256_add_pd(sum4, _mm256_mul_pd(xx, _mm256_mul_pd(yy, yy)));
	}
	sum4 = _mm256_add_pd(sum4, _mm256_permute2f128_pd(sum4, sum4, 0x01));
	// SSE
	__m128d sum2 = _mm256_castpd256_pd128(sum4);
	if (n >= 2)
	{
		__m128d xx = _mm_loadu_pd(x); x += 2;
		__m128d yy = _mm_loadu_pd(y); y += 2;
		sum2 = _mm_add_pd(sum2, _mm_mul_pd(xx, _mm_mul_pd(yy, yy)));
	}
	sum2 = _mm_add_pd(sum2, _mm_shuffle_pd(sum2, sum2, 0x01));
	// remainder
	double sum = _mm_cvtsd_f64(sum2);
	if (n > 0) sum += (*x) * (*y) * (*y);
	return sum;
}
#endif

/// sum_i x[i]*y[i]*y[i]
extern "C" double f64_dot_sp(size_t n, const double *x, const double *y)
{
	return d_dot_sp(n, x, y);
}


// ========================================================================= //

/// y[i] += alpha * x[i]
inline static void _axpy(size_t n, double alpha, const double *x, double *y)
{
	for (; n > 0; n--) (*y++) += alpha * (*x++);
}





