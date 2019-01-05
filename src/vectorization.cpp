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

#ifdef __GNUC__
#pragma GCC optimize("O3")
#endif

// Function multiversioning (requiring target_clones)
#if (defined(__GNUC__) && (__GNUC__ >= 6))
#   define COREARRAY_HAVE_TARGET
#   define COREARRAY_TARGET(opt)    __attribute__((target(opt)))
#   define COREARRAY_HAVE_TARGET_CLONES
#   define COREARRAY_TARGET_CLONES(opt)    __attribute__((target_clones(opt)))
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


#include <cstring>
#include <Rdefines.h>
#include <R.h>


using namespace std;


// ========================================================================= //

#ifdef COREARRAY_TARGET_DEFAULT
static COREARRAY_TARGET_DEFAULT const char *simd_version()
	{ return "generic"; }
#endif

#ifdef COREARRAY_TARGET_SSE2
static COREARRAY_TARGET_SSE2 const char *simd_version()
	{ return "SSE2"; }
#endif

#ifdef COREARRAY_TARGET_AVX
static COREARRAY_TARGET_AVX const char *simd_version()
	{ return "AVX"; }
#endif

/// SIMD version
extern "C" SEXP saige_simd_version()
{
	return mkString(simd_version());
}


// ========================================================================= //
// sum_i x[i]*y[i]

inline static COREARRAY_TARGET_CLONES("avx,sse2,default")
	double d_dot(size_t n, const double *x, const double *y)
{
	double sum = 0;
	for (size_t i=0; i < n; i++) sum += x[i] * y[i];
	return sum;
}

/// sum_i x[i]*y[i]
extern "C" double f64_dot(size_t n, const double *x, const double *y)
{
	return d_dot(n, x, y);
}


// ========================================================================= //
// sum_i x[i]*y[i]*y[i]

inline static COREARRAY_TARGET_CLONES("avx,sse2,default")
	double d_dot_sp(size_t n, const double *x, const double *y)
{
	double sum = 0;
	for (size_t i=0; i < n; i++) sum += x[i] * y[i] * y[i];
	return sum;
}

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





