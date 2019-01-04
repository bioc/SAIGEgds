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
#include <R_ext/BLAS.h>

// Streaming SIMD Extensions (SSE, SSE2)
#if (defined(__SSE__) && defined(__SSE2__))
#   include <xmmintrin.h>  // SSE
#   include <emmintrin.h>  // SSE2
#endif


/// sum_i x[i]*y[i]
inline static double _dot(size_t n, const double *x, const double *y)
{
	double sum = 0;

#if defined(__SSE__)
	__m128d sum2 = _mm_setzero_pd();
	for (; n >= 2; n-=2)
	{
		__m128d d1 = _mm_loadu_pd(x); x += 2;
		__m128d d2 = _mm_loadu_pd(y); y += 2;
		sum2 = _mm_add_pd(sum2, _mm_mul_pd(d1, d2));
	}
	sum2 = _mm_add_pd(sum2, _mm_shuffle_pd(sum2, sum2, 0x01));
	sum += _mm_cvtsd_f64(sum2);
#endif

	for (; n > 0; n--) sum += (*x++) * (*y++);
	return sum;
}


/// sum_i x[i]*y[i]
extern "C" double vec_dot(size_t n, const double *x, const double *y)
{
	// int nn=n, inc=1;
	// return F77_NAME(ddot)(&nn, x, &inc, y, &inc);
	return _dot(n, x, y);
}


/// sum_i x[i]*y[i]*y[i]
extern "C" double vec_dot_sp(size_t n, const double *x, const double *y)
{
	double sum = 0;
#if defined(__AVX__)

#elif defined(__SSE__)
	__m128d sum2 = _mm_setzero_pd();
	for (; n >= 2; n-=2)
	{
		__m128d d1 = _mm_loadu_pd(x); x += 2;
		__m128d d2 = _mm_loadu_pd(y); y += 2;
		__m128d dd = _mm_mul_pd(d1, _mm_mul_pd(d2, d2));
		sum2 = _mm_add_pd(sum2, dd);
	}
	sum2 = _mm_add_pd(sum2, _mm_shuffle_pd(sum2, sum2, 0x01));
	sum += _mm_cvtsd_f64(sum2);
#endif

	for (; n > 0; n--, x++, y++) sum += (*x) * (*y) * (*y);
	return sum;
}


/// y[i] += alpha * x[i]
inline static void _axpy(size_t n, double alpha, const double *x, double *y)
{
	for (; n > 0; n--) (*y++) += alpha * (*x++);
}





