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

#if defined(__clang__)
#pragma clang optimize on
#elif defined(__GNUC__) && ((__GNUC__>4) || (__GNUC__==4 && __GNUC_MINOR__>=4))
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
static COREARRAY_TARGET_DEFAULT const char *simd_version() { return "generic"; }
#endif

#ifdef COREARRAY_TARGET_SSE2
static COREARRAY_TARGET_SSE2 const char *simd_version() { return "SSE2"; }
#endif

#ifdef COREARRAY_TARGET_AVX
static COREARRAY_TARGET_AVX const char *simd_version() { return "AVX"; }
#endif

/// SIMD version
extern "C" SEXP saige_simd_version()
{
	const char *s = simd_version();
#ifdef COREARRAY_HAVE_TARGET_CLONES
	char buffer[256];
	stpncpy(buffer, s, sizeof(buffer));
	strcpy(buffer+strlen(s), " (FMV)");
	s = buffer;
#endif
	return mkString(s);
}


// ========================================================================= //
// get the index of each nonzero value in x and return the number of nonzeros

inline static COREARRAY_TARGET_CLONES("avx,sse2,default")
	size_t d_nonzero_index(size_t n, const double *x, int *i)
{
	size_t n_i = 0;
	for (size_t j=0; j < n; j++)
		if (x[j] != 0) i[n_i++] = j;
	return n_i;
}

/// get the index of each nonzero value in x and return the number of nonzeros
extern "C" size_t f64_nonzero_index(size_t n, const double *x, int *i)
{
	return d_nonzero_index(n, x, i);
}


// ========================================================================= //
// y[i] = x - y[i]

inline static COREARRAY_TARGET_CLONES("avx,sse2,default")
	void d_sub(size_t n, double x, double *y)
{
	for (size_t i=0; i < n; i++) y[i] = x - y[i];
}

/// y[i] = x - y[i]
extern "C" void f64_sub(size_t n, double x, double *y)
{
	d_sub(n, x, y);
}


// ========================================================================= //
// y[i] = x * y[i]

inline static COREARRAY_TARGET_CLONES("avx,sse2,default")
	void d_mul(size_t n, double x, double *y)
{
	for (size_t i=0; i < n; i++) y[i] *= x;
}

/// y[i] = x - y[i]
extern "C" void f64_mul(size_t n, double x, double *y)
{
	d_mul(n, x, y);
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
// out1 = sum_i x1[i]*y[i], out2 = sum_i x2[i]*y[i]*y[i]

inline static COREARRAY_TARGET_CLONES("avx,sse2,default")
	void d_dot_sp(size_t n, const double *x1, const double *x2, const double *y,
		double &out1, double &out2)
{
	double sum1=0, sum2=0;
	for (size_t i=0; i < n; i++)
	{
		sum1 += x1[i] * y[i];
		sum2 += x2[i] * y[i] * y[i];
	}
	out1 = sum1; out2 = sum2;
}

/// out1 = sum_i x1[i]*y[i], out2 = sum_i x2[i]*y[i]*y[i]
extern "C" void f64_dot_sp(size_t n, const double *x1, const double *x2,
	const double *y, double &out1, double &out2)
{
	d_dot_sp(n, x1, x2, y, out1, out2);
}


// ========================================================================= //
// vec(p_m) = mat(x_{m*n}) * vec(y_n)
// vec(p_m) = mat(x_{m*n}) * vec(y_n), y is a sparse vector with indices
// vec(p_n) = t(mat(x_{m*n})) * vec(y_m), with a subset

inline static COREARRAY_TARGET_CLONES("avx,sse2,default")
	void d_mul_m_v(size_t n, size_t m, const double *x, const double *y, double *p)
{
	memset(p, 0, sizeof(double)*m);
	for (; n > 0; n--)
	{
		double alpha = *y++;
		if (alpha != 0) // sparse vector
			for (size_t i=0; i < m; i++) p[i] += alpha * x[i];
		x += m;
	}
}

inline static COREARRAY_TARGET_CLONES("avx,sse2,default")
	void d_mul_mat_vec_sp(size_t n, const int *idx, size_t m, const double *x,
		const double *y, double *p)
{
	memset(p, 0, sizeof(double)*m);
	for (; n > 0; n--)
	{
		size_t i = *idx++;
		double alpha = y[i];
		const double *xx = &x[m*i];
		for (size_t j=0; j < m; j++) p[j] += alpha * xx[j];
	}
}

inline static COREARRAY_TARGET_CLONES("avx,sse2,default")
	void d_mul_mat_vec_sub(size_t n, const int *idx, size_t m,
		const double *x, const double *y, double *p)
{
	for (size_t i=0; i < n; i++)
	{
		size_t k = idx[i];
		const double *xx = &x[m*k];
		double sum = 0;
		for (size_t j=0; j < m; j++) sum += y[j] * xx[j];
		p[i] = sum;
	}
}


/// vec(p_m) = mat(x_{m*n}) * vec(y_n), y is a sparse vector
extern "C" void f64_mul_mat_vec(size_t n, size_t m, const double *x,
	const double *y, double *p)
{
	d_mul_m_v(n, m, x, y, p);
}

// vec(p_m) = mat(x_{m*n}) * vec(y_n), y is a sparse vector with indices
extern "C" void f64_mul_mat_vec_sp(size_t n, const int *idx, size_t m,
	const double *x, const double *y, double *p)
{
	d_mul_mat_vec_sp(n, idx, m, x, y, p);
}

/// vec(p_n) = t(mat(x_{m*n})) * vec(y_m), with a subset
extern "C" void f64_mul_mat_vec_sub(size_t n, const int *idx, size_t m,
	const double *x, const double *y, double *p)
{
	d_mul_mat_vec_sub(n, idx, m, x, y, p);
}



// ========================================================================= //
// vec(p_n) = vec(x_n) - t(mat(y_{m*n})) * vec(z_m)

inline static COREARRAY_TARGET_CLONES("avx,sse2,default")
	void d_sub_mul_mat_vec(size_t n, size_t m,
		const double *x, const double *y, const double *z, double *p)
{
	switch (m)
	{
	case 1:
		for (size_t i=0; i < n; i++) p[i] = x[i] - z[0]*y[i];
		break;
	case 2:
		for (size_t i=0; i < n; i++, y+=2)
			p[i] = x[i] - (z[0]*y[0] + z[1]*y[1]);
		break;
	case 3:
		for (size_t i=0; i < n; i++, y+=3)
			p[i] = x[i] - (z[0]*y[0] + z[1]*y[1] + z[2]*y[2]);
		break;
	case 4:
		for (size_t i=0; i < n; i++, y+=4)
			p[i] = x[i] - (z[0]*y[0] + z[1]*y[1] + z[2]*y[2] + z[3]*y[3]);
		break;
	case 5:
		for (size_t i=0; i < n; i++, y+=5)
			p[i] = x[i] - (z[0]*y[0] + z[1]*y[1] + z[2]*y[2] + z[3]*y[3] +
				z[4]*y[4]);
		break;
	case 6:
		for (size_t i=0; i < n; i++, y+=6)
			p[i] = x[i] - (z[0]*y[0] + z[1]*y[1] + z[2]*y[2] + z[3]*y[3] +
				z[4]*y[4] + z[5]*y[5]);
		break;
	case 7:
		for (size_t i=0; i < n; i++, y+=7)
			p[i] = x[i] - (z[0]*y[0] + z[1]*y[1] + z[2]*y[2] + z[3]*y[3] +
				z[4]*y[4] + z[5]*y[5] + z[6]*y[6]);
		break;
	case 8:
		for (size_t i=0; i < n; i++, y+=8)
			p[i] = x[i] - (z[0]*y[0] + z[1]*y[1] + z[2]*y[2] + z[3]*y[3] +
				z[4]*y[4] + z[5]*y[5] + z[6]*y[6] + z[7]*y[7]);
		break;
	default:
		for (; n > 0; n--)
		{
			double sum = 0;
			for (size_t i=0; i < m; i++) sum += y[i] * z[i];
			y += m;
			*p++ = (*x++) - sum;
		}
	}
}

/// vec(p_n) = vec(x_n) - t(mat(y_{m*n})) * vec(z_m)
extern "C" void f64_sub_mul_mat_vec(size_t n, size_t m,
	const double *x, const double *y, const double *z, double *p)
{
	d_sub_mul_mat_vec(n, m, x, y, z, p);
}



// ========================================================================= //
// t(vec(y)) * mat(x) * vec(y)

inline static COREARRAY_TARGET_CLONES("avx,sse2,default")
	double d_sum_mat_vec(size_t n, const double *x, const double *y)
{
	double sum = 0;
	for (size_t i=0; i < n; i++)
	{
		const double *xx = &x[n*i], a = y[i];
		for (size_t j=0; j < n; j++)
			sum += a * y[j] * xx[j];
	}
	return sum;
}

/// t(vec(y)) * mat(x) * vec(y)
extern "C" double f64_sum_mat_vec(size_t n, const double *x, const double *y)
{
	return d_sum_mat_vec(n, x, y);
}
