// ===========================================================
//
// vectorization.h: optimization with vectorization
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

#include <string.h>

extern "C"
{
	/// return allele frequency and impute genotype using the mean
	void f64_af_ac_impute(double *ds, size_t n, double &AF, double &AC,
		int &Num, int buf_idx[]);
	/// get the index of each nonzero value in x and return the number of nonzeros
	size_t f64_nonzero_index(size_t n, const double *x, int *i);
	/// y[i] = x - y[i]
	void f64_sub(size_t n, double x, double *y);
	/// y[i] = x * y[i]
	void f64_mul(size_t n, double x, double *y);
	/// sum_i x[i]*y[i]
	double f64_dot(size_t n, const double *x, const double *y);
	/// out1 = sum_i x1[i]*y[i], out2 = sum_i x2[i]*y[i]*y[i]
	void f64_dot_sp(size_t n, const double *x1, const double *x2,
		const double *y, double &out1, double &out2);
	/// vec(p_m) = mat(x_{m*n}) * vec(y_n), y is a sparse vector
	void f64_mul_mat_vec(size_t n, size_t m, const double *x,
		const double *y, double *p);
	/// vec(p_m) = mat(x_{m*n}) * vec(y_n), y is a sparse vector with indices
	void f64_mul_mat_vec_sp(size_t n, const int *idx, size_t m,
		const double *x, const double *y, double *p);
	/// vec(p_n) = t(mat(x_{m*n})) * vec(y_m), with a subset
	void f64_mul_mat_vec_sub(size_t n, const int *idx, size_t m,
		const double *x, const double *y, double *p);
	/// vec(p_n) = vec(x_n) - t(mat(y_{m*n})) * vec(z_m)
	void f64_sub_mul_mat_vec(size_t n, size_t m,
		const double *x, const double *y, const double *z, double *p);
	/// t(vec(y)) * mat(x) * vec(y)
	double f64_sum_mat_vec(size_t n, const double *x, const double *y);
}

