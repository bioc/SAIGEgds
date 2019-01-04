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


/// sum_i p[i]*s[i]
extern "C" double vec_dot(size_t n, const double *p, const double *s)
{
	int nn=n, inc=1;
	return F77_NAME(ddot)(&nn, p, &inc, s, &inc);
}


/// sum_i p[i]*s[i]*s[i]
extern "C" double vec_dot_sp(size_t n, const double *p, const double *s)
{
	double sum = 0;
	for (; n > 0; n--, p++, s++) sum += (*p) * (*s) * (*s);
	return sum;
}

