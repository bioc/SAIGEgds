// ===========================================================
//
// saige_fitnull.cpp: C++ implementation of fitting the null model
//
// Copyright (C) 2019-2022    Xiuwen Zheng / AbbVie-ComputationalGenomics
//
// This file is part of SAIGEgds. It was created based on the original SAIGE
// C++ and R codes in the SAIGE package. Compared with the original SAIGE,
// all single-precision floating-point numbers are changed to double precision,
// and a more efficient algorithm with packed 2-bit and sparse genotypes is
// implemented to calculate the cross product of genetic relationship matrix.
//
// SAIGEgds is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License Version 3 as published
// by the Free Software Foundation.
//
// SAIGEgds is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along
// with SAIGEgds.
// If not, see <http://www.gnu.org/licenses/>.

#include "vectorization.h"
#include "vec_ext.h"
#include <RcppArmadillo.h>
#include "saige.h"
#include <climits>
#include <vector>
#include <algorithm>

using namespace std;
using namespace Rcpp;
using namespace arma;
using namespace RcppParallel;
using namespace vectorization;
using namespace SAIGE;


// ========================================================================= //
/// SPAtest

extern "C" double Saddle_Prob(double q, double m1, double var1, size_t n_g,
	const double mu[], const double g[], double cutoff, bool &converged, double *p_noadj);

extern "C" double Saddle_Prob_Fast(double q, double m1, double var1, size_t n_g,
	const double mu[], const double g[], size_t n_nonzero, const int nonzero_idx[],
	double cutoff, bool &converged, double buf_spa[], double *p_noadj);


// ========================================================================= //
// Internal functions

inline static double ds_nan(BYTE g)
{
	return (g < 3) ? g : R_NaN;
}

inline static NumericVector SEXP_VEC(const dvec &x)
{
	return NumericVector(x.begin(), x.end());
}

inline static void set_seed(unsigned int seed)
{
	Environment base_env("package:base");
	Function set_seed_r = base_env["set.seed"];
	set_seed_r(seed);
}

/// Print a numeric vector
inline static void print_vec(const char *s, const dvec &x,
	const char *indent=NULL, bool nl=true)
{
	if (!indent) indent = "";
	Rprintf("%s%s(", indent, s);
	for (size_t i=0; i < x.n_elem; i++)
	{
		if (i > 0) Rprintf(", ");
		Rprintf("%0.7g", x[i]);
	}
	if (nl) Rprintf(")\n"); else Rprintf(")");
}


// ========================================================================= //

RcppExport SEXP saige_set_numthread(SEXP nthread)
{
	const int n = Rf_asInteger(nthread);
	if (n > 0)
		SAIGE_NumThread = n;
	else
		Rf_error("Invalid number of threads: %d.", n);
	return R_NilValue;
}


// ========================================================================= //
// GRM data structure

struct TypeGRM
{
	int NumSamp;        //< the number of samples
	int NumVariant;     //< the number of variants in GRM
	int PackedNumSamp;  //< the number of bytes for 2-bit packed samples

	// ====  Full GRM  ====
	/// pointer to the 2-bit packed genotypes for full GRM
	BYTE *PackedG;
	/// list vector of sparse structure of genotypes for full GRM
	/// integers: c(n1, n2, n3, n1-length int vector, n2-length, n3-length)
	/// n1 for genotype 1, n2 for genotype 2, n3 for missing genotype
	SEXP SparseG;
	// buffer for full GRM
	double *buf_std_geno;   //< a 4-by-n_variant look-up matrix
	double *buf_crossprod;  //< NumThread-by-n_samp matrix
	double *buf_diag_grm;   //< n_samp-length, sigma_i = avg_j adj.g[i,j]^2

	// ====  Approximate dense or sparse GRM  ====
	/// user-defined approximate dense matrix
	Type_Matrix ApproxDMat;
	/// user-defined approximate sparse matrix
	Type_dgCMatrix ApproxSpMat;

	TypeGRM() { Reset(); }
	void Reset()
	{
		NumSamp = NumVariant = PackedNumSamp = 0;
		PackedG = NULL; SparseG = NULL;
		buf_std_geno = buf_crossprod = buf_diag_grm = NULL;
		ApproxDMat = NULL;
		ApproxSpMat.reset(NULL);
	}
};

static TypeGRM GRM;


/// Initialize GRM storage
RcppExport SEXP saige_init_fit_grm()
{
	GRM.Reset();
	return R_NilValue;
}


/// List for the 2-bit packed genotypes in random SNP markers for variance ratio
static SEXP Geno_PackedRawRandVR = NULL;


// ========================================================================= //
// Store 2-bit packed SNP genotypes

// Internal 2-bit genotype lookup tables
static bool lookup_has_init = false;
static BYTE num_valid[256], num_sum[256];

static void init_lookup_table()
{
	if (lookup_has_init) return;
	for (int i=0; i < 256; i++)
	{
		int b0 = i & 0x03;
		int b1 = (i >> 2) & 0x03;
		int b2 = (i >> 4) & 0x03;
		int b3 = (i >> 6) & 0x03;
		num_valid[i] = (b0<3 ? 1:0) + (b1<3 ? 1:0) + (b2<3 ? 1:0) + (b3<3 ? 1:0);
		num_sum[i] = (b0<3 ? b0:0) + (b1<3 ? b1:0) + (b2<3 ? b2:0) + (b3<3 ? b3:0);
	}
	lookup_has_init = true;
}

/// Store SNP genotype in the 2-bit packed format
RcppExport SEXP saige_store_2b_geno(SEXP rawgeno, SEXP raw_g_RandVR,
	SEXP num_samp, SEXP r_buf_geno, SEXP r_buf_sigma, SEXP r_buf_crossprod)
{
BEGIN_RCPP

	// initialize the basic genotype variables
	RawMatrix RawG(rawgeno);
	GRM.PackedG = RAW(rawgeno);  // use packed 2-bit genotypes
	GRM.SparseG = NULL;          // not use sparse genotype list
	GRM.NumSamp = Rf_asInteger(num_samp);
	GRM.PackedNumSamp = RawG.nrow();
	GRM.NumVariant = RawG.ncol();

	// genotype matrix list for estimating variance ratios
	Geno_PackedRawRandVR = raw_g_RandVR;
	// check
	for (int i=0; i < Rf_length(raw_g_RandVR); i++)
	{
		RawMatrix G(VECTOR_ELT(raw_g_RandVR, i));
		if (RawG.nrow() != G.nrow())
			throw std::invalid_argument("Invalid argument in saige_store_2b_geno()");
	}

	// set the buffer for get_crossprod_b_grm()
	NumericMatrix mat(r_buf_crossprod);
	GRM.buf_crossprod = REAL(r_buf_crossprod);
	if (SAIGE_NumThread > GRM.NumSamp)
		SAIGE_NumThread = GRM.NumSamp;
	if (SAIGE_NumThread > GRM.NumVariant)
		SAIGE_NumThread = GRM.NumVariant;
	if (SAIGE_NumThread < 1) SAIGE_NumThread = 1;

	// build the look-up table of standardized genotypes
	init_lookup_table();
	GRM.buf_std_geno = REAL(r_buf_geno);

	PARALLEL_THREAD_BLOCK
	PARALLEL_FOR(i, GRM.NumVariant, true)
	{
		BYTE *g = GRM.PackedG + GRM.PackedNumSamp*i;
		// calculate allele frequency
		int n_valid=0, sum=0;
		for (int j=0; j < GRM.PackedNumSamp; j++)
		{
			n_valid += num_valid[g[j]];
			sum += num_sum[g[j]];
		}
		double af = double(sum) / (2*n_valid);
		double inv = 1 / sqrt(2*af*(1-af));
		if (!R_FINITE(af) || !R_FINITE(inv))
			af = inv = 0;
		double *p = &GRM.buf_std_geno[4*i];
		p[0] = (0 - 2*af) * inv; p[1] = (1 - 2*af) * inv;
		p[2] = (2 - 2*af) * inv; p[3] = 0;
	}
	PARALLEL_END
	PARALLEL_THREAD_BLOCK_END

	// calculate diag(grm)
	GRM.buf_diag_grm = REAL(r_buf_sigma);
	memset(GRM.buf_diag_grm, 0, sizeof(double)*GRM.NumSamp);
	for (int i=0; i < GRM.NumVariant; i++)
	{
		BYTE *g = GRM.PackedG + GRM.PackedNumSamp*i;
		const double *base = GRM.buf_std_geno + 4*i;
		size_t n = GRM.NumSamp;
		double *p = GRM.buf_diag_grm;
		for (; n >= 4; n-=4, p+=4)
		{
			BYTE gg = *g++;
			p[0] += sq(base[gg & 0x03]);
			p[1] += sq(base[(gg >> 2) & 0x03]);
			p[2] += sq(base[(gg >> 4) & 0x03]);
			p[3] += sq(base[gg >> 6]);
		}
		for (BYTE gg = (n>0 ? *g : 0); n > 0; n--)
		{
			(*p++) += sq(base[gg & 0x03]);
			gg >>= 2;
		}
	}
	f64_mul(GRM.NumSamp, 1.0 / GRM.NumVariant, GRM.buf_diag_grm);

END_RCPP
}


// ========================================================================= //
// Store sparse structure of SNP genotypes

/// sparse genotype buffer
static int *tmp_sp_geno_b256_buf = NULL;
static BYTE *tmp_sp_geno_b1_buf = NULL;

RcppExport SEXP saige_init_sparse(SEXP nsamp, SEXP buf_b256, SEXP buf_b1)
{
	GRM.NumSamp = Rf_asInteger(nsamp);
	tmp_sp_geno_b256_buf = INTEGER(buf_b256);
	tmp_sp_geno_b1_buf = (BYTE*)RAW(buf_b1);
	return R_NilValue;
}

RcppExport SEXP saige_get_sparse_info(SEXP x)
{
	const size_t n = Rf_xlength(x);
	const BYTE *s = (const BYTE *)RAW(x);
	// get #
	int m = 0;
	for (size_t p = *((unaligned_uint32 *)s); p < n; p+=s[p] + 2, m++);
	// fill
	SEXP rv = NEW_INTEGER(m);
	int *out = INTEGER(rv);
	m = 0;
	for (size_t p = *((unaligned_uint32 *)s); p < n; )
	{
		out[m++] = s[p] + 1;
		p += s[p] + 2;
	}
	// output
	return rv;
}

/// convert to BYTE(0, 1, 2, ...) if needed
static void *get_geno_byte(SEXP geno)
{
	switch (TYPEOF(geno))
	{
	case RAWSXP:
		return RAW(geno);
	case INTSXP:
		{
			const size_t num = Rf_xlength(geno);
			void *base_ptr = INTEGER(geno);
			BYTE *p = (BYTE*)base_ptr;
			const int *s = (int*)base_ptr;
			for (size_t i=0; i < num; i++)
				if (0<=s[i] && s[i]<=2) p[i]=s[i]; else p[i]=3;
			return base_ptr;
		}
		break;
	case REALSXP:
		{
			const size_t num = Rf_xlength(geno);
			void *base_ptr = REAL(geno);
			BYTE *p = (BYTE*)base_ptr;
			const double *s = (double*)base_ptr;
			for (size_t i=0; i < num; i++)
			{
				if (R_FINITE(s[i]))
				{
					int g = (int)round(s[i]);
					if (0<=g && g<=2) p[i]=g; else p[i]=3;
				} else
					p[i] = 3;
			}
			return base_ptr;
		}
		break;
	default:
		throw std::invalid_argument("Invalid data type.");
	}
}

static int find_nonzero(int i, int num, const BYTE *s)
{
	for (int n=num-i; n >= 8; n-=8, i+=8)
		if (*((const unaligned_uint64*)&s[i]) != 0) break;
	for (; (i < num) && (s[i] == 0); ) i++;
	return i;
}

static BYTE* get_geno_sp_idx(const BYTE *p_g, BYTE geno, size_t n, BYTE *p_idx)
{
	// n should be <= 256
	for (size_t i=0; i < n; i++)
		if (p_g[i] == geno) *p_idx++ = i;
	return p_idx;
}

static BYTE* get_geno_sp_idx2(const BYTE *p_g, BYTE geno, size_t n, BYTE *p_idx)
{
	// n should be <= 256
	for (size_t i=0; i < n; i++)
		if (p_g[i] >= geno) *p_idx++ = i;
	return p_idx;
}

/// Get sparse structure from the input genotypes
RcppExport SEXP saige_get_sparse(SEXP geno)
{
BEGIN_RCPP
	const int num = GRM.NumSamp;
	if (num > Rf_length(geno))
		throw std::invalid_argument("No enough genotypes.");

	// convert to BYTE(0, 1, 2, ...) if needed
	void *base_ptr = get_geno_byte(geno);
	// determine whether need to flip
	BYTE *gs = (BYTE*)base_ptr;
	int n=0, sum=0;
	for (int i=0; i < num; i++)
		if (gs[i] < 3) { sum += gs[i]; n++; }
	if (sum > n)  // flip
	{
		for (int i=0; i < num; i++)
			if (gs[i] < 3) gs[i] = 2 - gs[i];
	}

	// construct sparse structure
	int m, *p_b256 = tmp_sp_geno_b256_buf + 1;
	BYTE *p = tmp_sp_geno_b1_buf;
	// genotype 1
	{
		int &n256 = p_b256[0]; n256 = 0; p_b256 ++;
		for (int i=0; i < num; i+=256)
		{
			i = find_nonzero(i, num, gs);  // find the first non-zero
			int n = std::min(256, num-i);  // n should be <= 256
			BYTE *p_next = get_geno_sp_idx(&gs[i], 1, n, p+1);
			if ((m = p_next - (p+1)) > 0)  // # of nonzero
			{
				*p = m - 1;  // the first is # of nonzero - 1
				p = p_next; n256 ++; *p_b256++ = i;
			}
		}
	}
	// genotype 2
	{
		int &n256 = p_b256[0]; n256 = 0; p_b256 ++;
		for (int i=0; i < num; i+=256)
		{
			i = find_nonzero(i, num, gs);  // find the first non-zero
			int n = std::min(256, num-i);  // n should be <= 256
			BYTE *p_next = get_geno_sp_idx(&gs[i], 2, n, p+1);
			if ((m = p_next - (p+1)) > 0)  // # of nonzero
			{
				*p = m - 1;  // the first is # of nonzero - 1
				p = p_next; n256 ++; *p_b256++ = i;
			}
		}
	}
	// genotype 3: missing value
	{
		int &n256 = p_b256[0]; n256 = 0; p_b256 ++;
		for (int i=0; i < num; i+=256)
		{
			i = find_nonzero(i, num, gs);  // find the first non-zero
			int n = std::min(256, num-i);  // n should be <= 256
			BYTE *p_next = get_geno_sp_idx2(&gs[i], 3, n, p+1);
			if ((m = p_next - (p+1)) > 0)  // # of nonzero
			{
				*p = m - 1;  // the first is # of nonzero - 1
				p = p_next; n256 ++; *p_b256++ = i;
			}
		}
	}

	// finally
	int &nb_256 = tmp_sp_geno_b256_buf[0];  // in bytes
	nb_256 = (p_b256 - tmp_sp_geno_b256_buf) * 4;
	int nb_b1 = p - tmp_sp_geno_b1_buf;

	// output
	SEXP rv_ans = NEW_RAW(nb_256 + nb_b1);
	PROTECT(rv_ans);
	memcpy(RAW(rv_ans), tmp_sp_geno_b256_buf, nb_256);
	memcpy(RAW(rv_ans)+nb_256, tmp_sp_geno_b1_buf, nb_b1);
	UNPROTECT(1);
	return rv_ans;
END_RCPP
}


/// Initialize sparse structure of genotypes
RcppExport SEXP saige_store_sp_geno(SEXP sp_geno_list, SEXP raw_g_RandVR,
	SEXP num_samp, SEXP r_buf_geno, SEXP r_buf_sigma, SEXP r_buf_crossprod)
{
BEGIN_RCPP

	// initialize the basic genotype variables
	GRM.SparseG = sp_geno_list;  // use sparse genotype list
	GRM.PackedG = NULL;          // not use packed 2-bit genotypes
	GRM.NumSamp = Rf_asInteger(num_samp);
	GRM.NumVariant = Rf_length(sp_geno_list);

	// genotype matrix list for estimating variance ratios
	Geno_PackedRawRandVR = raw_g_RandVR;
	RawMatrix RawG(VECTOR_ELT(raw_g_RandVR, 0));
	GRM.PackedNumSamp = RawG.nrow();
	// check
	for (int i=1; i < Rf_length(raw_g_RandVR); i++)
	{
		RawMatrix G(VECTOR_ELT(raw_g_RandVR, i));
		if (RawG.nrow() != G.nrow())
			throw std::invalid_argument("Invalid argument in saige_store_sp_geno()");
	}

	// set the buffer for get_crossprod_b_grm()
	NumericMatrix mat(r_buf_crossprod);
	GRM.buf_crossprod = REAL(r_buf_crossprod);
	if (SAIGE_NumThread > GRM.NumSamp)
		SAIGE_NumThread = GRM.NumSamp;
	if (SAIGE_NumThread > GRM.NumVariant)
		SAIGE_NumThread = GRM.NumVariant;
	if (SAIGE_NumThread < 1) SAIGE_NumThread = 1;

	// build the look-up table of standardized genotypes
	GRM.buf_std_geno = REAL(r_buf_geno);

	PARALLEL_THREAD_BLOCK
	PARALLEL_FOR(i, GRM.NumVariant, true)
	{
		const BYTE *pg = (const BYTE*)RAW(VECTOR_ELT(GRM.SparseG, i));
		const int *s256 = (const int *)pg;
		const BYTE *s = pg + (*s256++);
		int n[3] = { 0, 0, 0 };
		for (int k=0; k < 3; k++)
		{
			int n256 = *s256++;
			for (int j=0; j < n256; j++)
			{
				const int m = s[0] + 1;  // # of nonzero
				n[k] += m; s += m + 1;
			}
			s256 += n256;
		}
		const int n_valid = GRM.NumSamp - n[2];
		const int g_sum = n[0] + 2*n[1];
		double af = double(g_sum) / (2*n_valid);
		double inv = 1 / sqrt(2*af*(1-af));
		if (!R_FINITE(af) || !R_FINITE(inv))
			af = inv = 0;
		double *p = &GRM.buf_std_geno[4*i];
		p[0] = (0 - 2*af) * inv; p[1] = (1 - 2*af) * inv;
		p[2] = (2 - 2*af) * inv; p[3] = 0;
		p[1] -= p[0]; p[2] -= p[0]; p[3] -= p[0];  // adjustment
	}
	PARALLEL_END
	PARALLEL_THREAD_BLOCK_END

	// calculate diag(grm)
	GRM.buf_diag_grm = REAL(r_buf_sigma);
	memset(GRM.buf_diag_grm, 0, sizeof(double)*GRM.NumSamp);
	double adj_g0 = 0;
	for (int i=0; i < GRM.NumVariant; i++)
	{
		const double *p = &GRM.buf_std_geno[4*i];
		const BYTE *pg = (const BYTE*)RAW(VECTOR_ELT(GRM.SparseG, i));
		const int *s256 = (const int *)pg;
		const BYTE *s = pg + (*s256++);
		// g0^2
		const double g0_2 = sq(p[0]);
		adj_g0 += g0_2;
		// g1^2, g2^2, g3^2
		for (int k=0; k < 3; k++)
		{
			const double v = sq(p[k+1] + p[0]) - g0_2;
			for (int n256 = *s256++; n256 > 0; n256--)
			{
				const size_t bb = *s256++;
				for (int n = *s++; n >= 0; n--)  // for-each nonzero
					GRM.buf_diag_grm[bb + (*s++)] += v;
			}
		}
	}
	f64_add(GRM.NumSamp, adj_g0, GRM.buf_diag_grm);
	f64_mul(GRM.NumSamp, 1.0 / GRM.NumVariant, GRM.buf_diag_grm);

END_RCPP
}


// ========================================================================= //
// Store used-defined approximate dense or sparse GRM

/// Initialize approximate dense GRM matrix
RcppExport SEXP saige_store_dense_grm(SEXP num_samp, SEXP d_mat,
	SEXP r_buf_sigma)
{
BEGIN_RCPP
	GRM.NumSamp = Rf_asInteger(num_samp);
	if (SAIGE_NumThread > GRM.NumSamp)
		SAIGE_NumThread = GRM.NumSamp;

	GRM.ApproxDMat.reset(d_mat);  // use dense matrix
	Type_Matrix &M = GRM.ApproxDMat;
	if (GRM.NumSamp!=M.ncol() || GRM.NumSamp!=M.nrow())
		throw std::invalid_argument("Invalid GRM in saige_store_dense_grm().");
	GRM.ApproxSpMat.reset(NULL);  // not use sparse matrix

	if (!GRM.PackedG && !GRM.SparseG)
	{
		// no full GRM
		GRM.buf_diag_grm = REAL(r_buf_sigma);  // diagonal of GRM
		for (size_t i=0; i < (size_t)GRM.NumSamp; i++)
			GRM.buf_diag_grm[i] = GRM.ApproxDMat.val[i + i*GRM.NumSamp];
	}
END_RCPP
}


/// Initialize approximate sparse GRM matrix
RcppExport SEXP saige_store_sparse_grm(SEXP num_samp, SEXP sp_mat,
	SEXP r_buf_sigma)
{
BEGIN_RCPP
	GRM.NumSamp = Rf_asInteger(num_samp);
	if (SAIGE_NumThread > GRM.NumSamp)
		SAIGE_NumThread = GRM.NumSamp;

	GRM.ApproxSpMat.reset(sp_mat);  // use sparse matrix
	GRM.ApproxDMat.reset(NULL);     // not use dense matrix

	if (!GRM.PackedG && !GRM.SparseG)
	{
		// no full GRM
		GRM.buf_diag_grm = REAL(r_buf_sigma);  // diagonal of GRM
		for (int i=0; i < GRM.NumSamp; i++)
		{
			double d = 0;
			const int j_end = GRM.ApproxSpMat.p[i+1];
			for (int j=GRM.ApproxSpMat.p[i]; j < j_end; j++)
			{
				if (GRM.ApproxSpMat.i[j] == i)
					{ d = GRM.ApproxSpMat.x[j]; break; }
			}
			GRM.buf_diag_grm[i] = d;
		}
	}
END_RCPP
}


// ========================================================================= //

/// Cross-product of standardized genotypes (G) and a numeric vector
/// Input: b (n_samp-length)
/// Output: out_b (n_samp-length) = GRM * b = G' G b
static MATH_OFAST void get_crossprod_b_grm(const dcolvec &b, dvec &out_b)
{
	// dense or sparse genotypes
	if (GRM.SparseG || GRM.PackedG)
	{
		// initialize
		memset(GRM.buf_crossprod, 0, sizeof(double)*GRM.NumSamp*SAIGE_NumThread);
		double sum_b = f64_sum(GRM.NumSamp, &b[0]);
		dvec sum_cp_g0;
		sum_cp_g0.zeros(SAIGE_NumThread);

		// crossprod with b
		if (GRM.SparseG)
		{
			// sparse genotypes
			PARALLEL_FOR(i, GRM.NumVariant, true)
			{
				const double *p = &GRM.buf_std_geno[4*i];
				const BYTE *pg = (const BYTE*)RAW(VECTOR_ELT(GRM.SparseG, i));
				// calculate dot = g * b
				double d =
					sum_b * p[0] +                     // g0 * b
					(*fc_get_dot_sp_b)(p, &b[0], pg);  // g1 * b, g2 * b, g3 * b
				// update GRM.buf_crossprod += d .* std.geno
				double *pbb = GRM.buf_crossprod + GRM.NumSamp * th_idx;
				sum_cp_g0[th_idx] += d * p[0];      // g0 * d
				(*fc_set_dot_sp_b)(pbb, d, p, pg);  // g1 * d, g2 * d, g3 * d
			}
			PARALLEL_END
		} else {
			// dense packed genotypes
			PARALLEL_FOR(i, GRM.NumVariant, true)
			{
				const BYTE *g = GRM.PackedG + GRM.PackedNumSamp*i;
				const double *base = GRM.buf_std_geno + 4*i, *pb = &b[0];

				// get dot = sum(std.geno .* b)
				double dot = 0;
				size_t n = GRM.NumSamp;
				for (; n >= 4; n-=4, pb+=4)
				{
					BYTE gg = *g++;
					dot += base[gg & 0x03] * pb[0] + base[(gg >> 2) & 0x03] * pb[1] +
						base[(gg >> 4) & 0x03] * pb[2] + base[gg >> 6] * pb[3];
				}
				for (BYTE gg = (n>0 ? *g : 0); n > 0; n--)
				{
					dot += base[gg & 0x03] * (*pb++);
					gg >>= 2;
				}

				// update GRM.buf_crossprod += dot .* std.geno
				double *pbb = GRM.buf_crossprod + GRM.NumSamp * th_idx;
				g = GRM.PackedG + GRM.PackedNumSamp*i;
				n = GRM.NumSamp;
				for (; n >= 4; n-=4, pbb+=4)
				{
					BYTE gg = *g++;
					pbb[0] += dot * base[gg & 0x03];
					pbb[1] += dot * base[(gg >> 2) & 0x03];
					pbb[2] += dot * base[(gg >> 4) & 0x03];
					pbb[3] += dot * base[gg >> 6];
				}
				for (BYTE gg = (n>0 ? *g : 0); n > 0; n--)
				{
					(*pbb++) += dot * base[gg & 0x03];
					gg >>= 2;
				}
			}
			PARALLEL_END
		}

		// normalize out_b
		out_b.resize(GRM.NumSamp);
		const double sum_g0 = sum(sum_cp_g0);
		const double scalar = 1.0/GRM.NumVariant;
		PARALLEL_RANGE(st, ed, GRM.NumSamp, false)
		{
			size_t len = ed - st;
			const double *s = GRM.buf_crossprod + st;
			double *p = &out_b[st];
			memset(p, 0, sizeof(double)*len);
			for (int i=0; i < SAIGE_NumThread; i++, s += GRM.NumSamp)
				f64_add(len, s, p);
			if (!GRM.PackedG)
				f64_add(len, sum_g0, p);
			f64_mul(len, scalar, p);
		}
		PARALLEL_END

	} else if (!GRM.ApproxDMat.empty())
	{
		// user-defined dense matrix of GRM
		out_b.resize(GRM.NumSamp);
		PARALLEL_FOR(i, GRM.NumSamp, true)
		{
			const double *pb = &b[0];
			const double *p = &GRM.ApproxDMat.val[GRM.NumSamp*i];
			double sum = 0;
			for (int j=0; j < GRM.NumSamp; j++)
				sum += p[j] * pb[j];
			out_b[i] = sum;
		}
		PARALLEL_END
	
	} else if (!GRM.ApproxSpMat.empty())
	{
		// user-defined sparse matrix of GRM, dgCMatrix
		out_b.resize(GRM.NumSamp);
		PARALLEL_FOR(i, GRM.NumSamp, true)
		{
			const double *pb = &b[0];
			double sum = 0;
			const int C_end = GRM.ApproxSpMat.p[i+1];  // dgCMatrix
			for (int k=GRM.ApproxSpMat.p[i]; k < C_end; k++)
				sum += pb[GRM.ApproxSpMat.i[k]] * GRM.ApproxSpMat.x[k];
			// output
			out_b[i] = sum;
		}
		PARALLEL_END

	} else
		throw std::invalid_argument("Invalid internal GRM storage.");
}

  
/// Diagonal Sigma = tau[0] * diag(1/W) + tau[1] * diag(GRM)
/// Input: w, tau
/// Output: out_sigma
static void get_diag_sigma(const dvec& w, const dvec& tau, dvec &out_sigma)
{
	out_sigma.resize(GRM.NumSamp);
	PARALLEL_RANGE(st, ed, GRM.NumSamp, false)
	{
		const double tau0 = tau[0], tau1 = tau[1];
		double *out = &out_sigma[0];
		for (size_t i=st; i < ed; i++)
		{
			double v = tau0 / w[i] + tau1 * GRM.buf_diag_grm[i];
			if (v < 1e-4) v = 1e-4;
			out[i] = v;
		}
	}
	PARALLEL_END
}


/// Diagonal Sigma = tau[0] * b * diag(1/W) + tau[1] * diag(GRM, b)
/// Input: b, w, tau
/// Output: out_sigma
static dvec get_crossprod(const dcolvec &b, const dvec& w, const dvec& tau)
{
	const double tau0 = tau[0], tau1 = tau[1];
	if (tau1 == 0)
	{
		return(tau0 * (b % (1/w)));
	} else {
		dvec out_b;
		get_crossprod_b_grm(b, out_b);
		return(tau0 * (b % (1/w)) + tau1 * out_b);
	}
}


/// PCG algorithm for diagonal of Sigma = tau[0] * diag(1/W) + tau[1] * GRM
/// Input: w, tau, b, maxiterPCG, tolPCG
static dvec PCG_diag_sigma(const dvec &w, const dvec &tau, const dvec &b,
		int maxiterPCG, double tolPCG)
{
	dvec r = b, r1, minv;
	get_diag_sigma(w, tau, minv);
	minv = 1 / minv;

	dvec z = minv % r, z1;
	dvec p = z;
	dvec x;
	x.zeros(GRM.NumSamp);

	int iter = 0;
	while ((iter < maxiterPCG) && (sum(r % r) > tolPCG))
	{
		iter = iter + 1;
		dvec Ap = get_crossprod(p, w, tau);
		double a = sum(r % z) / sum(p % Ap);
		x += a * p;
		r1 = r - a * Ap;
		z1 = minv % r1;

		double bet = sum(z1 % r1) / sum(z % r);
		p = z1 + bet*p;
		z = z1;
		r = r1;
	}

	if (iter >= maxiterPCG)
		Rprintf("PCG does not converge (may need to increase 'maxiter').\n");

	return(x);
}


/// Calculate the coefficient of variation for mean of a vector
inline static double calcCV(const dvec &x)
{
	double x_mean = mean(x);
	double x_sd = stddev(x);
	return(x_sd / (x_mean * int(x.n_elem)));
}

inline static double calcCV(const dvec &x, int st)
{
	return(calcCV(x(span(st, x.size()-1))));
}


/// Calculate the trace of matrix for binary outcomes
static double get_trace(const dmat &Sigma_iX, const dmat& X, const dvec& w,
	const dvec& tau, const dmat& cov, int nrun, int maxiterPCG, double tolPCG,
	double traceCVcutoff, int seed)
{
	set_seed(seed);
	dmat Sigma_iXt = Sigma_iX.t();
	dvec Sigma_iu;  
	dcolvec Pu;
	dvec Au, u;

	int nrunStart = 0;
	int nrunEnd = nrun;
	double traceCV = traceCVcutoff + 0.1;
	dvec buf;
	buf.zeros(nrun);

	// Hutchinson's randomized trace estimator
	while (traceCV > traceCVcutoff)
	{
		for (int i=nrunStart; i < nrunEnd; i++)
		{
			// buf[i] = Pu Au = (u' P) (G' G u) = u' P G' G u
			u = 2 * as<dvec>(rbinom(GRM.NumSamp, 1, 0.5)) - 1;
			Sigma_iu = PCG_diag_sigma(w, tau, u, maxiterPCG, tolPCG);
			Pu = Sigma_iu - Sigma_iX * (cov *  (Sigma_iXt * u));
			get_crossprod_b_grm(u, Au);
			buf[i] = dot(Au, Pu);
		}
		traceCV = calcCV(buf);
		if (traceCV > traceCVcutoff)
		{
			nrunStart = nrunEnd;
			nrunEnd = nrunEnd + 10;
			buf.resize(nrunEnd);
			Rprintf("CV for trace random estimator using %d runs is %g > %g\n",
				nrun, traceCV, traceCVcutoff);
			Rprintf("try %d runs ...\n", nrunEnd);
		}
	}

	return(mean(buf));
}


/// Calculate the trace of matrix for quantitative outcomes
static void get_trace_q(const dmat &Sigma_iX, const dmat& X, const dvec& w,
	const dvec& tau, const dmat& cov, int nrun, int maxiterPCG, double tolPCG,
	double traceCVcutoff, int seed, double &outTrace0, double &outTrace1)
{
	set_seed(seed);
	dmat Sigma_iXt = Sigma_iX.t();
	dvec Sigma_iu;  
	dcolvec Pu;
	dvec Au, u;

	int nrunStart = 0;
	int nrunEnd = nrun;
	double traceCV, traceCV0;
	traceCV = traceCV0 = traceCVcutoff + 0.1;
	dvec buf, buf0;
	buf.zeros(nrun); buf0.zeros(nrun);

	// Hutchinson's randomized trace estimator
	while ((traceCV > traceCVcutoff) || (traceCV0 > traceCVcutoff))
	{
		for (int i=nrunStart; i < nrunEnd; i++)
		{
			// buf[i] = Pu Au = (u' P) (G' G u) = u' P G' G u
			u = 2 * as<dvec>(rbinom(GRM.NumSamp, 1, 0.5)) - 1;
			Sigma_iu = PCG_diag_sigma(w, tau, u, maxiterPCG, tolPCG);
			Pu = Sigma_iu - Sigma_iX * (cov *  (Sigma_iXt * u));
			get_crossprod_b_grm(u, Au);
			buf[i]  = dot(Au, Pu);
			buf0[i] = dot(u, Pu);
		}
		traceCV  = calcCV(buf);
		traceCV0 = calcCV(buf0);
		if ((traceCV > traceCVcutoff) || (traceCV0 > traceCVcutoff))
		{
			nrunStart = nrunEnd;
			nrunEnd = nrunEnd + 10;
			buf.resize(nrunEnd);
			buf0.resize(nrunEnd);
			Rprintf("CV for trace random estimator using %d runs is %g > %g\n",
				nrun, traceCV, traceCVcutoff);
			Rprintf("try %d runs ...\n", nrunEnd);
		}
	}

	outTrace0 = mean(buf0);
	outTrace1 = mean(buf);
}


/// matrix inverse
inline static dmat mat_inv(const dmat &m)
{
	dmat rv, xs = symmatu(m);
	if (!auxlib::inv_sympd(rv, xs))
	{
		// xs is singular or not positive definite (possibly due to rounding error)
		// try inv(), if inv() still fails throw an exception and stop fitting
		Rprintf("Warning: arma::inv_sympd(), matrix is singular or not positive definite, use arma::inv() instead.\n");
		rv = inv(xs);
	}
	return rv;
}


/// Calculate fixed and random effect coefficients given Y, X, w, tau
/// Input:  Y, X, w, tau, maxiterPCG, tolPCG
/// Output: Sigma_iY, Sigma_iX, cov, alpha, eta
static void get_coeff_w(const dvec &Y, const dmat &X, const dvec &w,
	const dvec &tau, int maxiterPCG, double tolPCG,
	dvec &Sigma_iY, dmat &Sigma_iX, dmat &cov, dvec &alpha, dvec &eta)
{
	int n_col_X = X.n_cols;
	Sigma_iY = PCG_diag_sigma(w, tau, Y, maxiterPCG, tolPCG);
	// Sigma_iX = (X' Sigma^-1)'
	Sigma_iX.resize(GRM.NumSamp, n_col_X);
	dvec xv_i;
	for(int i = 0; i < n_col_X; i++)
	{
		xv_i = X.col(i);
		Sigma_iX.col(i) = PCG_diag_sigma(w, tau, xv_i, maxiterPCG, tolPCG);
	}
	// cov = (X' Sigma^-1 X)^-1
	cov = mat_inv(X.t() * Sigma_iX);
	// alpha = (X' Sigma^-1 X)^-1 X' Sigma^-1 Y
	alpha = cov * (Sigma_iX.t() * Y);
	eta = Y - tau[0] * (Sigma_iY - Sigma_iX * alpha) / w;
}


/// Calculate Sigma^-1 X
static dmat get_sigma_X(dvec &w, dvec &tau, dmat &X, int maxiterPCG,
	double tolPCG)
{
	int ncol = X.n_cols;
	dmat Sigma_iX1(GRM.NumSamp, ncol);
	for(int i = 0; i < ncol; i++)
	{
		Sigma_iX1.col(i) =
			PCG_diag_sigma(w, tau, X.col(i), maxiterPCG, tolPCG);
	}
	return(Sigma_iX1);
}


// ========================================================================= //

/// Calculate fixed and random effect coefficients given Y, X, tau
/// Input:  Y, X, tau, ...
/// Output: alpha, eta, W, ...
static void get_coeff(const dvec &y, const dmat &X, const dvec &tau,
	const List &family, const dvec &alpha0, const dvec &eta0,
	const dvec &offset, int maxiterPCG, int maxiter, double tolPCG,
	bool verbose,
	dvec &Y, dvec &mu, dvec &alpha, dvec &eta, dvec &W, dmat &cov,
	dvec &Sigma_iY, dmat &Sigma_iX)
{
	// initialize
	const double tol_coef = 0.1;
	Function fc_linkinv = wrap(family["linkinv"]);
	Function fc_mu_eta = wrap(family["mu.eta"]);
	Function fc_variance = wrap(family["variance"]);

	mu = as<dvec>(fc_linkinv(eta0));
	dvec mu_eta = as<dvec>(fc_mu_eta(eta0));
	Y = eta0 - offset + (y - mu)/mu_eta;
	W = (mu_eta % mu_eta) / as<dvec>(fc_variance(mu));

	// iterate ...
	dvec a0 = alpha0;
	for (int i=0; i < maxiter; i++)
	{
		// calculate fixed and random effect coefficients given Y, X, w, tau
		get_coeff_w(
			Y, X, W, tau, maxiterPCG, tolPCG,    // input variables
			Sigma_iY, Sigma_iX, cov, alpha, eta  // output
		);

		eta += offset;
		mu = as<dvec>(fc_linkinv(eta));
		mu_eta = as<dvec>(fc_mu_eta(eta));
		Y = eta - offset + (y - mu)/mu_eta;
		W = (mu_eta % mu_eta) / as<dvec>(fc_variance(mu));

		if (max(abs(alpha - a0)/(abs(alpha) + abs(a0) + tol_coef)) < tol_coef)
			break;
		a0 = alpha;
	}
	if (verbose)
	{
		print_vec("    tau: ", tau);
		print_vec("    fixed coeff: ", alpha);
	}
}


// Get average information (AI) for binary outcomes
static void get_AI_score(const dvec &Y, const dmat &X, const dvec &w,
	const dvec &tau, const dvec &Sigma_iY, const dmat &Sigma_iX, const dmat &cov,
	int nrun, int maxiterPCG, double tolPCG, double traceCVcutoff, int seed,
	double &outYPAPY, double &outTrace, double &outAI)
{
	dmat Sigma_iXt = Sigma_iX.t();
	dvec PY = Sigma_iY - Sigma_iX * (cov * (Sigma_iXt * Y));
	dvec APY;
	get_crossprod_b_grm(PY, APY);
	// output
	outYPAPY = dot(PY, APY);
	outTrace = get_trace(Sigma_iX, X, w, tau, cov, nrun, maxiterPCG, tolPCG,
		traceCVcutoff, seed);
	dvec PAPY_1 = PCG_diag_sigma(w, tau, APY, maxiterPCG, tolPCG);
	dvec PAPY = PAPY_1 - Sigma_iX * (cov * (Sigma_iXt * PAPY_1));
	outAI = dot(APY, PAPY);
}

// Get average information (AI) for quantitative outcomes
static void get_AI_score_q(const dvec &Y, const dmat &X, const dvec &w,
	const dvec &tau, const dvec &Sigma_iY, const dmat &Sigma_iX, const dmat &cov,
	int nrun, int maxiterPCG, double tolPCG, double traceCVcutoff, int seed,
	double outYPAPY[], double outTrace[], dmat &outAI)
{
	dmat Sigma_iXt = Sigma_iX.t();
	dvec PY = Sigma_iY - Sigma_iX * (cov * (Sigma_iXt * Y));
	dvec A0PY = PY;
	dvec APY;
	get_crossprod_b_grm(PY, APY);
	// get YPAPY
	outYPAPY[0] = dot(PY, APY);   // YPAPY
	outYPAPY[1] = dot(PY, A0PY);  // YPA0PY
	// get Trace
	get_trace_q(Sigma_iX, X, w, tau, cov, nrun, maxiterPCG, tolPCG,
		traceCVcutoff, seed, outTrace[0], outTrace[1]);
	// get AI
	dmat AI(2,2);
	dvec PA0PY_1 = PCG_diag_sigma(w, tau, A0PY, maxiterPCG, tolPCG);
	dvec PA0PY = PA0PY_1 - Sigma_iX * (cov * (Sigma_iXt * PA0PY_1));
	AI(0,0) = dot(A0PY, PA0PY);
	dvec PAPY_1 = PCG_diag_sigma(w, tau, APY, maxiterPCG, tolPCG);
	dvec PAPY = PAPY_1 - Sigma_iX * (cov * (Sigma_iXt * PAPY_1));
	AI(1,1) = dot(APY, PAPY);
	AI(1,0) = AI(0,1) = dot(A0PY, PAPY);
	outAI = AI;
}


// Update tau for binary outcomes
static dvec fitglmmaiRPCG(const dvec &Y, const dmat &X, const dvec &w,
	const dvec &in_tau, const dvec &Sigma_iY, const dmat &Sigma_iX,
	const dmat &cov, int nrun, int maxiterPCG, double tolPCG, double tol,
	double traceCVcutoff, int seed)
{
	double YPAPY, Trace, AI;
	get_AI_score(Y, X, w, in_tau, Sigma_iY, Sigma_iX, cov, nrun,
		maxiterPCG, tolPCG, traceCVcutoff, seed,
		YPAPY, Trace, AI);
  	double score = YPAPY - Trace;
  	double Dtau = score / AI;
  	dvec tau = in_tau;
  	dvec tau0 = in_tau;
  	tau[1] = tau0[1] + Dtau;

  	for(size_t i=0; i < tau.n_elem; i++)
		if (tau[i] < tol) tau[i] = 0;

  	double step = 1.0;
  	while (tau[1] < 0.0)
  	{
		step *= 0.5;
		tau[1] = tau0[1] + step * Dtau;
  	}

  	for(size_t i=0; i < tau.n_elem; i++)
		if (tau[i] < tol) tau[i] = 0;

  	return tau;
}

// Update tau for quantitative outcomes
static dvec fitglmmaiRPCG_q(const dvec &Y, const dmat &X, const dvec &w,
	const dvec &in_tau, const dvec &Sigma_iY, const dmat &Sigma_iX,
	const dmat &cov, int nrun, int maxiterPCG, double tolPCG, double tol,
	double traceCVcutoff, int seed)
{
	uvec zero_v = (in_tau < tol);
	double YPAPY[2], Trace[2];
	dmat AI;
	get_AI_score_q(Y, X, w, in_tau, Sigma_iY, Sigma_iX, cov, nrun,
		maxiterPCG, tolPCG, traceCVcutoff, seed,
		YPAPY, Trace, AI);
	dvec score(2);
	score[0] = YPAPY[1] - Trace[0];
	score[1] = YPAPY[0] - Trace[1];
	dvec Dtau = solve(AI, score);

	dvec tau0 = in_tau;
	dvec tau = tau0 + Dtau;
	tau.elem( find(zero_v % (tau < tol)) ).zeros();

	double step = 1.0;
	while (tau[0] < 0.0 || tau[1]  < 0.0)
	{
		step *= 0.5;
		tau = tau0 + step * Dtau;
		tau.elem( find(zero_v % (tau < tol)) ).zeros();
	}
	tau.elem( find(tau < tol) ).zeros();

	return tau;
}


// ========================================================================= //

// Fitting the null model
RcppExport SEXP saige_fit_AI_PCG(SEXP r_fit0, SEXP r_X, SEXP r_tau,
	SEXP r_param)
{
BEGIN_RCPP

	// parameters for fitting the model
	List param(r_param);
	const int trait = param["trait"];
	const double tol = param["tol"];
	const double tol_inv_2 = 1 / (tol*tol);
	const double tolPCG = param["tolPCG"];
	const int seed = param["seed"];
	const int maxiter = param["maxiter"];
	const int maxiterPCG = param["maxiterPCG"];
	const bool no_iteration = Rf_asLogical(param["no_iteration"])==TRUE;
	const int nrun = param["nrun"];
	const double traceCVcutoff = param["traceCVcutoff"];
	const bool verbose = Rf_asLogical(param["verbose"])==TRUE;
	const char *indent = param["indent"];

	List fit0(r_fit0);
	dvec y = as<dvec>(fit0["y"]);
	const int n = y.size();
	dmat X = as<dmat>(r_X);
	dvec offset(y.size());
	if (Rf_isNull(fit0["offset"]))
		offset.zeros();
	else
		offset = as<dvec>(fit0["offset"]);

	List family = fit0["family"];
	Function fc_mu_eta = wrap(family["mu.eta"]);
	dvec eta = as<dvec>(fit0["linear.predictors"]);
	dvec eta0 = eta;
	dvec mu = as<dvec>(fit0["fitted.values"]);
	dvec mu_eta = as<dvec>(fc_mu_eta(eta0));
	dvec Y = eta - offset + (y - mu) / mu_eta;
	dvec alpha0 = as<dvec>(fit0["coefficients"]);
	dvec alpha = alpha0;
	dmat cov;
	dvec tau = as<dvec>(r_tau);
	dvec tau0 = tau;
	const double tau_init_sum = sum(tau);
	int ntry_Sigma_E_zero=0;
	int iter = 1;

	if (verbose && !no_iteration)
	{
		Rprintf(
			"%sInitial variance component estimates, tau (Sigma_E, Sigma_G):\n",
			indent);
	}

	PARALLEL_THREAD_BLOCK

	dvec re_Y, re_mu, re_alpha, re_eta, re_W, re_Sigma_iY;
	dmat re_cov, re_Sigma_iX;
	get_coeff(
		// input variables
		y, X, tau, family, alpha0, eta0, offset, maxiterPCG, maxiter,
		tolPCG, verbose,
		// output
		re_Y, re_mu, re_alpha, re_eta, re_W, re_cov, re_Sigma_iY, re_Sigma_iX
	);

/*
	if (no_iteration)
	{
		return List::create(
			_["coefficients"] = SEXP_VEC(re_alpha),
			_["tau"] = SEXP_VEC(tau),
			_["linear.predictors"] = SEXP_VEC(re_eta),
			_["fitted.values"] = SEXP_VEC(re_mu),
			_["residuals"] = SEXP_VEC(y - re_mu),
			_["cov"] = re_cov,
			_["converged"] = true);
	}
*/

	double YPAPY[2], Trace[2];
	switch (trait)
	{
	case TTrait::Quant:
		{
			dmat AI;
			get_AI_score_q(
				// input variables
				re_Y, X, re_W, tau, re_Sigma_iY, re_Sigma_iX, re_cov, nrun,
				maxiterPCG, tolPCG, traceCVcutoff, seed,
				// output
				YPAPY, Trace, AI
			);
			tau[0] = std::max(0.0, tau0[0] + sq(tau0[0])*(YPAPY[1]-Trace[0])/n);
			tau[1] = std::max(0.0, tau0[1] + sq(tau0[1])*(YPAPY[0]-Trace[1])/n);
			break;
		}
	case TTrait::Binary:
		{
			double AI;
			get_AI_score(
				// input variables
				re_Y, X, re_W, tau, re_Sigma_iY, re_Sigma_iX, re_cov, nrun, maxiterPCG,
				tolPCG, traceCVcutoff, seed,
				// output
				YPAPY[0], Trace[0], AI
			);
			tau[1] = std::max(0.0, tau0[1] + sq(tau0[1])*(YPAPY[0]-Trace[0])/n);
			break;
		}
	default:
		throw std::invalid_argument("Invalid trait.");
	}

	// iterate
	for (; iter <= maxiter; iter++)
	{
		if (verbose)
			Rprintf("%sIteration %d:\n", indent, iter);

		// save the old values
		alpha0 = re_alpha; tau0 = tau; eta0 = eta;

		// find the next tau, try multiple times if fails
		for (int itry=1; itry <= 11; itry++)
		{
			get_coeff(
				// input variables
				y, X, tau0, family, alpha0, eta0, offset, maxiterPCG, maxiter,
				tolPCG, verbose,
				// output
				re_Y, re_mu, re_alpha, re_eta, re_W, re_cov, re_Sigma_iY, re_Sigma_iX
			);
			// update tau
			switch (trait)
			{
			case TTrait::Quant:
				tau = fitglmmaiRPCG_q(re_Y, X, re_W, tau0, re_Sigma_iY,
					re_Sigma_iX, re_cov, nrun, maxiterPCG, tolPCG, tol,
					traceCVcutoff, seed);
				break;
			case TTrait::Binary:
				tau = fitglmmaiRPCG(re_Y, X, re_W, tau0, re_Sigma_iY,
					re_Sigma_iX, re_cov, nrun, maxiterPCG, tolPCG, tol,
					traceCVcutoff, seed);
				break;
			}
			// check new tau
			if (max(tau) > tol_inv_2)
			{
				if (itry <= 10)  // try at most 10 times
				{
					tau0[1] *= 0.5;
					if (verbose)
					{
						print_vec("    [Warning] tau: ", tau, indent, false);
						Rprintf(", large variance estimate observed, retry (%d) ...\n", itry);
						print_vec("    set new tau: ", tau0, indent);
					}
					continue;
				} else {
					if (verbose)
						print_vec("tau: ", tau, indent);
					throw std::overflow_error(
					"Large variance estimate observed in the iterations, model not converged!");
				}
			}
			break;  // <-- quit for-loop
		}

		// update coefficients
		cov = re_cov; alpha = re_alpha; eta = re_eta;
		Y = re_Y; mu = re_mu;

		if (tau[0] <= 0)
		{
			print_vec("    tau: ", tau, indent);
			const int ntry_max = 5;  // try at most 5 times
			if (ntry_Sigma_E_zero < ntry_max)
			{
				ntry_Sigma_E_zero ++;
				if (ntry_Sigma_E_zero < ntry_max)
				{
					tau[1] = tau_init_sum / (ntry_Sigma_E_zero+2);
					tau[0] = tau_init_sum - tau[1];
				} else {
					tau[0] = tau_init_sum;  // last try: no random effect
					tau[1] = 0;
					tau0 = tau;
				}
				Rprintf("    [Warning] Sigma_E = 0, retry (%d) using new Sigma_E & Sigma_G ...\n",
					ntry_Sigma_E_zero);
				print_vec("    tau: ", tau, indent);
				if (ntry_Sigma_E_zero < ntry_max)
					continue;
			} else {
				throw std::overflow_error("Sigma_E = 0, model not converged!");
			}
		}
		if (trait==TTrait::Binary && tau[1]==0)
			break;
		if (max(abs(tau-tau0) / (abs(tau)+abs(tau0)+tol)) < tol)
			break;
	}

	get_coeff(
		// input variables
		y, X, tau, family, alpha0, eta0, offset, maxiterPCG, maxiter,
		tolPCG, false,
		// output
		re_Y, re_mu, re_alpha, re_eta, re_W, re_cov, re_Sigma_iY, re_Sigma_iX
	);
	cov = re_cov; alpha = re_alpha; eta = re_eta;
	Y = re_Y; mu = re_mu;

	PARALLEL_THREAD_BLOCK_END

	if (verbose)
	{
		print_vec("Final tau: ", tau, indent);
		print_vec("    fixed coeff: ", alpha, indent);
	}

	return List::create(
		_["coefficients"] = SEXP_VEC(alpha),
		_["tau"] = SEXP_VEC(tau),
		_["linear.predictors"] = SEXP_VEC(eta),
		_["fitted.values"] = SEXP_VEC(mu),
		_["residuals"] = SEXP_VEC(y - mu),
		_["cov"] = cov,
		_["converged"] = bool(iter <= maxiter));

END_RCPP
}


// ========================================================================= //

// Get the diagnoal of the full GRM
RcppExport SEXP saige_get_grm_diag()
{
	if (GRM.PackedG || GRM.SparseG)
	{
		SEXP rv = NEW_NUMERIC(GRM.NumSamp);
		memcpy(REAL(rv), GRM.buf_diag_grm, sizeof(double)*GRM.NumSamp);
		return rv;
	} else {
		Rf_error("No full GRM.");
		return R_NilValue;
	}
}


/// Get a dosage vector for a specific SNP from a packed RAW matrix
static void get_geno_2b_ds(RawMatrix &G, int col_i, dvec &ds)
{
	// check
	if ((col_i < 0) || (col_i >= G.ncol()))
		throw std::invalid_argument("Invalid column index.");
	// dense genotypes
	ds.resize(GRM.NumSamp);
	const BYTE *g = &G[G.nrow() * col_i];
	double *p = &ds[0];
	size_t n = GRM.NumSamp;
	for (; n >= 4; n-=4, p+=4)
	{
		BYTE gg = *g++;
		p[0] = ds_nan(gg & 0x03);
		p[1] = ds_nan((gg >> 2) & 0x03);
		p[2] = ds_nan((gg >> 4) & 0x03);
		p[3] = ds_nan(gg >> 6);
	}
	for (BYTE gg = (n>0 ? *g : 0); n > 0; n--)
	{
		*p++ = ds_nan(gg & 0x03);
		gg >>= 2;
	}
}

/// Set genotypes according to indexing
RcppExport SEXP saige_set_geno2b_raw(SEXP gm_raw, SEXP row_ii, SEXP col_i)
{
static const BYTE mask[4] = { 0x03, 0x0C, 0x30, 0xC0 };
BEGIN_RCPP
	RawMatrix G(gm_raw);
	BYTE *p = &G(0, Rf_asInteger(col_i)-1);
	const int n = Rf_length(row_ii);
	int *ii = INTEGER(row_ii);
	for (int j=0; j < n; j++)
	{
		const int k = (ii[j] - 1) / 2;  // sample index
		const int i0 = k / 4, i1 = k % 4;
		const BYTE mask1 = mask[i1], mask2 = ~mask1;
		const BYTE new_v = ((p[i0] & mask1) >> (i1*2)) + 1;
		p[i0] = (p[i0] & mask2) | (new_v << (i1*2));
	}
END_RCPP
}


/// Calculate variance ratio for binary outcomes
RcppExport SEXP saige_calc_var_ratio(SEXP r_fit0, SEXP r_glmm, SEXP r_noK,
	SEXP r_param)
{
BEGIN_RCPP

	// output variables
	vector<int> lst_idx;
	vector<double> lst_maf, lst_mac, lst_var1, lst_var2, lst_ratio;
	const char *nm_var2;
	// multithreading blocks
	PARALLEL_THREAD_BLOCK

	List fit0(r_fit0);
	List glmm(r_glmm);
	List obj_noK(r_noK);
	List param(r_param);

	// parameters for fitting the model
	const int trait = param["trait"];
	const double tolPCG = param["tolPCG"];
	const int maxiterPCG = param["maxiterPCG"];
	const double ratioCVcutoff = param["ratioCVcutoff"];
	const int num_marker = param["num.marker"];
	const bool verbose = Rf_asLogical(param["verbose"])==TRUE;

	List family = fit0["family"];
	Function fc_mu_eta = wrap(family["mu.eta"]);
	Function fc_variance = wrap(family["variance"]);

	dvec eta = as<dvec>(glmm["linear.predictors"]);
	dvec mu = as<dvec>(glmm["fitted.values"]);
	dvec mu_eta = as<dvec>(fc_mu_eta(eta));
	dvec W = (mu_eta % mu_eta) / as<dvec>(fc_variance(mu));
	dvec tau = as<dvec>(glmm["tau"]);
	dmat X1 = as<dmat>(obj_noK["X1"]);
	dmat Sigma_iX = get_sigma_X(W, tau, X1, maxiterPCG, tolPCG);

	dvec y = as<dvec>(fit0["y"]);
	dmat noK_XXVX_inv = as<dmat>(obj_noK["XXVX_inv"]);
	dmat noK_XV = as<dmat>(obj_noK["XV"]);

	// user-defined GRM
	SEXP sigma_inv = obj_noK["Sigma_inv"];
	const bool user_def_GRM = Rf_isLogical(sigma_inv) == 0;
	nm_var2 = user_def_GRM ? "var2_u" : "var2";
	dmat Sigma_inv_dm;
	sp_mat Sigma_inv_sp;
	bool is_Sigma_inv_sp = false;
	if (user_def_GRM)
	{
		is_Sigma_inv_sp = Rf_isMatrix(sigma_inv) == 0;
		if (is_Sigma_inv_sp)
			Sigma_inv_sp = as<sp_mat>(sigma_inv);
		else
			Sigma_inv_dm = as<dmat>(sigma_inv);
	}

	dvec G0(GRM.NumSamp);
	vector<int> buf_idx(GRM.NumSamp);

	// for-loop each genotype matrix
	const int n_VR = Rf_length(Geno_PackedRawRandVR);
	int i_start = 0;
	for (int i_VR=0; i_VR < n_VR; i_VR++)
	{
		if (verbose && n_VR > 1)
			Rprintf("    #%d MAC category:\n", i_VR+1);
		RawMatrix G_VR(VECTOR_ELT(Geno_PackedRawRandVR, i_VR));
		const int n_tot_snp = G_VR.ncol();
		double ratioCV = ratioCVcutoff + 0.1;
		int i_num_marker = num_marker;
		int lst_st_idx = lst_ratio.size();

		for (int i=0; ratioCV > ratioCVcutoff && i < n_tot_snp; )
		{
			bool indent = true;
			while (i < i_num_marker && i < n_tot_snp)
			{
				// get the genotypes
				get_geno_2b_ds(G_VR, i, G0);
				double AF, AC; int Num;
				f64_af_ac_impute(&G0[0], GRM.NumSamp, AF, AC, Num, &buf_idx[0], 2);
				if (AF > 0.5)
				{
					f64_sub(GRM.NumSamp, 2, &G0[0]);
					AC = 2*Num - AC;
					AF = 1 - AF;
				}

				// adjusted genotypes
				dvec G = G0 - noK_XXVX_inv * (noK_XV * G0);
				dvec Sigma_iG = PCG_diag_sigma(W, tau, G, maxiterPCG, tolPCG);
				dvec adjG = Sigma_iX * mat_inv(X1.t() * Sigma_iX) * X1.t() * Sigma_iG;

				// variance ratio
				double var1 = sum(G % Sigma_iG) - sum(G % adjG);
				double var2;
				if (user_def_GRM)
				{
					if (is_Sigma_inv_sp)
						var2 = sum(G % (Sigma_inv_sp * G));
					else
						var2 = sum(G % (Sigma_inv_dm * G));
				} else {
					switch (trait)
					{
						case TTrait::Quant:
							var2 = sum(G % G) / tau[0]; break;
						case TTrait::Binary:
							var2 = sum(mu % (1 - mu) % G % G); break;
						default:
							throw std::invalid_argument("Invalid trait.");
					}
				}
				double ratio = var1 / var2;

				// not necessary, but correspond to the original SAIGE package
				const double InvAC = 1/AC;
				var1 *= InvAC; var2 *= InvAC;

				i++;
				lst_idx.push_back(i_start + i);
				lst_maf.push_back(AF);    lst_mac.push_back(AC);
				lst_var1.push_back(var1); lst_var2.push_back(var2);
				lst_ratio.push_back(ratio);
				if (verbose)
				{
					if (i <= 5)
					{
						Rprintf("    %d, maf: %0.5f, mac: %g,\tratio: %0.4f (var1: %.3g, %s: %.3g)\n",
							i, AF, AC, ratio, var1, nm_var2, var2);
					} else {
						Rprintf(indent ? "    ." : ".");
						indent = false;
					}
				}
			}
			if (verbose && (i > 5)) Rprintf("\n");

			// check the coefficient of variation (CV)
			ratioCV = calcCV(lst_ratio, lst_st_idx);
			if (ratioCV > ratioCVcutoff)
			{
				if (verbose)
				{
					Rprintf(
						"CV for variance ratio estimate using %d markers is %g > ratioCVcutoff (%g), try more markers ...\n",
						i_num_marker, ratioCV, ratioCVcutoff);
				}
				i_num_marker += 10;
			}
		}

		i_start += n_tot_snp;
		if (ratioCV > ratioCVcutoff)
		{
			Rprintf(
				"No more markers (n=%d) for variance ratio estimate (still CV %g > ratioCVcutoff %g)\n",
				n_tot_snp, ratioCV, ratioCVcutoff);
		}
	}
	PARALLEL_THREAD_BLOCK_END

	return DataFrame::create(
		_["id"]    = lst_idx,
		_["maf"]   = lst_maf,  _["mac"]   = lst_mac,
		_["var1"]  = lst_var1, _[nm_var2] = lst_var2, _["ratio"] = lst_ratio
	);

END_RCPP
}
