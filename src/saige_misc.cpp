// ===========================================================
//
// saige_misc.cpp: Miscellaneous functions
//
// Copyright (C) 2022-2023    Xiuwen Zheng / AbbVie-ComputationalGenomics
//
// This file is part of SAIGEgds.
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
#include <vector>
#include <algorithm>

using namespace std;
using namespace Rcpp;
using namespace arma;


namespace SAIGE
{

/// the number of threads used in the SAIGEgds package
int SAIGE_NumThread = 0;


void Type_Matrix::reset(SEXP mat)
{
	if (mat)
	{
		if (!Rf_isMatrix(mat))
			throw std::invalid_argument("Invalid argument in Type_Matrix::reset().");
		val = REAL(mat);
		NumericMatrix M(mat);
		m_nrow = M.nrow(); m_ncol = M.ncol();
	} else {
		val = NULL;
		m_nrow = m_ncol = 0;
	}
}

static void get_i_p_x(SEXP mat, const int *&i, const int *&p, const double *&x,
	int &nrow, int &ncol)
{
	S4 M(mat);
	i = INTEGER(M.slot("i"));
	p = INTEGER(M.slot("p"));
	x = REAL(M.slot("x"));
	IntegerVector dm = M.slot("Dim");
	nrow = dm[0]; ncol = dm[1];
}

void Type_dgCMatrix::reset(SEXP mat)
{
	if (mat)
	{
		if (!Rf_inherits(mat, "dgCMatrix"))
			throw std::invalid_argument("Invalid argument in Type_dgCMatrix::reset().");
		get_i_p_x(mat, i, p, x, m_nrow, m_ncol);
	} else {
		i = p = NULL; x = NULL;
		m_nrow = m_ncol = 0;
	}
}

void Type_dsCMatrix::reset(SEXP mat)
{
	if (mat)
	{
		if (!Rf_inherits(mat, "dsCMatrix"))
			throw std::invalid_argument("Invalid argument in Type_dsCMatrix::reset().");
		get_i_p_x(mat, i, p, x, m_nrow, m_ncol);
	} else {
		i = p = NULL; x = NULL;
		m_nrow = m_ncol = 0;
	}
}

}


namespace Misc
{

/// error information
static const char *ERR_INVALID_MAT = "Invalid input dense matrix!";
/// error information
static const char *ERR_INVALID_SP_MAT = "Invalid input sparse matrix!";


/// Calculate AF & MAC from a sparse genotype matrix, return the number of columns
int SummaryStat_Mat(SEXP mat, double out_af[], double out_mac[])
{
	// check
	if (!Rf_isMatrix(mat))
		Rf_error("The input matrix should be a dense matrix.");
	IntegerVector Dim = Rf_getAttrib(mat, R_DimSymbol);
	if (Dim.size() != 2) Rf_error(ERR_INVALID_MAT);
	// calculate AF & MAC for each SNV
	const int nrow = Dim[0];
	const int ncol = Dim[1];
	int i = 0;
	switch (TYPEOF(mat))
	{
	case RAWSXP:
		for (const Rbyte *p = RAW(mat); i < ncol; i++, p+=nrow)
		{
			int n=nrow, s=0;
			for (int j=0; j < nrow; j++)
				if (p[j] != Rbyte(0xFF)) s+=p[j]; else n--;
			out_af[i] = (n > 0) ? (double)s/(2*n) : R_NaN;
			out_mac[i] = std::min(s, 2*n-s);
		}
		break;
	case INTSXP:
		for (const int *p = INTEGER(mat); i < ncol; i++, p+=nrow)
		{
			int n=nrow, s=0;
			for (int j=0; j < nrow; j++)
				if (p[j] != NA_INTEGER) s+=p[j]; else n--;
			out_af[i] = (n > 0) ? (double)s/(2*n) : R_NaN;
			out_mac[i] = std::min(s, 2*n-s);
		}
		break;
	case REALSXP:
		for (const double *p = REAL(mat); i < ncol; i++, p+=nrow)
		{
			int n=nrow;
			double s=0;
			for (int j=0; j < nrow; j++)
				if (R_FINITE(p[j])) s+=p[j]; else n--;
			out_af[i] = (n > 0) ? s/(2*n) : R_NaN;
			out_mac[i] = std::min(s, 2*n-s);
		}
		break;
	default:
		Rf_error(ERR_INVALID_MAT);
	}
	// output
	return ncol;
}


/// Calculate AF & MAC from a sparse genotype matrix, return the number of columns
int SummaryStat_SpMat(SEXP mat, double out_af[], double out_mac[])
{
	// check
	if (!Rf_inherits(mat, "dgCMatrix"))
		Rf_error("The input matrix should be dgCMatrix.");
	S4 ObjM(mat);
	IntegerVector I = ObjM.slot("i");
	IntegerVector P = ObjM.slot("p");
	NumericVector X = ObjM.slot("x");
	IntegerVector Dim = ObjM.slot("Dim");
	if (Dim.size() != 2)
		Rf_error(ERR_INVALID_SP_MAT);
	if ((P.size() != Dim[1]+1) || (I.size() != X.size()))
		Rf_error(ERR_INVALID_SP_MAT);
	// calculate AF & MAC for each SNV
	const int nrow = Dim[0];
	const int ncol = Dim[1];
	const double *pX = &X[0];
	for (int i=0; i < ncol; i++)
	{
		const int st = P[i], ed = P[i+1];
		double sum = 0;
		int n = nrow;
		for (int j=st; j < ed; j++)
		{
			if (R_FINITE(pX[j]))
				sum += pX[j];
			else
				n--;
		}
		if (n > 0)
		{
			out_af[i] = sum / (2*n);
			out_mac[i] = std::min(sum, 2*n-sum);
		} else
			out_af[i] = out_mac[i] = R_NaN;
	}
	// output
	return ncol;
}


/// Return a sp_mat with imputed and flipped genotypes if needed
///    excluding monomorphic variants
sp_mat GetSp_Impute_SpMat(SEXP mat, double af[], double mac[],
	double mac_imp[])
{
	// check
	if (!Rf_inherits(mat, "dgCMatrix"))
		Rf_error("The input matrix should be dgCMatrix.");
	S4 ObjM(mat);
	IntegerVector I = ObjM.slot("i");
	IntegerVector P = ObjM.slot("p");
	NumericVector X = ObjM.slot("x");
	IntegerVector Dim = ObjM.slot("Dim");
	// calculate # of non-zero values
	const int nrow = Dim[0];
	const int ncol = Dim[1];
	const double *pX = &X[0];
	const int *pI = &I[0];
	size_t nnzero=0, new_ncol=0;
	for (int i=0; i < ncol; i++)
	{
		const double AF = af[i];
		if (R_FINITE(AF) && 0<AF && AF<1)
		{
			new_ncol ++;
			const int st = P[i], ed = P[i+1];
			if (AF <= 0.5)
			{
				nnzero += (ed - st);
			} else {
				// need flipping
				nnzero += nrow;
				for (int j=st; j < ed; j++)
				{
					if (R_FINITE(pX[j]) && pX[j]==2)
						nnzero--;
				}
			}
		}
	}
	// initialize the components of sp_mat
	uvec rowidx(nnzero);
	uvec colptr(new_ncol+1);
	dvec val(nnzero);
	size_t idx=0, new_idx_col=0;
	colptr[0] = idx;
	for (int i=0; i < ncol; i++)
	{
		const double AF = af[i];
		if (R_FINITE(AF) && 0<AF && AF<1)
		{
			const size_t old_idx = idx;
			const int st = P[i], ed = P[i+1];
			if (AF <= 0.5)
			{
				const double mean = AF*2;  // mean imputation
				for (int j=st; j < ed; j++)
				{
					rowidx[idx] = pI[j];
					val[idx] = R_FINITE(pX[j]) ? pX[j] : mean;
					idx ++;
				}
			} else {
				// need flipping
				af[i] = 1 - af[i];
				// check strictly increasing order
				for (int j=st+1; j < ed; j++)
				{
					if (pI[j-1] >= pI[j])
						Rf_error("Row indices in dgCMatrix should be strictly increasing.");
				}
				// fill
				const double mean = 2 - AF*2;  // mean imputation
				int i_j = 0, j = st;
				for (; i_j<nrow && j<ed; i_j++)
				{
					const int ii = pI[j];
					if (i_j < ii)
					{
						rowidx[idx] = i_j; val[idx] = 2; idx++;
					} else if (i_j == ii)
					{
						if (R_FINITE(pX[j]))
						{
							if (pX[j] != 2)
								{ rowidx[idx] = i_j; val[idx] = 2-pX[j]; idx++; }
						} else {
							rowidx[idx] = i_j; val[idx] = mean; idx++;
						}
						j++;  // next non-zero
					} else
						Rf_error("Internal error in GetSp_Impute_SpMat() using dgCMatrix.");
				}
				// fill end
				for (; i_j < nrow; i_j++)
					{ rowidx[idx] = i_j; val[idx] = 2; idx++; }
			}
			colptr[new_idx_col+1] = idx;
			af[new_idx_col] = af[i];
			mac[new_idx_col] = mac[i];
			// update MAC for missing genotypes
			double s = 0;
			for (size_t j=old_idx; j < idx; j++) s += val[j];
			mac_imp[new_idx_col] = s;
			// next column
			new_idx_col ++;
		}
	}
	if ((idx != nnzero) || (new_idx_col != new_ncol))
		Rf_error("Internal error in GetSp_Impute_SpMat().");
	// output
	return sp_mat(rowidx, colptr, val, nrow, new_ncol);
}


/// Collapse ultra rare variants with imputed genotype matrix
///     collapse_method = 1, presence or absence (PA)
///     collapse_method = 2, presence or absence (PA_int)
///     collapse_method = 3, sum up rare genotype (SumG)
sp_mat GetSp_CollapseGenoMat(const sp_mat &mat, double collapse_mac,
	int collapse_method,
	const double mac[], double inout_maf[], int &out_n_collapse)
{
	const int nrow = mat.n_rows;
	const int ncol = mat.n_cols;
	// check whether has any collapsing
	int n_collapse=0;
	for (int i=0; i < ncol; i++)
		if (mac[i] <= collapse_mac) n_collapse ++;
	out_n_collapse = n_collapse;
	if (n_collapse == 0) return mat;
	// need collapsing for n_collapse > 0
	uvec icol(ncol - n_collapse);
	int icol_st=0, n_g_c1=0;
	dvec g_c1;
	g_c1.zeros(nrow);
	// for-loop
	for (int i=0; i < ncol; i++)
	{
		if (mac[i] <= collapse_mac)
		{
			if (mac[i] > 0)
			{
				n_g_c1 ++;
				sp_mat::const_iterator it = mat.begin_col(i);
				sp_mat::const_iterator ed = mat.end_col(i);
				switch (collapse_method)
				{
				case 1:  // presence or absence (PA)
				case 2:  // presence or absence (PA_int), get rowMax
					for (; it != ed; ++it)
					{
						double g = *it;
						if (g > g_c1[it.row()]) g_c1[it.row()] = g;
					}
					break;
				case 3:  // sum up genotype (SumG)
					for (; it != ed; ++it) g_c1[it.row()] += *it;
					break;
				default:
					Rf_error("Invalid 'collapse_method: %d'.", collapse_method);
				}
			}
		} else
			icol[icol_st++] = i;
	}
	// set maf output
	for (size_t i=0; i < icol.size(); i++)
		inout_maf[i] = inout_maf[icol[i]];
	// output
	if (n_g_c1 > 0)
	{
		switch (collapse_method)
		{
		case 1:  // presence or absence (PA)
			for (size_t i=0; i < g_c1.n_elem; i++)
			{
				double &v = g_c1[i];
				if (v >= 1.5) v = 2; else if (v >= 0.5) v = 1;
			}
			break;
		case 2:  // presence or absence (PA_int: 0, 1, 2)
			for (size_t i=0; i < g_c1.n_elem; i++)
			{
				double &v = g_c1[i];
				if (v >= 1.5) v = 2; else if (v >= 0.5) v = 1; else v = 0;
			}
			break;
		}
		inout_maf[icol.size()] = sum(g_c1) / (2*g_c1.size());
		sp_mat sp = mat.cols(icol);
		return join_rows(sp, sp_mat(g_c1));
	} else
		return mat.cols(icol);
}

}


// ========================================================================= //
// Sparse GRM calculation in seqFitSparseGRM()

using namespace SAIGE;

RcppExport SEXP saige_grm_sp_reraw(SEXP g_pack, SEXP idx, SEXP r_buf)
{
BEGIN_RCPP
	RawMatrix G_Pack(g_pack);
	const size_t nrow = G_Pack.nrow();
	const size_t ncol = G_Pack.ncol();
	const size_t nsnp = nrow * 4;
	const int *p_idx = INTEGER(idx);
	BYTE *p_buf = RAW(r_buf);
	// for each sample
	for (size_t i=0; i < ncol; i++)
	{
		BYTE *p = &G_Pack(0, i);
		memset(p_buf, 0, nrow);
		for (size_t j=0; j < nsnp; j++)
		{
			size_t k = p_idx[j];
			BYTE v = (p[k >> 2] >> ((k & 0x03) * 2)) & 0x03;
			p_buf[j >> 2] |= v << ((j & 0x03) * 2);
		}
		memcpy(p, p_buf, nrow);
	}
	// return
	return R_NilValue;
END_RCPP
}


/// Indexing the non-zero value in a sparse matrix
struct t_sp_i_j
{
	int i, j;
	t_sp_i_j() {}
	t_sp_i_j(int vi, int vj): i(vi), j(vj) {}
};

#define GRM_USE_FLOAT_FOR_REDUCED_SET

#ifdef GRM_USE_FLOAT_FOR_REDUCED_SET
#   define USE_FLOAT_BOOL false
#   define USE_FLOAT_TYPE float
#else
#   define USE_FLOAT_BOOL true
#   define USE_FLOAT_TYPE double
#endif


/// Initialize the normalized genotype look-up table
inline static void grm_sp_init_lookup(SEXP g_pack, SEXP g_lookup, bool use_f64)
{
	RawMatrix G_Pack(g_pack);
	const size_t nrow = G_Pack.nrow();
	const size_t ncol = G_Pack.ncol();
	NumericMatrix G_Lookup(g_lookup);
	void *p_G = &G_Lookup[0];

	PARALLEL_THREAD_BLOCK
	PARALLEL_FOR(i, nrow, true)
	{
		int n[4] = { 0, 0, 0, 0 }, sum[4] = { 0, 0, 0, 0 };
		BYTE *p = &G_Pack[i];
		// get n & sum
		for (size_t j=0; j < ncol; j++)
		{
			BYTE g = *p; p += nrow;
			for (size_t k=0; k < 4; k++)
			{
				BYTE b = g & 0x03; g >>= 2;
				n[k] += (b < 3);
				sum[k] += (b < 3) ? b : 0;
			}
		}
		// output
		if (use_f64)
		{
			double *F = ((double *)p_G) + i*4*8;
			for (size_t k=0; k < 4; k++, F+=8)
			{
				double f = (n[k] > 0) ? double(sum[k])/(2*n[k]) : R_NaN;
				double d = 1 / sqrt(2*f*(1-f));
				if (!R_FINITE(d)) d = f = 0;
				const double f2 = 2*f;
				// normalized genotypes
				const double g0=(0-f2)*d, g1=(1-f2)*d, g2=(2-f2)*d;
				// output
				F[0] = g0 * g0; F[1] = g0 * g1; F[2] = 0; F[3] = g1 * g1;
				F[4] = g0 * g2; F[5] = g1 * g2; F[6] = g2 * g2; F[7] = 0;
			}
		} else {
			float *F = ((float *)p_G) + i*4*8;
			for (size_t k=0; k < 4; k++, F+=8)
			{
				double f = (n[k] > 0) ? double(sum[k])/(2*n[k]) : R_NaN;
				double d = 1 / sqrt(2*f*(1-f));
				if (!R_FINITE(d)) d = f = 0;
				const double f2 = 2*f;
				// normalized genotypes
				const double g0=(0-f2)*d, g1=(1-f2)*d, g2=(2-f2)*d;
				// output
				F[0] = g0 * g0; F[1] = g0 * g1; F[2] = 0; F[3] = g1 * g1;
				F[4] = g0 * g2; F[5] = g1 * g2; F[6] = g2 * g2; F[7] = 0;
			}
		}
	}
	PARALLEL_END
	PARALLEL_THREAD_BLOCK_END
}


/// Calculate GRM in the block (i,j)
static void grm_sp_calc_block(const double rel,
	const int i_st, const int i_n, const int j_st, const int j_n,
	const size_t n_byte, const RawMatrix &G_Pack, const NumericMatrix &G_Lookup,
	PARALLEL_VECTOR<t_sp_i_j> &sp)
{
	// for each block (i_n, j_n)
	PARALLEL_FOR(i_b, i_n, true)
	{
		const size_t nrow = G_Pack.nrow();
		const Rbyte *p_i = &G_Pack[(i_b + i_st)*nrow];
		const Rbyte *p_j = &G_Pack[(0   + j_st)*nrow];
		// look-up table
		const USE_FLOAT_TYPE *p_G = (const USE_FLOAT_TYPE*)&G_Lookup[0];
		for (int j_b=0; j_b < j_n; j_b++, p_j += nrow)
		{
			if (j_b+j_st >= (int)i_b+i_st)
			{
				int miss_n = 0;
				double sum = 0;
				// for-each SNP (*p_i, *p_j)
				(*fc_grm_calc_update_f32)(p_i, p_j, n_byte, p_G, miss_n, sum);
				// set the relatedness in the block
				const int sum_n = 4*n_byte - miss_n;
				double v = (sum_n > 0) ? (sum / sum_n) : 0;
				if (v >= rel)
					sp.push_back(t_sp_i_j(i_b+i_st, j_b+j_st));
			}
		}
	}
	PARALLEL_END
}


/// Create a sparse GRM by return a list of (i,j,x)
RcppExport SEXP saige_grm_sp_calc(SEXP nVariant, SEXP g_pack, SEXP g_lookup,
	SEXP rel_cutoff, SEXP bl_size, SEXP prog, SEXP prog_func)
{
BEGIN_RCPP
	// initialize values
	const int n_snp = Rf_asInteger(nVariant);
	RawMatrix G_Pack(g_pack);
	const int nSamp = G_Pack.ncol();
	NumericMatrix G_Lookup(g_lookup);
	const double rel = Rf_asReal(rel_cutoff);
	const int bs = Rf_asInteger(bl_size);
	Function prog_fc_r(prog_func);
	const bool verbose = !Rf_isNull(prog);

	// number of threading
	if (SAIGE_NumThread > nSamp) SAIGE_NumThread = nSamp;

	// fill g_lookup using float or double
	grm_sp_init_lookup(g_pack, g_lookup, USE_FLOAT_BOOL);
	PARALLEL_VECTOR<t_sp_i_j> sp;
	sp.reserve(nSamp*4);

	PARALLEL_THREAD_BLOCK
	// get numbers
	const size_t n_byte = (n_snp / 4) + (n_snp % 4 > 0);
	const int n_block = (nSamp / bs) + (nSamp % bs > 0);
	// for-loop blocks
	for (int i=0; i < n_block; i++)
	{
		// block i
		const int i_st = i * bs;
		const int i_n  = std::min(i_st+bs, nSamp) - i_st;
		for (int j=i; j < n_block; j++)
		{
			// block j
			const int j_st = j * bs;
			const int j_n  = std::min(j_st+bs, nSamp) - j_st;
			// calculation in the block (i,j)
			grm_sp_calc_block(rel, i_st, i_n, j_st, j_n, n_byte,
				G_Pack, G_Lookup, sp);
			// update progress
			if (verbose) prog_fc_r(prog, 1);
		}
	}
	PARALLEL_THREAD_BLOCK_END

	// output
	const size_t n = sp.size();
	IntegerVector r_i(n), r_j(n);
	for (size_t i=0; i < n; i++)
	{
		t_sp_i_j &v = sp[i];
		r_i[i] = v.i; r_j[i] = v.j;
	}
	return List::create(_["i"] = r_i, _["j"] = r_j);
END_RCPP
}


/// Update the non-zero entries (i,j,x) in GRM using the full variant set
RcppExport SEXP saige_grm_sp_calc_ijx(SEXP I, SEXP J, SEXP nVariant,
	SEXP g_pack, SEXP g_lookup, SEXP bl_size, SEXP prog, SEXP prog_func)
{
BEGIN_RCPP
	// initialize values
	const int n_snp = Rf_asInteger(nVariant);
	RawMatrix G_Pack(g_pack);
	NumericMatrix G_Lookup(g_lookup);
	const size_t bs = Rf_asInteger(bl_size);
	Function prog_fc_r(prog_func);
	const bool verbose = !Rf_isNull(prog);

	// number of non-zero entries
	const size_t nnzero = Rf_xlength(I);
	// number of blocks
	const size_t n_block = (nnzero / bs) + ((nnzero % bs) ? 1 : 0);
	// number of threading
	if (SAIGE_NumThread > n_block) SAIGE_NumThread = n_block;
	// fill g_lookup
	grm_sp_init_lookup(g_pack, g_lookup, true);
	// resulting numeric vector
	NumericVector X(nnzero);

	PARALLEL_THREAD_BLOCK
	// sparse matrix structure
	const int *ptrI = INTEGER(I);
	const int *ptrJ = INTEGER(J);
	double *ptrX = REAL(X);
	const double *p_G = &G_Lookup[0];  // look-up table
	const size_t nrow = G_Pack.nrow();
	const size_t n_byte = (n_snp / 4) + (n_snp % 4 > 0);
	// for-loop blocks
	for (size_t i_b=0; i_b < n_block; i_b++)
	{
		// block i
		const size_t i_st = i_b * bs;
		const size_t i_n  = std::min(bs, nnzero-i_st);
		// parallel for-loop
		PARALLEL_FOR(i, i_n, false)
		{
			const size_t k = i_st + i;
			const Rbyte *p_i = &G_Pack[ptrI[k]*nrow];
			const Rbyte *p_j = &G_Pack[ptrJ[k]*nrow];
			int miss_n = 0;
			double sum = 0;
			(*fc_grm_calc_update_f64)(p_i, p_j, n_byte, p_G, miss_n, sum);
			// set the relatedness in the block
			const int sum_n = 4*n_byte - miss_n;
			ptrX[k] = sum / sum_n;
		}
		PARALLEL_END
		// update progress
		if (verbose) prog_fc_r(prog, 1);
	}
	PARALLEL_THREAD_BLOCK_END

	// no output
	return X;
END_RCPP
}


// ========================================================================= //

/// Create a sparse GRM
RcppExport SEXP saige_grm_ds_calc(SEXP nVariant, SEXP g_pack, SEXP g_lookup,
	SEXP use_double, SEXP bl_size, SEXP prog, SEXP prog_func)
{
BEGIN_RCPP
	// initialize values
	const int n_snp = Rf_asInteger(nVariant);
	RawMatrix G_Pack(g_pack);
	const int nSamp = G_Pack.ncol();
	NumericMatrix G_Lookup(g_lookup);
	const bool use_f64 = Rf_asLogical(use_double) == TRUE;
	const int bs = Rf_asInteger(bl_size);
	Function prog_fc_r(prog_func);
	const bool verbose = !Rf_isNull(prog);

	// number of threading
	if (SAIGE_NumThread > nSamp) SAIGE_NumThread = nSamp;

	// fill g_lookup
	grm_sp_init_lookup(g_pack, g_lookup, use_f64);
	// resulting matrix
	NumericMatrix grm(nSamp, nSamp);

	PARALLEL_THREAD_BLOCK
	// get numbers
	const size_t n_byte = (n_snp / 4) + (n_snp % 4 > 0);
	const int n_block = (nSamp / bs) + (nSamp % bs > 0);
	// for-loop blocks
	for (int i=0; i < n_block; i++)
	{
		// block i
		const int i_st = i * bs;
		const int i_n  = std::min(i_st+bs, nSamp) - i_st;
		for (int j=i; j < n_block; j++)
		{
			// block j
			const int j_st = j * bs;
			const int j_n  = std::min(j_st+bs, nSamp) - j_st;
			// calculation in the block (i,j)
			PARALLEL_FOR(i_b, i_n, true)
			{
				const size_t nrow = G_Pack.nrow();
				const Rbyte *p_i = &G_Pack[(i_b + i_st)*nrow];
				const Rbyte *p_j = &G_Pack[(0   + j_st)*nrow];
				const void *p_G = &G_Lookup[0];  // look-up table
				const int m_i = i_b + i_st;
				for (int j_b=0; j_b < j_n; j_b++, p_j += nrow)
				{
					const int m_j = j_b + j_st;
					if (m_j >= m_i)
					{
						int miss_n = 0;
						double sum = 0;
						// for-each SNP (*p_i, *p_j)
						if (use_f64)
						{
							(*fc_grm_calc_update_f64)(p_i, p_j, n_byte,
								(const double *)p_G, miss_n, sum);
						} else {
							(*fc_grm_calc_update_f32)(p_i, p_j, n_byte,
								(const float *)p_G, miss_n, sum);
						}
						// set the relatedness in the block
						const int sum_n = 4*n_byte - miss_n;
						grm(m_i, m_j) = grm(m_j, m_i) =
							(sum_n > 0) ? (sum / sum_n) : R_NaN;
					}
				}
			}
			PARALLEL_END
			// update progress
			if (verbose) prog_fc_r(prog, 1);
		}
	}
	PARALLEL_THREAD_BLOCK_END

	// output
	return grm;
END_RCPP
}
