// ===========================================================
//
// saige_main.cpp: SAIGE association analysis
//
// Copyright (C) 2019-2022    Xiuwen Zheng / AbbVie-ComputationalGenomics
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

// To ensure atanpi, cospi, sinpi, tanpi are defined, used in ACAT calculation
#ifndef __STDC_WANT_IEC_60559_FUNCS_EXT__
#   define __STDC_WANT_IEC_60559_FUNCS_EXT__    1
#endif

#include "vectorization.h"
#include "vec_ext.h"
#include <RcppArmadillo.h>
#include "saige.h"
#include <Rconfig.h>
#include <Rmath.h>

#include <math.h>
#include <vector>
#include <algorithm>


using namespace Rcpp;
using namespace vectorization;
using namespace arma;
using namespace SAIGE;


// ========================================================================= //

// disable timing
// #define TIMING

#ifdef TIMING

#include <chrono>

using namespace std::chrono;

struct auto_timing
{
	auto_timing(nanoseconds &tm): time(tm)
		{ tp = steady_clock::now(); }
	~auto_timing()
		{ time += duration_cast<nanoseconds>(steady_clock::now() - tp); }
private:
	steady_clock::time_point tp;
	nanoseconds &time;
};

static const int n_run_time = 4 + 5;
static const int n_st_sv = 0;
static const int n_st_skat = 4;
static nanoseconds run_time[n_run_time];

static inline void init_run_time()
{
	for (int i=0; i < n_run_time; i++)
		run_time[i] = nanoseconds::zero();
}

RcppExport SEXP saige_timing()
{
BEGIN_RCPP
	NumericVector rv(n_run_time);
	for (int i=0;i < n_run_time; i++)
		rv[i] = run_time[i].count() / 1000000000.0;
	return rv;
END_RCPP
}

#endif


// ========================================================================= //
// Internal functions

// whether always use the fast SPA or not
static bool SPA_always_use_fastSPA = true;
// the default Cutoff value used in SPAtest
static const double SPA_default_cutoff = 2;
// the default p-value cutoff according to SPA_default_cutoff
static const double SPA_default_pval_cutoff = 0.05;

// ==== saige_misc.cpp ====
namespace Misc
{
	extern int SummaryStat_Mat(SEXP mat, double out_af[], double out_mac[]);
	extern int SummaryStat_SpMat(SEXP mat, double out_af[], double out_mac[]);
	extern sp_mat GetSp_Impute_SpMat(SEXP mat, double af[], double mac[],
		double mac_imp[]);
	extern sp_mat GetSp_CollapseGenoMat(const sp_mat &mat, double collapse_mac,
		int collapse_method,
		const double mac[], double inout_maf[], int &out_n_collapse);
}


/// minor allele frequency
inline double MAF(double af) { return std::min(af, 1-af); }


// ========================================================================= //

static TTrait mod_trait = TTrait::Unknown;  //< outcome type (Quant, ...)

static double threshold_maf = 0;  //< the threshold of MAF filter
static double threshold_mac = 0;  //< the threshold of MAC filter
static double threshold_missing = 1;      //< the threshold of missing proportion per variant
static double threshold_pval_spa = 0.05;  //< the threshold of p-value filter for SPA

static int mod_NSamp = 0;   //< the number of samples
static int mod_NCoeff = 0;  //< the number of beta coefficients

static int mod_ploidy = 2;  //< ploidy or <=0 if it is not a genotype

                                         // K is # of beta coefficients
static const double *mod_tau = NULL;     //< variance components: tau[0], tau[1]
static const double *mod_y = NULL;       //< a n_samp-length vector
static const double *mod_mu = NULL;      //< a n_samp-length vector
static const double *mod_y_mu = NULL;    //< a n_samp-length vector, y-mu (residuals)
static const double *mod_mu2 = NULL;     //< a n_samp-length vector, mu*(1-mu) for binary outcomes
static const double *mod_t_XXVX_inv = NULL;    //< a K-by-n_samp matrix
static const double *mod_XV = NULL;      //< a K-by-n_samp matrix
static const double *mod_t_XVX_inv_XV = NULL;  //< a K-by-n_samp matrix
static const double *mod_XVX = NULL;           //< a K-by-K matrix
static const double *mod_t_X = NULL;           //< a K-by-n_samp matrix
static const double *mod_S_a = NULL;           //< a K-length vector

static const double *mod_Si_X = NULL;               //< a K-by-n_samp matrix
static const double *mod_XVX_inv_XV_X_Si_X = NULL;  //< a K-by-n_samp matrix

static double *mod_varRatio = NULL;         //< the variance ratio; NaN for not using variance ratio
static double *mod_varRatioSqrt = NULL;     //< sqrt(mod_varRatio)
static double *mod_varRatioInvSqrt = NULL;  //< 1/sqrt(mod_varRatio)
static double *mod_varRatioMAC = NULL;      //< MAC categories
static int mod_varRatioMACLen = 0;          //< length(mod_varRatioMAC)

static const double *mod_sigma_inv_val = NULL;  //< the inverse of Sigma matrix
static const int *mod_sigma_inv_i = NULL;       //< if sparse matrix (dsCMatrix)
static const int *mod_sigma_inv_p = NULL;       //< if sparse matrix (dsCMatrix)
static const double *mod_chol_inv_X_Sigma = NULL;

static double *buf_dosage = NULL;    //< temporary buffer for real dosages
static double *buf_coeff  = NULL;     // beta coefficients
static double *buf_coeff2 = NULL;     // beta coefficients
static double *buf_adj_g = NULL;     //< genotype after adjusting for fixed effects
static int *buf_index = NULL;
static double *buf_B = NULL;
static double *buf_g_tilde = NULL;   //< est. adj. genotypes
static double *buf_X1 = NULL;        //< ncol(X1)
static double *buf_spa = NULL;       //< buffer for SPA calculation

// aggregate tests
static double threshold_summac = 0;  //< the threshold of weighted sum MAC
static int num_wbeta = 0;            //< # of beta parameters
static double *buf_wbeta = NULL;     //< beta parameters
static int num_unitsz = 0;           //< max unit size
static double *buf_unitsz = NULL;    //< length() = max unit size

/// the MAC threshold for collapsing ultra rare variants in ACAT-V
static double threshold_acatv_mac = 0;
/// the MAC threshold for collapsing ultra rare variants in SKAT
static double threshold_skat_mac = 0;

// the index of non-zero stored in buf_index
#define IDX_i    buf_index[i]


/// initialize internal parameters from the model object
RcppExport SEXP saige_score_test_init(SEXP model)
{
BEGIN_RCPP

	List M(model);

	// trait
	mod_trait = (TTrait)Rf_asInteger(M["trait"]);
	if (mod_trait!=TTrait::Quant && mod_trait!=TTrait::Binary)
		Rf_error("Invalid trait index: %d.", (int)mod_trait);

	// threshold setting
	threshold_maf = Rf_asReal(M["maf"]);
	if (!R_FINITE(threshold_maf)) threshold_maf = -1;
	threshold_mac = Rf_asReal(M["mac"]);
	if (!R_FINITE(threshold_mac)) threshold_mac = -1;
	threshold_missing = Rf_asReal(M["missing"]);
	if (!R_FINITE(threshold_missing)) threshold_missing = 1;
	threshold_pval_spa = Rf_asReal(M["spa.pval"]);
	if (!R_FINITE(threshold_pval_spa))
		threshold_pval_spa = SPA_default_pval_cutoff;

	// model parameters
	mod_NSamp = Rf_length(M["y"]);
	mod_NCoeff = NumericMatrix(wrap(M["XV"])).nrow();
	mod_ploidy = Rf_asInteger(M["geno.ploidy"]);
	mod_tau = REAL(M["tau"]);
	mod_y = REAL(M["y"]);
	mod_mu = REAL(M["mu"]);
	mod_y_mu = REAL(M["y_mu"]);
	mod_mu2 = REAL(M["mu2"]);
	mod_t_XXVX_inv = REAL(M["t_XXVX_inv"]);
	mod_XV = REAL(M["XV"]);
	mod_t_XVX_inv_XV = REAL(M["t_XVX_inv_XV"]);
	mod_XVX = REAL(M["XVX"]);
	mod_t_X = REAL(M["t_X"]);
	mod_S_a = REAL(M["S_a"]);

	mod_varRatio = REAL(M["var.ratio"]);
	mod_varRatioSqrt = REAL(M["vr_sqrt"]);
	mod_varRatioInvSqrt = REAL(M["vr_inv_sqrt"]);
	mod_varRatioMAC = REAL(M["vr_cateMAC"]);
	mod_varRatioMACLen = Rf_length(M["vr_cateMAC"]);
	if (Rf_length(M["var.ratio"]) != mod_varRatioMACLen+1)
		Rf_error("Length of var.ratio is invalid.");

	SEXP sigma = M["Sigma_inv"];
	if (!Rf_isNull(sigma))
	{
		mod_Si_X = REAL(M["Si_X"]);
		mod_XVX_inv_XV_X_Si_X = REAL(M["XVX_inv_XV_X_Si_X"]);
		// save Sigma inv
		if (Rf_isMatrix(sigma))
		{
			mod_sigma_inv_val = REAL(sigma);
			mod_sigma_inv_i = NULL;
			mod_sigma_inv_p = NULL;
		} else {
			if (!Rf_inherits(sigma, "dsCMatrix"))
				Rf_error("'Sigma_inv' should be a dsCMatrix matrix.");
			S4 mat(sigma);
			mod_sigma_inv_i = INTEGER(mat.slot("i"));
			mod_sigma_inv_p = INTEGER(mat.slot("p"));
			mod_sigma_inv_val = REAL(mat.slot("x"));
		}
		if (!Rf_isMatrix(M["chol_inv_X_Sigma"]))
			Rf_error("'chol_inv_X_Sigma' should be a matrix.");
		mod_chol_inv_X_Sigma = REAL(M["chol_inv_X_Sigma"]);
	} else {
		mod_Si_X = mod_XVX_inv_XV_X_Si_X = NULL;
		mod_sigma_inv_val = NULL;
		mod_sigma_inv_i = mod_sigma_inv_p = NULL;
		mod_chol_inv_X_Sigma = NULL;
	}

	// buffer
	buf_dosage = REAL(M["buf_dosage"]);
	buf_coeff  = REAL(M["buf_coeff"]);
	buf_coeff2 = REAL(M["buf_coeff"]) + mod_NCoeff;
	buf_adj_g = REAL(M["buf_adj_g"]);
	buf_index = INTEGER(M["buf_index"]);
	buf_B = REAL(M["buf_B"]);
	buf_g_tilde = REAL(M["buf_g_tilde"]);
	buf_X1 = REAL(M["buf_X1"]);
	buf_spa = REAL(M["buf_spa"]);

	// buffer for aggregate tests
	threshold_summac = Rf_asReal(M["summac"]);
	if (!R_FINITE(threshold_summac)) threshold_summac = -1;
	threshold_acatv_mac = Rf_asReal(M["acatv_mac"]);
	if (!R_FINITE(threshold_acatv_mac)) threshold_acatv_mac = 10;
	threshold_skat_mac = Rf_asReal(M["skat_mac"]);
	if (!R_FINITE(threshold_skat_mac))
		Rf_error("'skat.collapse.mac' should be a finite number.");

	num_wbeta = Rf_length(M["buf_wbeta"]) / 2;  // # of columns
	buf_wbeta = REAL(M["buf_wbeta"]);
	num_unitsz = M["num_unitsz"];
	buf_unitsz = REAL(M["buf_unitsz"]);

#ifdef TIMING
	init_run_time();
#endif

END_RCPP
}



// ========================================================================= //
// Calculate the variance of Tstat

inline static int get_var_ratio_index(double mac)
{
	for (int i=0; i < mod_varRatioMACLen; i++)
	{
		if (mac < mod_varRatioMAC[i])
			return i;
	}
	return mod_varRatioMACLen;
}


// ========================================================================= //
// Calculate the variance of Tstat

// G -- raw genotype vector (0,1,2), get G' Sigma_inv G
MATH_OFAST static double calc_G_Si_G(
	size_t nnzero, int nnz_idx[], const double G[])
{
#ifdef TIMING
	auto_timing tm3(run_time[n_st_sv+3]);
#endif
	const size_t n = mod_NSamp;
	double var_s = 0;
	if (mod_sigma_inv_i)
	{
		// sparse matrix, it should be dsCMatrix
		for (size_t i=0; i < nnzero; i++)
		{
			const int I = nnz_idx[i];
			const double G_I = G[I];
			const size_t st = mod_sigma_inv_p[I];
			const size_t ed = mod_sigma_inv_p[I+1];
			// much fast if only a few nonzero values in the column I
			for (size_t j=st; j < ed; j++)
			{
				int j_i = mod_sigma_inv_i[j];
				double v = G_I * G[j_i] * mod_sigma_inv_val[j];
				var_s += (I != j_i) ? (v+v) : v;
			}
		}
	} else if (mod_sigma_inv_val)
	{
		// dense matrix
		for (size_t i=0; i < nnzero; i++)
		{
			const size_t I = nnz_idx[i];
			const double *p = &mod_sigma_inv_val[n*I];
			const double G_i = G[I];
			var_s += G_i * G_i * p[I];
			const double G2_i = G_i + G_i;
			for (size_t j=i+1; j < nnzero; j++)
			{
				const size_t J = nnz_idx[j];
				var_s += G2_i * G[J] * p[J];
			}
		}
	}
	// output
	return var_s;
}

// G -- adjusted genotypes (G tilde), get G' Sigma_inv G
COREARRAY_TARGET_CLONES MATH_OFAST
	static double calc_varT(size_t nnzero, int nnz_idx[], const double G[])
{
	double var_s = calc_G_Si_G(nnzero, nnz_idx, G);

	// buf_coeff = XVX_inv_XV * G
	f64_mul_mat_vec_sp(nnzero, buf_index, mod_NCoeff, mod_t_XVX_inv_XV, G,
		buf_coeff);

	// buf_coeff2 = Si_X * G
	f64_mul_mat_vec_sp(nnzero, buf_index, mod_NCoeff, mod_Si_X, G,
		buf_coeff2);
	double s = 0;
	for (int i=0; i < mod_NCoeff; i++) s += buf_coeff[i] * buf_coeff2[i];
	var_s -= s + s;

	// buf_coeff2 = XVX_inv_XV_X_Si_X * G
	f64_mul_mat_vec_sp(nnzero, buf_index, mod_NCoeff, mod_XVX_inv_XV_X_Si_X, G,
		buf_coeff2);
	s = 0;
	for (int i=0; i < mod_NCoeff; i++) s += buf_coeff[i] * buf_coeff2[i];
	var_s += s;

	// output
	return var_s;
}


// Gmat -- raw genotype matrix (0,1,2), get Gmat' Sigma_inv Gmat
MATH_OFAST static void calc_mat_G_Si_G(
	const sp_mat &Gmat, const Type_dgCMatrix &Si, dmat &out_GPG)
{
	const size_t ncol = Gmat.n_cols;
	out_GPG.set_size(ncol, ncol);
	dvec G0(Gmat.n_rows);
	double *G0_p = &G0[0];
	const size_t g0_sz = sizeof(double)*Gmat.n_rows;
	// sort by the number of non-zero entries, decreasing
	// faster calculation after sorting
	std::vector< std::pair<int, int> > Idx(ncol);
	for (size_t i=0; i < ncol; i++)
	{
		std::pair<int, int> &v = Idx[i];
		v.first = -(Gmat.end_col(i).pos() - Gmat.begin_col(i).pos());
		v.second = i;
	}
	std::sort(Idx.begin(), Idx.end());
	// for each pair <i,j>
	for (size_t i=0; i < ncol; i++)
	{
		const size_t I = Idx[i].second;
		// G0 = Gmat.col(i)
		memset(G0_p, 0, g0_sz);
		sp_mat::const_iterator i_it = Gmat.begin_col(I);
		sp_mat::const_iterator i_ed = Gmat.end_col(I);
		for (; i_it != i_ed; ++i_it)
			G0_p[i_it.row()] = *i_it;
		// pair(i, j)
		for (size_t j=i; j < ncol; j++)
		{
			const size_t J = Idx[j].second;
			sp_mat::const_iterator it = Gmat.begin_col(J);
			sp_mat::const_iterator ed = Gmat.end_col(J);
			double sum = 0;
			for (; it != ed; ++it)
			{
				// much fast if only a few nonzero values in the row/column k (Si is symmetric)
				const size_t k = it.row();
				const size_t k_ed = Si.p[k+1];
				double s = 0;
				for (size_t h=Si.p[k]; h < k_ed; h++)
					s += G0_p[Si.i[h]] * Si.x[h];
				sum += (*it) * s;
			}
			out_GPG.at(I,J) = out_GPG.at(J,I) = sum;
		}
	}
}


// ========================================================================= //

static const char *ERR_DS_TYPE = "Invalid type of dosages.";
static const char *ERR_DS_LEN  = "Invalid length of dosages: %d.";

/// get numeric dosages from an R object
static double *get_ds(SEXP ds, R_xlen_t n, R_xlen_t start, double *ds_buf=NULL)
{
	const R_xlen_t ntot = Rf_xlength(ds);
	if (start < 0) start = 0;
	if (n - start > ntot) Rf_error(ERR_DS_LEN, ntot);

	R_xlen_t i = 0;
	switch (TYPEOF(ds))
	{
	case REALSXP:
		if (ds_buf)
		{
			memcpy(ds_buf, REAL(ds) + start, sizeof(double)*n);
			return ds_buf;
		} else
			return REAL(ds) + start;
	case INTSXP:
		if (!ds_buf) ds_buf = buf_dosage;
		for (const int *p = INTEGER(ds) + start; i < n; i++)
			ds_buf[i] = (p[i] != NA_INTEGER) ? p[i] : R_NaN;
		return ds_buf;
	case RAWSXP:
		if (!ds_buf) ds_buf = buf_dosage;
		for (const Rbyte *p = RAW(ds) + start; i < n; i++)
			ds_buf[i] = (p[i] != Rbyte(0xFF)) ? p[i] : R_NaN;
		return ds_buf;
	}
	Rf_error(ERR_DS_TYPE);
	return NULL;
}


// ====================================

/// single variant test with score statistics
///   assuming no missing genotype in G and AF <= 0.5 when mod_ploidy > 0
static size_t g_score_test(const double G[], double mac,
	double *out_beta, double *out_SE, double *out_pval, double *out_pval_noadj,
	bool *out_converged, double *Tstat, bool always_use_fast)
{
	// get the number of nonzeros and the nonzero indices
	const size_t nnzero = f64_nonzero_index(mod_NSamp, &G[0], buf_index);

	// get MAC and variance ratio index
	if (mod_ploidy > 0)
	{
		if (!R_FINITE(mac))
		{
			mac = 0;
			for (size_t i=0; i < nnzero; i++) mac += G[buf_index[i]];
		}
	} else {
		mac = nnzero;
	}
	const int v_idx = get_var_ratio_index(mac);
	const double c_varRatio = (mac>=0 ? mod_varRatio[v_idx] : 1);
	const double c_varRatioInvSqrt = (mac>=0 ? mod_varRatioInvSqrt[v_idx] : 1);

	// buf_coeff = XVX_inv_XV * G
	f64_mul_mat_vec_sp(nnzero, buf_index, mod_NCoeff, mod_t_XVX_inv_XV,
		&G[0], buf_coeff);
	// buf_B = t(X) * buf_coeff
	f64_mul_mat_vec_sub(nnzero, buf_index, mod_NCoeff, mod_t_X, buf_coeff,
		buf_B);
	// g_tilde = G - B
	for (size_t i=0; i < nnzero; i++)
		buf_g_tilde[i] = G[buf_index[i]] - buf_B[i];

	// for Tstat
	// S1 = sum(y_mu .* g_tilde)
	double S1 = 0;
	for (size_t i=0; i < nnzero; i++)
		S1 += mod_y_mu[buf_index[i]] * buf_g_tilde[i];
	// buf_X1 = t(X1) * (y-mu)
	f64_mul_mat_vec_sp(nnzero, buf_index, mod_NCoeff, mod_t_X,
		mod_y_mu, buf_X1);
	// S2 = sum((buf_X1 - mod_S_a) .* buf_coeff)
	double S2 = 0;
	for (int i=0; i < mod_NCoeff; i++)
		S2 += (buf_X1[i] - mod_S_a[i]) * buf_coeff[i];
	double S = S1 + S2;
	if (Tstat) *Tstat = S;

	// for SE
	double var2;
	if (!mod_sigma_inv_val)
	{
		// t(buf_coeff) %*% XVX %*% buf_coeff
		var2 = f64_sum_mat_vec(mod_NCoeff, mod_XVX, buf_coeff);
		// no user-defined GRM
		if (mod_trait == TTrait::Binary)
		{
			// var2 = t(buf_coeff) %*% XVX %*% buf_coeff - sum(B^2 .* mu2) +
			//        sum(g_tilde^2 .* mu2)
			for (size_t i=0; i < nnzero; i++)
				var2 += (sq(buf_g_tilde[i]) - sq(buf_B[i])) * mod_mu2[IDX_i];
		} else {
			// var2 = t(buf_coeff) %*% XVX %*% buf_coeff - sum(B^2) +
			//        sum(g_tilde^2)
			for (size_t i=0; i < nnzero; i++)
				var2 += sq(buf_g_tilde[i]) - sq(buf_B[i]);
			var2 *= mod_tau[0];
		}
	} else {
	#ifdef TIMING
		auto_timing tm2(run_time[n_st_sv+2]);
	#endif
		// use user-defined GRM
		// variance according to Tstat (needing tau[0])
		var2 = calc_varT(nnzero, buf_index, G) * sq(mod_tau[0]);
	}
	const double var2_u_grm = var2;
	const double var1 = var2 * c_varRatio;
	double beta = (mod_trait == TTrait::Binary) ?
		(S / var1) : (S / var1 * mod_tau[0]);

	double pval_noadj = ::Rf_pchisq(S*S/var1, 1, FALSE, FALSE);
	double pval = pval_noadj;
	bool converged = R_FINITE(pval_noadj) != 0;

	// need further SPAtest or not, if binary outcome
	if ((mod_trait == TTrait::Binary) &&
		converged && (pval_noadj <= threshold_pval_spa))
	{
		// need adjusted genotypes, adj_g = G - XXVX_inv * (XV * G)
		// buf_coeff = XV * G
		f64_mul_mat_vec_sp(nnzero, buf_index, mod_NCoeff, mod_XV, &G[0],
			buf_coeff);
		// buf_adj_g = G - XXVX_inv * buf_coeff
		f64_sub_mul_mat_vec(mod_NSamp, mod_NCoeff, &G[0], mod_t_XXVX_inv,
			buf_coeff, buf_adj_g);
		// q = sum(y .* adj_g)
		double q = f64_dot(mod_NSamp, mod_y, buf_adj_g);
		// m1 = sum(mu .* adj_g)
		// var2 = sum(mu*(1-mu) .* adj_g .* adj_g)
		double m1, var2;
		f64_dot_sp2(mod_NSamp, mod_mu, mod_mu2, buf_adj_g, m1, var2);
		// var with mixed effects
		double var1, vr_inv_sqrt;
		if (!mod_sigma_inv_val)
		{
			var1 = var2 * c_varRatio;
			vr_inv_sqrt = c_varRatioInvSqrt;
		} else {
			var1 = var2_u_grm * c_varRatio;
			vr_inv_sqrt = sqrt(var2 / var1);
		}

		double Tstat = q - m1;
		double qtilde = Tstat * vr_inv_sqrt + m1;

	#ifdef TIMING
		auto_timing tm1(run_time[n_st_sv+1]);
	#endif
		// call Saddle_Prob in SPAtest
		if (always_use_fast || (2*nnzero <= (size_t)mod_NSamp))
		{
			pval = Saddle_Prob_Fast(qtilde, m1, var2, mod_NSamp, mod_mu,
				buf_adj_g, nnzero, buf_index, SPA_default_cutoff,
				converged, buf_spa, NULL);
		} else {
			pval = Saddle_Prob(qtilde, m1, var2, mod_NSamp, mod_mu,
				buf_adj_g, SPA_default_cutoff, converged, NULL);
		}
		if (pval==0 && pval_noadj>0)
		{
			pval = pval_noadj;
			converged = false;
		}

		// effect size
		beta = Tstat / var1;
	}

	double SE = fabs(beta / ::Rf_qnorm5(pval/2, 0, 1, TRUE, FALSE));

	// output
	if (out_beta) *out_beta = beta;
	if (out_SE) *out_SE = SE;
	if (out_pval) *out_pval = pval;
	if (out_pval_noadj) *out_pval_noadj = pval_noadj;
	if (out_converged) *out_converged = converged;
	return nnzero;
}

/// single variant test with score statistics
static bool single_score_test(double G[],
	double &oAF, double &omac, int &onum, double &obeta, double &oSE,
	double &opval, double &opval_noadj, bool &oconverged, double *Tstat=NULL)
{
	// calculate allele freq, and impute geno using the mean
	double AF, AC;
	int Num;  // # of valid SNVs
	f64_af_ac_impute(&G[0], mod_NSamp, AF, AC, Num, buf_index, mod_ploidy);
	// check thresholds
	const double maf =
		(mod_ploidy > 0) ? std::min(AF, 1-AF) : AF;
	const double mac =
		(mod_ploidy > 0) ? std::min(AC, mod_ploidy*Num - AC) : AC;
	const double missing =
		double(mod_NSamp - Num) / mod_NSamp;
	const bool b = (mod_ploidy > 0) ?
		((maf > 0) && (maf >= threshold_maf) && (mac >= threshold_mac)) : true;

	if ((Num > 0) && b && (missing <= threshold_missing))
	{
		const bool minus = (mod_ploidy>0) ? (AF > 0.5) : false;
		if (mod_ploidy>0 && minus)
			f64_sub(mod_NSamp, mod_ploidy, &G[0]);
		// MAC after mean imputation
		const double mac_g = (Num < mod_NSamp) ? (2 * mod_NSamp * maf) : mac;
		// output
		oAF = AF; omac = mac; onum = Num;
		size_t nnz = g_score_test(G, mac_g, &obeta, &oSE, &opval, &opval_noadj,
			&oconverged, Tstat, SPA_always_use_fastSPA);
		if (minus) obeta = -obeta;
		// nnzero if not genotype
		if (mod_ploidy <= 0)
		{
			if (nnz == 0) return false;
			omac = nnz;
		}
		return true;
	} else
		return false;
}


// ====================================

/// calculate single-variant p-values for binary outcomes
RcppExport SEXP saige_score_test_pval(SEXP dosage)
{
BEGIN_RCPP
#ifdef TIMING
	auto_timing tm(run_time[n_st_sv+0]);
#endif

	// get numeric dosage
	if (Rf_xlength(dosage) != mod_NSamp)
		Rf_error(ERR_DS_LEN, Rf_xlength(dosage));
	double *G = get_ds(dosage, mod_NSamp, 0);

	int num = 0;
	bool converged = false;
	double AF, mac, beta, SE, pval, pval_noadj;
	AF = mac = beta = SE = pval = pval_noadj = R_NaN;

	if (single_score_test(G, AF, mac, num, beta, SE, pval, pval_noadj,
		converged))
	{
		const bool binary = (mod_trait == TTrait::Binary);
		// output
		SEXP rv_ans = NEW_NUMERIC(binary ? 8 : 6);
		double *ans = REAL(rv_ans);
		ans[0] = AF;    ans[1] = mac;   ans[2] = num;
		ans[3] = beta;  ans[4] = SE;    ans[5] = pval;
		if (binary)
		{
			ans[6] = pval_noadj;
			ans[7] = converged ? 1 : 0;
		}
		return rv_ans;
	} else
		return R_NilValue;
END_RCPP
}


// ========================================================================= //
// Aggregate Tests

static const int AGGR_INDEX_START = 9;

static void summary_maf_mac(NumericVector &ans, int n_snv,
	const double maf[], const double mac[])
{
	ans[0] = n_snv;
	f64_mean_sd_maxmin(maf, n_snv, ans[1], ans[2], ans[4], ans[3]);
	f64_mean_sd_maxmin(mac, n_snv, ans[5], ans[6], ans[8], ans[7]);
}

static sp_mat get_G0_flipped_impute(SEXP dosage, double maf[],
	double mac[], double mac_imp[])
{
	if (Rf_isMatrix(dosage))
	{
		// it could be a RAW, integer or numeric matrix
		Misc::SummaryStat_Mat(dosage, maf, mac);
		Rf_error("not support!");
		return sp_mat();
	} else {
		// it should be a sparse matrix (dgCMatrix)
		// maf stores AF here for imputation
		Misc::SummaryStat_SpMat(dosage, maf, mac);
		// monomorphic variants are removed
		return Misc::GetSp_Impute_SpMat(dosage, maf, mac, mac_imp);
		// maf stores MAF since flipping
	}
}

inline static void add_g_w(dvec &G, const sp_mat &GMat, int col_i, double w)
{
	sp_mat::const_iterator it = GMat.begin_col(col_i);
	sp_mat::const_iterator ed = GMat.end_col(col_i);
	for (; it != ed; ++it)
		G[it.row()] += (*it) * w;
}

inline static double sum_col(const sp_mat &GMat, int col_i)
{
	sp_mat::const_iterator it = GMat.begin_col(col_i);
	sp_mat::const_iterator ed = GMat.end_col(col_i);
	double sum = 0;
	for (; it != ed; ++it) sum += *it;
	return sum;
}


// ====================================

static void gmat_burden_test(const sp_mat &G0, double beta_b1, double beta_b2,
	const double maf[], const double mac[], double w_burden[], double out_ans[])
{
	const int n_snv = G0.n_cols;
	int n_valid = 0;
	double summac = 0;  // the sum of MAC
	for (int i=0; i < n_snv; i++)
	{
		const double F = maf[i];
		if (R_FINITE(F) && (F > 0))
		{
			// set weights
			w_burden[i] = Rf_dbeta(F, beta_b1, beta_b2, FALSE);
			n_valid ++;
		} else
			w_burden[i] = R_NaN;
		if (R_FINITE(mac[i])) summac += mac[i];
	}

	// collapse SNVs with weights
	dvec G;
	G.zeros(mod_NSamp);
	f64_normalize(n_snv, w_burden);
	for (int i=0; i < n_snv; i++) add_g_w(G, G0, i, w_burden[i]);

	// p-value calculation
	bool converged = false;
	double beta, SE, pval, pval_noadj;
	beta = SE = pval = pval_noadj = R_NaN;
	if ((n_valid > 0) && (summac >= threshold_summac) && (summac > 0))
	{
		g_score_test(&G[0], summac, &beta, &SE, &pval, &pval_noadj,
			&converged, NULL, false);
	}

	// set the output
	out_ans[0] = summac;
	out_ans[1] = beta; out_ans[2] = SE;
	out_ans[3] = pval;
	if (mod_trait == TTrait::Binary)
	{
		out_ans[4] = pval_noadj;
		out_ans[5] = converged ? 1 : 0;
	}
}


/// calculate burden p-values
RcppExport SEXP saige_burden_test_pval(SEXP dosage)
{
BEGIN_RCPP
	// buffer
	double *maf = buf_unitsz;
	double *mac = buf_unitsz + num_unitsz;
	double *mac_imp = buf_unitsz + 2*num_unitsz;
	double *w_burden = buf_unitsz + 3*num_unitsz;

	// get genotype matrix
	sp_mat G0 = get_G0_flipped_impute(dosage, maf, mac, mac_imp);
	const int n_snv = G0.n_cols;
	const bool binary = (mod_trait == TTrait::Binary);
	const int n_each_wb = binary ? 5 : 3;

	// summarize maf & mac
	const int st_idx = AGGR_INDEX_START + 1;
	NumericVector ans(st_idx + n_each_wb*num_wbeta);
	summary_maf_mac(ans, n_snv, maf, mac);
	ans[AGGR_INDEX_START] = f64_sum(n_snv, mac);

	// for each beta weight
	for (int i=0; i < num_wbeta; i++)
	{
		// get weights, beta function
		const double b1 = buf_wbeta[2*i+0], b2 = buf_wbeta[2*i+1];
		// calculation
		double v[6];
		gmat_burden_test(G0, b1, b2, maf, mac, w_burden, v);
		// set the output
		const int k = st_idx + n_each_wb*i;
		if (i == 0) ans[AGGR_INDEX_START] = v[0];  // summac
		ans[k+0] = v[1];  // beta
		ans[k+1] = v[2];  // SE
		ans[k+2] = v[3];  // pval
		if (binary)
		{
			ans[k+3] = v[4];  // p.norm
			ans[k+4] = v[5];  // converged
		}
	}

	// output
	return ans;
END_RCPP
}


// ====================================

// Data structure for SKAT method
struct Struct_SKAT
{
	dmat Sigma_inv_d;
	Type_dgCMatrix Sigma_inv_s;
	bool is_sparse;
	dmat XVX_inv_XV, Si_X, XVX_inv_XV_X_Si_X;
	Function chisq_pval;
	// collapse_method = 1, presence or absence (PA)
	// collapse_method = 2, presence or absence (PA_int)
	// collapse_method = 3, sum up rare genotype (SumG)
	int collapse_method;

	// constructor
	Struct_SKAT(SEXP sigma_inv, SEXP m1, SEXP m2, SEXP m3, Function f, int cm):
		chisq_pval(f)
	{
		if (Rf_isNull(sigma_inv))
			Rf_error("Sigma_inv should not be NULL.");
		this->is_sparse = !Rf_isMatrix(sigma_inv);
		if (this->is_sparse)
			this->Sigma_inv_s.reset(sigma_inv);
		else
			this->Sigma_inv_d = as<dmat>(sigma_inv);
		this->XVX_inv_XV = as<dmat>(m1);
		this->Si_X = as<dmat>(m2);
		this->XVX_inv_XV_X_Si_X = as<dmat>(m3);
		this->collapse_method = cm;
	}
} *p_struct_skat = NULL;


RcppExport SEXP saige_skat_test_init(SEXP sigma_inv, SEXP XVX_inv_XV,
	SEXP Si_X, SEXP XVX_inv_XV_X_Si_X, SEXP collapse_method)
{
BEGIN_RCPP
	Environment pkg = Environment::namespace_env("SAIGEgds");
	Function f_r = pkg[".skat_eig_chiq"];
	// initialize Struct_SKAT
	p_struct_skat = new Struct_SKAT(sigma_inv, XVX_inv_XV, Si_X,
		XVX_inv_XV_X_Si_X, f_r, Rf_asInteger(collapse_method));
	return R_NilValue;
END_RCPP
}

RcppExport SEXP saige_skat_test_reset()
{
	p_struct_skat = NULL;
	return R_NilValue;
}

RcppExport SEXP saige_skat_test_done()
{
BEGIN_RCPP
	if (p_struct_skat)
	{
		delete p_struct_skat;
		p_struct_skat = NULL;
	}
	return R_NilValue;
END_RCPP
}

/// initialize Ts and GPG in SKAT
static void gmat_skat_test_p1(const sp_mat &G0, double var_ratio[],
	dvec &Ts, dmat &GPG)
{
#ifdef TIMING
	auto_timing tm1(run_time[n_st_skat+1]);
#endif
	const int ncol_g = G0.n_cols;
	Ts.set_size(ncol_g);
	// Tstat for each G_i in ds
	for (int i=0; i < ncol_g; i++)
	{
		const double mac = sum_col(G0, i);  // minor allele count
		const int vr_idx = get_var_ratio_index(mac);
		const double vr_sqrt = mod_varRatioSqrt[vr_idx];
		// genotypes
		dvec G(G0.col(i));
		Ts[i] = 0;  // initialize
		double pval_noadj=1, pval=pval_noadj;
		g_score_test(&G[0], mac, NULL, NULL, &pval, &pval_noadj, NULL,
			&Ts[i], false);
		// Ts[i] should get the correct value
		// get variance ratio from p-values for case-control imbalance
		var_ratio[i] = vr_sqrt;
		if (R_FINITE(pval_noadj) && R_FINITE(pval) && pval_noadj!=pval)
		{
			double v0 = ::Rf_qchisq(pval_noadj, 1, FALSE, FALSE);
			double v1 = ::Rf_qchisq(pval, 1, FALSE, FALSE);

// Rprintf("VarS/VarS_org (%d/%d): %g [MAC:%g, p0:%.3g, p1:%.3g]\n", i+1, ncol_g, v0/v1, mac, pval_noadj, pval);

			double r = sqrt(v0 / v1) * vr_sqrt;
			if (R_FINITE(r)) var_ratio[i] = r;
		}
	}
	// G_tilde' P G_tilde
	// G'SiG - G'SiXUG - (G'SiXUG)' + G'U'X'SiXUG, where U=XVX_inv_XV
	if (p_struct_skat->is_sparse)
	{
	#ifdef TIMING
		auto_timing tm3(run_time[n_st_skat+3]);
	#endif
		// GPG = G0.t() * p_struct_skat->Sigma_inv_s * G0;
		calc_mat_G_Si_G(G0, p_struct_skat->Sigma_inv_s, GPG);
	} else {
		GPG = G0.t() * p_struct_skat->Sigma_inv_d * G0;
	}
	{
	#ifdef TIMING
		auto_timing tm4(run_time[n_st_skat+4]);
	#endif
		dmat UG(p_struct_skat->XVX_inv_XV * G0);
		{
			// - G'SiXUG - (G'SiXUG)'
			dmat GSiX(p_struct_skat->Si_X * G0);
			dmat M = GSiX.t() * UG;
			GPG -= M;
			GPG -= M.t();
		}
		dmat GUXSiX(p_struct_skat->XVX_inv_XV_X_Si_X * G0);
		GPG += GUXSiX.t() * UG;
	}
}

/// calculation SKAT p-value with beta parameters
static void gmat_skat_test_p2(const sp_mat &G0, const dvec &Ts, const dmat &GPG,
	double beta_b1, double beta_b2,
	const double maf[], const double var_ratio[],
	double w_skat[], double out_ans[])
{
#ifdef TIMING
	auto_timing tm2(run_time[n_st_skat+2]);
#endif
	const int ncol_g = G0.n_cols;
	if (ncol_g <= 0)
	{
		out_ans[0] = R_NaN;
		return;
	}
	// set weights (w_skat)
	double w_sum = 0;
	for (int i=0; i < ncol_g; i++)
	{
		double v = w_skat[i] = Rf_dbeta(maf[i], beta_b1, beta_b2, FALSE);
		w_sum += v;
	}
	f64_mul(ncol_g, 1/w_sum, w_skat);  // normalize

	// Q stat
	double Q = 0;
	for (int i=0; i < ncol_g; i++) Q += sq(w_skat[i] * Ts[i]);
	Q /= sq(mod_tau[0]);  // adjust according to variance Y and GPG

	// update w_skat with variance ratio
	for (int i=0; i < ncol_g; i++) w_skat[i] *= var_ratio[i];

	// UG'PGU (U = w_skat)
	dmat D(size(GPG));
	{
		double *p = &D[0];
		const double *s = &GPG[0];
		for (int i=0; i < ncol_g; i++)
		{
			const double w = w_skat[i];
			for (int j=0; j < ncol_g; j++)
				p[j] = w * w_skat[j] * s[j];
			p += ncol_g; s += ncol_g;
		}
	}

	// check r_min (additional adjustment for case-control imbalance)
	if ((mod_trait == TTrait::Binary) && (threshold_pval_spa > 0))
	{
		// need a burden test
		dvec g_b;
		g_b.zeros(G0.n_rows);
		// collapse genotypes, using the original weights for SKAT
		for (int i=0; i < ncol_g; i++)
			add_g_w(g_b, G0, i, w_skat[i]/var_ratio[i]);
		// get Score and variance for g_b
		double S=0, pval=1;
		{
			// not need sparse GRM here
			const double *old = mod_sigma_inv_val;
			mod_sigma_inv_val = NULL;
			g_score_test(&g_b[0], -1, NULL, NULL, &pval, NULL, NULL, &S, false);
			mod_sigma_inv_val = old;
		}
// Rprintf("pval: %g, noabj: %g, S^2: %g\n", pval, pval_noabj, sq(S));  // debug

		if (R_FINITE(pval) && (pval > 0) && (pval < 1))
		{
			double V_sum = sq(S) / ::Rf_qchisq(pval, 1, FALSE, FALSE);
			double varQ = f64_sum(size_t(ncol_g)*ncol_g, &D[0]);
			double r_min = varQ / V_sum;
			if (r_min < 1) D *= 1/r_min;
// Rprintf("V_sum: %g, varQ: %g\n", V_sum, varQ);  // debug
// Rprintf("r_min: %g\n", r_min);  // debug
		}
	}

	// get eigenvalues
	dvec ev = (ncol_g > 1) ? eig_sym(D) : D;

	// p-value
	double pval;
	if (ncol_g == 1)
	{
		// regular chi-square distribution
		pval = ::Rf_pchisq(Q/ev[0], 1, FALSE, FALSE);
	} else {
		// find the first positive value
		int i_st=0;
		for (; i_st < ncol_g; i_st++) if (ev[i_st] > 0) break;
		// from largest to smallest, and remove any negative value
		NumericVector e(ncol_g - i_st);
		for (int i=0; i < ncol_g-i_st; i++)
			e[i] = ev[ncol_g-i-1];
		// call CompQuadForm::davies()
		pval = Rf_asReal(p_struct_skat->chisq_pval(Q, e));
	}
	out_ans[0] = pval;
}


/// calculate SKAT p-values
RcppExport SEXP saige_skat_test_pval(SEXP dosage)
{
BEGIN_RCPP
#ifdef TIMING
	auto_timing tm(run_time[n_st_skat+0]);
#endif
	// buffer
	double *maf = buf_unitsz;
	double *mac = buf_unitsz + num_unitsz;
	double *mac_imp = buf_unitsz + 2*num_unitsz;
	double *var_ratio = buf_unitsz + 3*num_unitsz;
	double *w_skat = buf_unitsz + 4*num_unitsz;

	// get genotype matrix
	sp_mat G0 = get_G0_flipped_impute(dosage, maf, mac, mac_imp);
	const int n_snv = G0.n_cols;

	// returned values
	const int st_idx = AGGR_INDEX_START + 3;
	NumericVector ans(st_idx + num_wbeta);
	summary_maf_mac(ans, n_snv, maf, mac);

	// collapse ultra rare variants
	// maf could be revised according to the new G0
	int n_collapse;
	G0 = Misc::GetSp_CollapseGenoMat(G0, threshold_skat_mac,
		p_struct_skat->collapse_method, mac_imp, maf, n_collapse);
	const int g_ncol = G0.n_cols;
	ans[9]  = n_collapse;
	ans[10] = g_ncol;
	ans[11] = min(sum(G0, 0));

	// Tstat for each G_i in G0
	dvec Ts;
	// G_tilde' P G_tilde
	// G'SiG - G'SiXUG - (G'SiXUG)' + G'U'X'SiXUG, where U=XVX_inv_XV
	dmat GPG;
	gmat_skat_test_p1(G0, var_ratio, Ts, GPG);

	// for each beta weight
	for (int w_i=0; w_i < num_wbeta; w_i++)
	{
		// get weights, beta function
		const double b1 = buf_wbeta[2*w_i+0], b2 = buf_wbeta[2*w_i+1];
		// calculation
		double pval;
		gmat_skat_test_p2(G0, Ts, GPG, b1, b2, maf, var_ratio, w_skat, &pval);
		// output
		ans[st_idx + w_i] = pval;
	}

	// output
	return ans;
END_RCPP
}


// ====================================

static double acat_pval(R_xlen_t n, const double pval[], const double w[],
	bool throw_error);

static void gmat_acatv_test(const sp_mat &G0, double beta_b1, double beta_b2,
	const double maf[], const double mac[], const double mac_imp[],
	double pvals[], double w_pval[], double out_ans[])
{
	const int n_snv = G0.n_cols;
	if (n_snv <= 0)
	{
		out_ans[0] = out_ans[1] = 0;
		out_ans[2] = out_ans[3] = out_ans[4] = out_ans[5] = R_NaN;
		return;
	}
	// initialize
	int n_test   = 0;  // # of tests including single variant tests
	int n_burden = 0;  // # of ultra rare SNVs for burden test

	// for-loop for each variant
	for (int i=0; i < n_snv; i++)
	{
		const double C = mac[i];
		if (R_FINITE(C) && (C > 0))
		{
			if (C > threshold_acatv_mac)
			{
				if (!R_FINITE(pvals[i]))
				{
					dvec G(G0.col(i));
					// p-value calculation
					double pval=R_NaN, pval_noadj=pval;
					g_score_test(&G[0], mac_imp[i], NULL, NULL,
						&pval, &pval_noadj, NULL, NULL, SPA_always_use_fastSPA);
					pvals[i] = pval;  // p-value for this SNV
				}
				// save
				const double p = maf[i];
				w_pval[i] = sq(Rf_dbeta(p, beta_b1, beta_b2, FALSE)) * p * (1-p);
				n_test ++;
			} else {
				// burden test for ultra rare variants
				n_burden ++;
			}
		}
	}

	// if collapsed SNVs for burden test
	int ultra_var_idx = -1;
	if (n_burden > 0)
	{
		dvec G;
		G.zeros(mod_NSamp);
		double sum_w=0, summaf=0, summac=0;
		// get G and sum
		for (int i=0; i < n_snv; i++)
		{
			const double C = mac[i];
			if (R_FINITE(C) && (C > 0) && (C <= threshold_acatv_mac))
			{
				double w = Rf_dbeta(maf[i], beta_b1, beta_b2, FALSE);
				add_g_w(G, G0, i, w);
				sum_w += w;         // add weight
				summaf += maf[i]; summac += mac[i];
				ultra_var_idx = i;  // save the position
			}
		}
		// normalize G
		G *= 1 / sum_w;
		// burden test
		double pval=R_NaN, pval_noadj=pval;
		if ((summac > 0) && (summac >= threshold_summac))
		{
			// p-value calculation
			g_score_test(&G[0], summac, NULL, NULL, &pval, &pval_noadj,
				NULL, NULL, false);
		}
		if (R_FINITE(pval))
		{
			const double p = summaf / n_burden;
			w_pval[ultra_var_idx] =
				sq(Rf_dbeta(p, beta_b1, beta_b2, FALSE)) * p * (1-p);
			pvals[ultra_var_idx] = pval;
			n_test ++;
		}
	}

	// set the output
	out_ans[0] = n_test - (n_burden > 0);  // n.single
	out_ans[1] = n_burden;  // n.burden
	out_ans[2] =
		(n_test > 0) ? acat_pval(n_snv, pvals, w_pval, false) : R_NaN;
	f64_medmaxmin(pvals, n_snv, out_ans[3], out_ans[4], out_ans[5]);

	// clear
	if (ultra_var_idx >= 0) pvals[ultra_var_idx] = R_NaN;
}
	

/// calculate ACAT-V p-values
RcppExport SEXP saige_acatv_test_pval(SEXP dosage)
{
BEGIN_RCPP
	// buffer
	double *maf = buf_unitsz;
	double *mac = buf_unitsz + num_unitsz;
	double *mac_imp = buf_unitsz + 2*num_unitsz;
	double *w_pval = buf_unitsz + 3*num_unitsz;
	double *pvals  = buf_unitsz + 4*num_unitsz;

	// get genotype matrix
	sp_mat G0 = get_G0_flipped_impute(dosage, maf, mac, mac_imp);
	const int n_snv = G0.n_cols;
	// initialize pvals with NaN
	for (int i=0; i < n_snv; i++) pvals[i] = R_NaN;

	// summarize maf & mac
	const int st_idx = AGGR_INDEX_START + 2;
	NumericVector ans(st_idx + num_wbeta*4);
	summary_maf_mac(ans, n_snv, maf, mac);

	// for each beta weight
	for (int w_i=0; w_i < num_wbeta; w_i++)
	{
		// get weights, beta function
		const double b1 = buf_wbeta[2*w_i+0], b2 = buf_wbeta[2*w_i+1];
		// calculate individual p-value, if pvals[j]==NaN
		double v[6];
		gmat_acatv_test(G0, b1, b2, maf, mac, mac_imp, pvals, w_pval, v);
		// set the output
		if (w_i == 0)
		{
			ans[9]  = v[0];  // n.single
			ans[10] = v[1];  // n.burden
		}
		const int k = st_idx + w_i * 4;
		ans[k+0] = v[2];  // pval
		ans[k+1] = v[3];  // p.median
		ans[k+2] = v[4];  // p.min
		ans[k+3] = v[5];  // p.max
	}

	// output
	return ans;
END_RCPP
}



// ====================================

/// calculate ACAT-O p-values
RcppExport SEXP saige_acato_test_pval(SEXP dosage)
{
BEGIN_RCPP
	// buffer
	double *maf = buf_unitsz;
	double *mac = buf_unitsz + num_unitsz;
	double *mac_imp = buf_unitsz + 2*num_unitsz;
	double *ws  = buf_unitsz + 3*num_unitsz;
	double *pvals = buf_unitsz + 4*num_unitsz;
	double *var_ratio = buf_unitsz + 5*num_unitsz;
	double *maf_s = buf_unitsz + 6*num_unitsz;

	// get genotype matrix
	sp_mat G0 = get_G0_flipped_impute(dosage, maf, mac, mac_imp);
	const int n_snv = G0.n_cols;
	// initialize pvals with NaN for ACAT-V
	for (int i=0; i < n_snv; i++) pvals[i] = R_NaN;

	// summarize maf & mac
	const int st_idx = AGGR_INDEX_START + 1;
	const int n_each_wb = p_struct_skat ? 3 : 2;
	NumericVector ans(st_idx + n_each_wb*num_wbeta);
	summary_maf_mac(ans, n_snv, maf, mac);

	// initialize SKAT structure
	sp_mat G0_s;  // genotype matrix for SKAT
	dvec Ts;      // Tstat for each G_i in G0
	dmat GPG;     // G_tilde' P G_tilde
	if (p_struct_skat)
	{
		memcpy(maf_s, maf, sizeof(double)*n_snv);
		int n_collapse = 0;
		G0_s = Misc::GetSp_CollapseGenoMat(G0, threshold_skat_mac,
			p_struct_skat->collapse_method,
			mac_imp, maf_s, n_collapse);
		// maf_s could be revised according to collapsed genotypes
		gmat_skat_test_p1(G0_s, var_ratio, Ts, GPG);
	}

	// for each beta weight
	for (int i=0; i < num_wbeta; i++)
	{
		// get weights, beta function
		const double b1 = buf_wbeta[2*i+0], b2 = buf_wbeta[2*i+1];
		// calculation
		const int k = st_idx + n_each_wb*i;
		double v[6];
		// burden p-value
		gmat_burden_test(G0, b1, b2, maf, mac, ws, v);
		ans[k+0] = v[3];
		// calculate individual ACAT-V p-value, if pvals[j]==NaN
		gmat_acatv_test(G0, b1, b2, maf, mac, mac_imp, pvals, ws, v);
		ans[k+1] = v[2];
		// SKAT p-value
		if (p_struct_skat)
		{
			double pval;
			gmat_skat_test_p2(G0_s, Ts, GPG, b1, b2, maf_s, var_ratio, ws,
				&pval);
			ans[k+2] = pval;
		}
	}

	// combined p-value
	const int n_pval = n_each_wb*num_wbeta;
	ws = buf_unitsz;
	if (n_pval > 5*num_unitsz)
		ws = REAL(NEW_NUMERIC(n_pval));
	for (int i=0; i < n_pval; i++) ws[i] = 1;
	ans[st_idx-1] = acat_pval(n_pval, &ans[st_idx], ws, false);

	// output
	return ans;
END_RCPP
}



// ========================================================================= //

#ifdef HAVE_ATANPI
extern double atanpi(double x);
#else
inline static double atanpi(double x) { return atan(x) / M_PI; }
#endif

static const double ROUND_ZERO = 1e-300;
static const double ROUND_ONE  = 1 - 1e-16;

/// p-value from ACAT combination method
static double acat_pval(R_xlen_t n, const double pval[], const double w[],
	bool throw_error)
{
	// get the weight sum
	double sumw = 0;
	for (R_xlen_t i=0; i < n; i++)
		if (R_FINITE(pval[i]) && R_FINITE(w[i])) sumw += w[i];
	if (sumw <= 0)
	{
		if (throw_error)
			Rf_error("the sum of weights should be > 0.");
		else
			return R_NaN;
	}
	// get statistic
	double Tstat = 0;
	for (R_xlen_t i=0; i < n; i++)
	{
		double p = pval[i];
		if (R_FINITE(p) && R_FINITE(w[i]))
		{
			// check p-value
			if (p < 0 || p > 1)
			{
				if (throw_error)
					Rf_error("Invalid input p-value: %g.", p);
				else
					return R_NaN;
			}
			if (p < ROUND_ZERO)
				p = ROUND_ZERO;  // almost the smallest number > 0
			else if (p > ROUND_ONE)
				p = ROUND_ONE;  // almost the closest number around 1
			// calc stat
			if (p >= 1e-15)
			{
				Tstat += w[i] * tanpi(0.5 - p);
			} else {
				// based on the taylor series expansion,
				//   Series[tan[(1/2-x)*pi], {x, 0, 5}]
				//   1/(pi*x) - pi/3*x - pi^3/45*x^3 - 2*pi^5/945*x^5 + O(x^7)
				Tstat += w[i] / p / M_PI;
			}
		}
	}
	Tstat /= sumw;
	// get p-value from Tstat, and return
	if (Tstat <= 5e+14)
		return 0.5 - atanpi(Tstat);
	else
		return 1.0 / Tstat * M_1_PI;
}

RcppExport SEXP saige_acat_p(SEXP pval, SEXP weight)
{
	const R_xlen_t n = Rf_xlength(pval);
	// check pval
	if (n <= 0)
		Rf_error("the number of p-values should be > 0.");
	else if (n == 1)
		return pval;
	// check weight
	if (Rf_isNull(weight))
	{
		weight = NEW_NUMERIC(n);
		double *w = REAL(weight);
		for (R_xlen_t i=0; i < n; i++) w[i] = 1;
	}
	// check
	if (n != Rf_xlength(weight))
		Rf_error("weights should have the same length as p-values.");
	if (TYPEOF(pval) != REALSXP)
		Rf_error("p-values should be numeric.");
	if (TYPEOF(weight) != REALSXP)
		Rf_error("weights should be numeric.");
	// calculate
	double v = acat_pval(n, REAL(pval), REAL(weight), true);
	// output
	return Rf_ScalarReal(v);
}


// ========================================================================= //

inline static const char *b2s(bool v) { return v ? "true" : "false"; }

RcppExport SEXP saige_set_option(SEXP val, SEXP use_avx, SEXP Rverbose)
{
	const int avx = Rf_asInteger(use_avx);
	const bool verbose = (Rf_asLogical(Rverbose)==TRUE);

	const bool old_use_fastSPA = SPA_always_use_fastSPA;
	SPA_always_use_fastSPA = (Rf_asLogical(val) == TRUE);

	extern bool fc_use_avx512f;
	extern bool fc_use_avx2;
	const bool old_fc_use_avx512f = fc_use_avx512f;
	const bool old_fc_use_avx2 = fc_use_avx2;
	switch (avx)
	{
		case 1:
			fc_use_avx512f = true; fc_use_avx2 = true; break;
		case 2:
			fc_use_avx512f = false; fc_use_avx2 = true; break;
		case 3:
			fc_use_avx512f = fc_use_avx2 = false; break;
	}
	vec_init_function();

	if (verbose)
	{
		Rprintf("SPA_always_use_fastSPA: %s => %s\n",
			b2s(old_use_fastSPA), b2s(SPA_always_use_fastSPA));
		Rprintf("fc_use_avx512f: %s => %s\n",
			b2s(old_fc_use_avx512f), b2s(fc_use_avx512f));
		Rprintf("fc_use_avx2: %s => %s\n",
			b2s(old_fc_use_avx2), b2s(fc_use_avx2));
	}
	return R_NilValue;
}

RcppExport SEXP saige_simd_version();
RcppExport SEXP saige_store_2b_geno(SEXP, SEXP, SEXP, SEXP, SEXP);
RcppExport SEXP saige_store_sp_geno(SEXP, SEXP, SEXP, SEXP, SEXP);

/// initialize the package
RcppExport void R_init_SAIGEgds(DllInfo *info)
{
	#define CALL(name, num)	   { #name, (DL_FUNC)&name, num }

	static R_CallMethodDef callMethods[] =
	{
		CALL(saige_score_test_init, 1),
		CALL(saige_simd_version, 0),
		CALL(saige_store_2b_geno, 5),
		CALL(saige_store_sp_geno, 5),
		CALL(saige_acat_p, 2),
		{ NULL, NULL, 0 }
	};

	R_registerRoutines(info, NULL, callMethods, NULL, NULL);
	vec_init_function();
}
