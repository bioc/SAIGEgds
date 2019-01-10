// ===========================================================
//
// saige_main.cpp: SAIGE association analysis for each variant
//
// Copyright (C) 2019    Xiuwen Zheng
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
// with SNPRelate.
// If not, see <http://www.gnu.org/licenses/>.


#include <Rcpp.h>
#include <algorithm>

using namespace Rcpp;


/// y[i] = x - y[i]
extern "C" void f64_sub(size_t n, double x, double *y);
/// y[i] = x * y[i]
extern "C" void f64_mul(size_t n, double x, double *y);
/// sum_i x[i]*y[i]
extern "C" double f64_dot(size_t n, const double *x, const double *y);
/// out1 = sum_i x1[i]*y[i], out2 = sum_i x2[i]*y[i]*y[i]
extern "C" void f64_dot_sp(size_t n, const double *x1, const double *x2,
	const double *y, double &out1, double &out2);
/// vec(p_m) = mat(x_{m*n}) * vec(y_n), y is a sparse vector
extern "C" void f64_mul_mat_vec(size_t n, size_t m, const double *x, const double *y, double *p);
/// vec(p_n) = vec(x_n) - t(mat(y_{m*n})) * vec(z_m)
extern "C" void f64_sub_mul_mat_vec(size_t n, size_t m,
	const double *x, const double *y, const double *z, double *p);

/// SPAtest
extern "C" double Saddle_Prob(double q, double m1, double var1, size_t n_g,
	const double mu[], const double g[], double cutoff, bool &converged);


/// Get the list element named str, or return NULL
static SEXP get_item(SEXP list, const char *str)
{
	SEXP names = Rf_getAttrib(list, R_NamesSymbol);
	R_xlen_t n = XLENGTH(list);
	for (R_xlen_t i=0; i < n; i++)
	{
		if (strcmp(CHAR(STRING_ELT(names, i)), str) == 0)
			return VECTOR_ELT(list, i);
	}
	return R_NilValue;
}


static double threshold_maf = 0;  //< the threshold of MAF filter
static double threshold_mac = 0;  //< the threshold of MAC filter

static int model_NSamp = 0;   //< the number of samples
static int model_NCoeff = 0;  //< the number of beta coefficients

static double *model_y = NULL;
static double *model_mu = NULL;
static double *model_y_mu = NULL;
static double *model_mu2 = NULL;
static double *model_t_XXVX_inv = NULL;
static double *model_XV = NULL;

static double model_varRatio = 0;

static double *buf_coeff = NULL;
static double *buf_adj_g = NULL;


/// return allele frequency and impute genotype using the mean
static void SummaryAndImputeGeno(double *ds, size_t n, double &AF, double &AC,
	int &Num)
{
	double sum = 0; int num = 0;
	for (size_t i=0; i < n; i++)
		if (R_FINITE(ds[i])) { sum += ds[i]; num ++; }
	AF = (num > 0) ? (sum/(2*num)) : R_NaN;
	AC = sum; Num = num;
	if (num < (int)n)
	{
		double d = AF * 2;
		for (size_t i=0; i < n; i++)
			if (!R_FINITE(ds[i])) ds[i] = d;
	}
}



/// initialize model objects
RcppExport SEXP saige_score_test_init(SEXP model)
{
BEGIN_RCPP
	// threshold setting
	threshold_maf = Rf_asReal(get_item(model, "maf"));
	if (!R_FINITE(threshold_maf)) threshold_maf = -1;
	threshold_mac = Rf_asReal(get_item(model, "mac"));
	if (!R_FINITE(threshold_mac)) threshold_mac = -1;
	// model parameters
	model_NSamp = Rf_length(get_item(model, "y"));
	model_y = REAL(get_item(model, "y"));
	model_mu = REAL(get_item(model, "mu"));
	model_y_mu = REAL(get_item(model, "y_mu"));
	model_mu2 = REAL(get_item(model, "mu2"));
	model_t_XXVX_inv = REAL(get_item(model, "t_XXVX_inv"));
	model_XV = REAL(get_item(model, "XV"));
	model_varRatio = Rf_asReal(get_item(model, "var.ratio"));
	NumericMatrix m(get_item(model, "XV"));
	model_NCoeff = m.nrow();
	// buffer
	buf_coeff = REAL(get_item(model, "buf1"));
	buf_adj_g = REAL(get_item(model, "buf2"));
	return R_NilValue;
END_RCPP
}


/// 
RcppExport SEXP saige_score_test_quant(SEXP Dosage)
{
BEGIN_RCPP
	return R_NilValue;
END_RCPP
}



/// return AF, AC, num, beta, se, pval
RcppExport SEXP saige_score_test_bin(SEXP Dosage)
{
BEGIN_RCPP

	// dosages
	NumericVector ds(Dosage);
	const size_t num_samp = ds.size();
	// calc allele freq, and impute geno using the mean
	double AF, AC;
	int Num;
	SummaryAndImputeGeno(&ds[0], num_samp, AF, AC, Num);

	double maf = std::min(AF, 1-AF);
	double mac = std::min(AC, 2*Num - AC);
	if (Num>0 && maf>=threshold_maf && mac>=threshold_mac)
	{
		bool minus = (AF > 0.5);
		if (minus) f64_sub(model_NSamp, 2, &ds[0]);

		// adj_g = G - XXVX_inv * (XV * G), adjusted genotypes
		// coeff = XV * G
		f64_mul_mat_vec(model_NSamp, model_NCoeff, model_XV, &ds[0], buf_coeff);
		// adj_g = G - XXVX_inv * coeff
		f64_sub_mul_mat_vec(model_NSamp, model_NCoeff,
			&ds[0], model_t_XXVX_inv, buf_coeff, buf_adj_g);

		// inner product
		double S, var;
		// S = sum((y - mu) .* adj_g)
		// var = sum(mu*(1-mu) .* adj_g .* adj_g)
		f64_dot_sp(model_NSamp, model_y_mu, model_mu2, buf_adj_g, S, var);
		var *= model_varRatio;

		// p-value
		double pval_noadj = ::Rf_pchisq(S*S/var, 1, FALSE, FALSE);
		double pval = pval_noadj;
		double beta = (minus ? -1 : 1) * S / var;
		bool converged = true;

		// need SPAtest or not?
		if (R_FINITE(pval_noadj) && pval_noadj <= 0.05)
		{
			double AC2 = minus ? (2*Num - AC) : AC;
			// adj_g = adj_g / sqrt(AC2)
			f64_mul(model_NSamp, 1/sqrt(AC2), buf_adj_g);
			// q = sum(y .* adj_g)
			double q = f64_dot(model_NSamp, model_y, buf_adj_g);
			double m1, var2;
			// m1 = sum(mu .* adj_g)
			// var2 = sum(mu*(1-mu) .* adj_g .* adj_g)
			f64_dot_sp(model_NSamp, model_mu, model_mu2, buf_adj_g, m1, var2);
			double var1 = var2 * model_varRatio;
			double Tstat = q - m1;
			double qtilde = Tstat/sqrt(var1) * sqrt(var2) + m1;
			// call Saddle_Prob in SPAtest
			pval = Saddle_Prob(qtilde, m1, var2, model_NSamp, model_mu,
				buf_adj_g, 2, converged);
			beta = (Tstat / var1) / sqrt(AC2);
		}
		double SE = abs(beta/::Rf_qnorm5(pval/2, 0, 1, TRUE, FALSE));

		NumericVector ans(8);
		ans[0] = AF;    ans[1] = AC;    ans[2] = Num;
		ans[3] = beta;  ans[4] = SE;    ans[5] = pval;
		ans[6] = pval_noadj;
		ans[7] = converged ? 1 : 0;
		return ans;
	} else {
		return R_NilValue;
	}

END_RCPP
}


/// initialize the package
RcppExport void R_init_SAIGEgds(DllInfo *info)
{
	#define CALL(name, num)	   { #name, (DL_FUNC)&name, num }

	static R_CallMethodDef callMethods[] =
	{
		CALL(saige_score_test_init, 1),
		{ NULL, NULL, 0 }
	};

	R_registerRoutines(info, NULL, callMethods, NULL, NULL);
}
