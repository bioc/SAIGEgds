// ===========================================================
//
// saige_main.cpp: Identity by state (IBS) Analysis on GWAS
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


/// sum_i x[i]*y[i]
extern "C" double f64_dot(size_t n, const double *x, const double *y);
/// sum_i x[i]*y[i]*y[i]
extern "C" double f64_dot_sp(size_t n, const double *x, const double *y);
/// vec(p_m) = mat(x_{m*n}) * vec(y_n), y is a sparse vector
extern "C" void f64_mul_mat_vec(size_t n, size_t m, const double *x, const double *y, double *p);
/// vec(p_n) = vec(x_n) - t(mat(y_{m*n})) * vec(z_m)
extern "C" void f64_sub_mul_mat_vec(size_t n, size_t m,
	const double *x, const double *y, const double *z, double *p);



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


static double threshold_maf = 0;
static double threshold_mac = 0;

static int model_num_samp = 0;
static int model_num_coeff = 0;
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
	double sum = 0;
	int num = 0;
	for (size_t i=0; i < n; i++)
	{
		if (R_FINITE(ds[i]))
		{
			sum += ds[i];
			num ++;
		}
	}
	AC = sum;
	AF = (num > 0) ? (sum/(2*num)) : R_NaN;
	Num = num;
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
	model_num_samp = Rf_length(get_item(model, "y"));
	model_y = REAL(get_item(model, "y"));
	model_mu = REAL(get_item(model, "mu"));
	model_y_mu = REAL(get_item(model, "y_mu"));
	model_mu2 = REAL(get_item(model, "mu2"));
	model_t_XXVX_inv = REAL(get_item(model, "t_XXVX_inv"));
	model_XV = REAL(get_item(model, "XV"));
	model_varRatio = Rf_asReal(get_item(model, "var.ratio"));
	NumericMatrix m(get_item(model, "XV"));
	model_num_coeff = m.nrow();
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
		// adj_g = G - XXVX_inv * (XV * G), adjusted genotypes
		// coeff = XV * G
		f64_mul_mat_vec(model_num_samp, model_num_coeff, model_XV, &ds[0], buf_coeff);
		// adj_g = G - XXVX_inv * coeff
		f64_sub_mul_mat_vec(model_num_samp, model_num_coeff,
			&ds[0], model_t_XXVX_inv, buf_coeff, buf_adj_g);

		// inner product
		// S = sum(model_y_mu .* adj_g)
		double S = f64_dot(model_num_samp, model_y_mu, buf_adj_g);
		// var = sum(model_mu2 .* adj_g .* adj_g) * varRatio
		double var = f64_dot_sp(model_num_samp, model_mu2, buf_adj_g) *
			model_varRatio;

		// p-value
		double pval_noadj = ::Rf_pchisq(S*S/var, 1, FALSE, FALSE);
		double pval = pval_noadj;
		double beta = S / var;
		double se   = abs(beta/::Rf_qnorm5(pval_noadj/2, 0, 1, TRUE, FALSE));

		// need SPAtest or not?
		if (pval_noadj <= 0.05)
		{
		}

		NumericVector ans(7);
		ans[0] = AF;    ans[1] = AC;    ans[2] = Num;
		ans[3] = beta;  ans[4] = se;    ans[5] = pval;
		ans[6] = pval_noadj;
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
