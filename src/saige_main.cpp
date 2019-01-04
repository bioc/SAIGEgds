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


#include <RcppArmadillo.h>
#include <algorithm>

using namespace Rcpp;
using namespace arma;


/// sum_i x[i]*y[i]
extern "C" double vec_dot(size_t n, const double *x, const double *y);
/// sum_i x[i]*y[i]*y[i]
extern "C" double vec_dot_sp(size_t n, const double *x, const double *y);



/// Get the list element named str, or return NULL
static SEXP GetListElement(SEXP list, const char *str)
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
static double *model_y = NULL;
static double *model_mu = NULL;
static double *model_y_mu = NULL;
static double *model_mu2 = NULL;
static SEXP model_null_XXVX_inv = NULL;
static SEXP model_null_XV = NULL;
static double model_varRatio = 0;



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
	threshold_maf = Rf_asReal(GetListElement(model, "maf"));
	if (!R_FINITE(threshold_maf)) threshold_maf = -1;
	threshold_mac = Rf_asReal(GetListElement(model, "mac"));
	if (!R_FINITE(threshold_mac)) threshold_mac = -1;
	// model parameters
	SEXP y = GetListElement(model, "y");
	model_num_samp = Rf_length(y);
	model_y = REAL(y);
	model_mu = REAL(GetListElement(model, "mu"));
	model_y_mu = REAL(GetListElement(model, "y_mu"));
	model_mu2 = REAL(GetListElement(model, "mu2"));
	model_null_XXVX_inv = GetListElement(model, "XXVX_inv");
	model_null_XV = GetListElement(model, "XV");
	model_varRatio = Rf_asReal(GetListElement(model, "var.ratio"));
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
		// genotype vector, reuse memory and avoid extra copy
		colvec G(ds.begin(), num_samp, false);
		// XV matrix, reuse memory
		NumericMatrix mt1(model_null_XV);
		mat XV(mt1.begin(), mt1.nrow(), mt1.ncol(), false);
		// XXVX_inv matrix, reuse memory
		NumericMatrix mt2(model_null_XXVX_inv);
		mat XXVX_inv(mt2.begin(), mt2.nrow(), mt2.ncol(), false);
		// adjusted genotypes
		colvec g = G - XXVX_inv * (XV * G);

		// inner product
		double S = vec_dot(num_samp, model_y_mu, &g[0]);
		double var = vec_dot_sp(num_samp, model_mu2, &g[0]) * model_varRatio;

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
