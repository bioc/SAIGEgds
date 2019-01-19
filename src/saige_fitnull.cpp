#include <RcppArmadillo.h>
#include "vectorization.h"

#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

using namespace Rcpp;
using namespace arma;
using namespace std;


// ========================================================================= //
// Store packed SNP genotypes

static unsigned char *Geno_PackedRaw = NULL;
static size_t Geno_NumSamp = 0;
static size_t Geno_PackedNumSamp = 0;
static size_t Geno_NumVariant = 0;

static double *buf_std_geno = NULL;  //< a 4-by-n_variant look-up matrix
static double *buf_diag_grm = NULL;  //< n_samp-length, sigma_i = sum_j adj.g[i,j]^2


// Internal 2-bit genotype lookup tables
static bool lookup_has_init = false;
static unsigned char num_valid[256], num_sum[256];

static void init_lookup_table()
{
	if (!lookup_has_init)
	{
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
}


inline static double sq(double v) { return v*v; }

/// store SNP genotype in the 2-bit packed format
RcppExport SEXP saige_store_geno(SEXP rawgeno, SEXP num_samp, SEXP buf_geno,
	SEXP buf_sigma)
{
BEGIN_RCPP

	// initialize the basic genotype variables
	RawMatrix RawGeno(rawgeno);
	Geno_PackedRaw = (unsigned char*)&RawGeno[0];
	Geno_NumSamp = Rf_asInteger(num_samp);
	Geno_PackedNumSamp = RawGeno.nrow();
	Geno_NumVariant = RawGeno.ncol();

	// build the look-up table of standardized genotypes
	init_lookup_table();
	buf_std_geno = REAL(buf_geno);
	for (size_t i=0; i < Geno_NumVariant; i++)
	{
		unsigned char *g = Geno_PackedRaw + Geno_PackedNumSamp*i;
		// calculate allele frequency
		int n_valid=0, sum=0;
		for (size_t j=0; j < Geno_PackedNumSamp; j++)
		{
			n_valid += num_valid[g[j]];
			sum += num_sum[g[j]];
		}
		double af = double(sum) / (2*n_valid);
		double inv = 1 / sqrt(2*af*(1-af));
		if (!R_FINITE(af) || !R_FINITE(inv))
			af = inv = 0;
		double *p = &buf_std_geno[4*i];
		p[0] = (0 - 2*af) * inv; p[1] = (1 - 2*af) * inv;
		p[2] = (2 - 2*af) * inv; p[3] = 0;
	}

	// calculate diag(grm)
	buf_diag_grm = REAL(buf_sigma);
	memset(buf_sigma, 0, sizeof(double)*Geno_NumSamp);
	for (size_t i=0; i < Geno_NumVariant; i++)
	{
		unsigned char *g = Geno_PackedRaw + Geno_PackedNumSamp*i;
		const double *base = buf_std_geno + 4*i;
		size_t n = Geno_NumSamp;
		double *p = buf_diag_grm;
		for (; n >= 4; n-=4, p+=4)
		{
			unsigned char gg = *g++;
			p[0] += sq(base[gg & 0x03]);
			p[1] += sq(base[(gg >> 2) & 0x03]);
			p[2] += sq(base[(gg >> 4) & 0x03]);
			p[3] += sq(base[gg >> 6]);
		}
		for (unsigned char gg = (n>0 ? *g : 0); n > 0; n--)
		{
			(*p++) += sq(base[gg & 0x03]);
			gg >>= 2;
		}
	}
	f64_mul(Geno_NumSamp, 1.0 / Geno_NumVariant, buf_diag_grm);

END_RCPP
}



// ========================================================================= //

/// Cross-product of standardized genotypes and a numeric vector
/// Input: b (n_samp-length)
/// Output: out_b (n_samp-length)
static void get_crossprod_b_grm(const dcolvec &b, dvec &out_b)
{
	out_b.resize(Geno_NumSamp);
	memset(&out_b[0], 0, sizeof(double)*Geno_NumSamp);

	for (size_t i=0; i < Geno_NumVariant; i++)
	{
		unsigned char *g = Geno_PackedRaw + Geno_PackedNumSamp*i;
		const double *base = buf_std_geno + 4*i;

		// get dot = sum(std.geno .* b)
		double dot = 0;
		const double *pb = &b[0];
		size_t n = Geno_NumSamp;
		for (; n >= 4; n-=4, pb+=4)
		{
			unsigned char gg = *g++;
			dot += base[gg & 0x03] * pb[0] + base[(gg >> 2) & 0x03] * pb[1] +
				base[(gg >> 4) & 0x03] * pb[2] + base[gg >> 6] * pb[3];
		}
		for (unsigned char gg = (n>0 ? *g : 0); n > 0; n--)
		{
			dot += base[gg & 0x03] * (*pb++);
			gg >>= 2;
		}

		// update the output rv += dot .* std.geno
		double *pbb = &out_b[0];
		g = Geno_PackedRaw + Geno_PackedNumSamp*i;
		n = Geno_NumSamp;
		for (; n >= 4; n-=4, pbb+=4)
		{
			unsigned char gg = *g++;
			pbb[0] += dot * base[gg & 0x03];
			pbb[1] += dot * base[(gg >> 2) & 0x03];
			pbb[2] += dot * base[(gg >> 4) & 0x03];
			pbb[3] += dot * base[gg >> 6];
		}
		for (unsigned char gg = (n>0 ? *g : 0); n > 0; n--)
		{
			(*pbb++) += dot * base[gg & 0x03];
			gg >>= 2;
		}
	}

	// normalize
	f64_mul(Geno_NumSamp, 1.0/Geno_NumVariant, &out_b[0]);
}

  
/// Sigma = tau[0] * diag(1/W) + tau[1] * diag(grm)
/// Input: w, tau
/// Output: out_sigma
static void get_diag_sigma(const dvec& w, const dvec& tau, dvec &out_sigma)
{
	const double tau0 = tau[0];
	const double tau1 = tau[1];
	out_sigma.resize(Geno_NumSamp);
	for (size_t i=0; i < Geno_NumSamp; i++)
	{
		double v = tau0 / w[i] + tau1 * buf_diag_grm[i];
		if (v < 1e-4) v = 1e-4;
		out_sigma[i] = v;
	}
}


/// Sigma = tau[0] * b * diag(1/W) + tau[1] * diag(grm, b)
/// Input: w, tau
/// Output: out_sigma
static dcolvec get_crossprod(const dcolvec &b, const dvec& w, const dvec& tau)
{
	const double tau0 = tau[0];
	const double tau1 = tau[1];
Rprintf("tau0: %g, tau1: %g\n", tau0, tau1);
	if (tau1 == 0)
	{
		dcolvec rv = tau0 * (b % (1/w));
		return(rv);
	} else {
		dvec out_b;
		get_crossprod_b_grm(b, out_b);
		dcolvec rv = tau0 * (b % (1/w)) + tau1 * out_b;
		return(rv);
	}
}


/// Sigma = tau[0] * diag(1/W) + tau[1] * grm
/// Input: wVec, tauVec, bVec, maxiterPCG, tolPCG
static dvec get_PCG_diag_sigma(const dvec &w, const dvec &tau, const dvec &b,
	int maxiterPCG, double tolPCG)
{
	dvec r = b, r1;

	dvec crossProdVec(Geno_NumSamp);
	dvec minvVec;
	get_diag_sigma(w, tau, minvVec);
	minvVec = 1 / minvVec;

Rprintf("sum of minvVec: %g\n", sum(minvVec));

	double sumr2 = sum(r % r);

	dvec zVec = minvVec % r;
	dvec z1Vec;
	dvec pVec = zVec;

	dvec xVec(Geno_NumSamp);
	xVec.zeros();

	int iter = 0;
	while (sumr2 > tolPCG && iter < maxiterPCG)
	{
		iter = iter + 1;
		dcolvec ApVec = get_crossprod(pVec, w, tau);
		dvec preA = (r.t() * zVec)/(pVec.t() * ApVec);

		double a = preA(0);
		xVec = xVec + a * pVec;
		r1 = r - a * ApVec;

		z1Vec = minvVec % r1;
		dvec Prebet = (z1Vec.t() * r1)/(zVec.t() * r);
		double bet = Prebet(0);
		pVec = z1Vec+ bet*pVec;
		zVec = z1Vec;
		r = r1;

		sumr2 = sum(r % r);
	}

	if (iter >= maxiterPCG)
	{
		cout << "pcg did not converge. You may increase maxiter number." << endl;
	}
	// cout << "iter from get_PCG_diag_sigma " << iter << endl;
	return(xVec);
}


// [[export]]
NumericVector nb(int n) {
  	return(rbinom(n,1,0.5));
}

//This function calculates the coefficients of variation for mean of a vector
// [[export]]
float calCV(dvec& xVec){
  int veclen = xVec.n_elem;
  float vecMean = mean(xVec);
  float vecSd = stddev(xVec);
  float vecCV = (vecSd/vecMean)/veclen;
  return(vecCV);
}



// [[export]]
static double get_trace(const dmat &Sigma_iX, const dmat& Xmat,
	const dvec& wVec, const dvec& tauVec, const dmat& cov1,
	int nrun, int maxiterPCG, double tolPCG, double traceCVcutoff)
{
	// set_seed(200);
	dmat Sigma_iXt = Sigma_iX.t();
	dvec Sigma_iu;  
	dcolvec Pu;
	dvec Au;
	dvec uVec;

	int nrunStart = 0;
	int nrunEnd = nrun;
	double traceCV = traceCVcutoff + 0.1;
	dvec tempVec(nrun);
	tempVec.zeros();

	while(traceCV > traceCVcutoff)
	{
		//dvec tempVec(nrun);
		//tempVec.zeros();
		for(int i = nrunStart; i < nrunEnd; i++)
		{
			NumericVector uVec0;
			uVec0 = nb(Geno_NumSamp);
			uVec = as<dvec>(uVec0);
			uVec = uVec*2 - 1;
			Sigma_iu = get_PCG_diag_sigma(wVec, tauVec, uVec, maxiterPCG, tolPCG);
			Pu = Sigma_iu - Sigma_iX * (cov1 *  (Sigma_iXt * uVec));
			get_crossprod_b_grm(uVec, Au);
			tempVec(i) = dot(Au, Pu);
			Au.clear();
			Pu.clear();
			Sigma_iu.clear();
			uVec.clear();
		}
		traceCV = calCV(tempVec);
		if(traceCV > traceCVcutoff)
		{
			nrunStart = nrunEnd;
			nrunEnd = nrunEnd + 10;
			tempVec.resize(nrunEnd);
			cout << "CV for trace random estimator using "<< nrun << " runs is " << traceCV <<  " > " << traceCVcutoff << endl;
			cout << "try " << nrunEnd << " runs" << endl;
		}
	}

	double tra = mean(tempVec);
	tempVec.clear();
	return(tra);
}


/// Calculate fixed and random effect coefficients
/// Input:  Y, X, w, tau, maxiterPCG, tolPCG
/// Output: Sigma_iY, Sigma_iX, cov, alpha, eta
static void get_coefficients(const dvec &Y, const dmat &X, const dvec &w,
	const dvec &tau, int maxiterPCG, double tolPCG,
	dvec &Sigma_iY, dmat &Sigma_iX, dmat &cov, dvec &alpha, dvec &eta)
{
	int n_col_X = X.n_cols;
	Sigma_iY = get_PCG_diag_sigma(w, tau, Y, maxiterPCG, tolPCG);
	Sigma_iX.resize(Geno_NumSamp, n_col_X);
	dvec xv_i;
	for(int i = 0; i < n_col_X; i++)
	{
		xv_i = X.col(i);
		Sigma_iX.col(i) = get_PCG_diag_sigma(w, tau, xv_i, maxiterPCG, tolPCG);
	}
	cov = inv_sympd(X.t() * Sigma_iX);
	alpha = cov * (Sigma_iX.t() * Y);
	eta = Y - tau(0) * (Sigma_iY - Sigma_iX * alpha) / w;
}


// ========================================================================= //

static List get_coeff(const dvec &y, const dmat &X, const dvec &tau,
	const List &family, const dvec &alpha0, const dvec &eta0,
	const dvec &offset, int maxiterPCG, int maxiter, double tolPCG,
	bool verbose)
{
	// initialize
	const double tol_coef = 0.1;
	Function fc_linkinv = wrap(family["linkinv"]);
	Function fc_mu_eta = wrap(family["mu.eta"]);
	Function fc_variance = wrap(family["variance"]);
	dvec mu = as<dvec>(fc_linkinv(eta0));
	dvec mu_eta = as<dvec>(fc_mu_eta(eta0));
	dvec Y = eta0 - offset + (y - mu)/mu_eta;
	dvec v = as<dvec>(fc_variance(mu));
	dvec W = (mu_eta % mu_eta) / v;

	dvec Sigma_iY, alpha, eta, a0=alpha0;
	dmat Sigma_iX, cov;

	// iteration
	for (int i=0; i < maxiter; i++)
	{
		get_coefficients(Y, X, W, tau, maxiterPCG, tolPCG,
			Sigma_iY, Sigma_iX, cov, alpha, eta);
		eta += offset;
		mu = as<dvec>(fc_linkinv(eta));
		mu_eta = as<dvec>(fc_mu_eta(eta));
		Y = eta0 - offset + (y - mu)/mu_eta;
		v = as<dvec>(fc_variance(mu));
		W = (mu_eta % mu_eta) / v;

		if (max(abs(alpha - a0)/(abs(alpha) + abs(a0) + tol_coef)) < tol_coef)
			break;
		a0 = alpha;

Rprintf("iter: %d\n", i);
	}

	// output
	return List::create(
		_["Y"]   = Y,   _["mu"] = mu, _["alpha"] = alpha,
		_["eta"] = eta, _["W"]  = W,  _["cov"]   = cov,
		_["Sigma_iY"] = Sigma_iY,
		_["Sigma_iX"] = Sigma_iX);
}


// Functon to get working vector and fixed & random coefficients
// Run iterations to get converged alpha and eta
RcppExport SEXP saige_get_coeff(SEXP r_y, SEXP r_X, SEXP r_tau, SEXP r_family,
	SEXP r_alpha0, SEXP r_eta0, SEXP r_offset, SEXP r_maxiterPCG,
	SEXP r_maxiter, SEXP r_tolPCG, SEXP r_verbose)
{
BEGIN_RCPP
	RObject rcpp_result_gen;
	RNGScope rcpp_rngScope_gen;
	traits::input_parameter< const dvec& >::type y(r_y);
	traits::input_parameter< const dmat& >::type X(r_X);
	traits::input_parameter< const dvec& >::type tau(r_tau);
	traits::input_parameter< const List& >::type family(r_family);
	traits::input_parameter< const dvec& >::type alpha0(r_alpha0);
	traits::input_parameter< const dvec& >::type eta0(r_eta0);
	traits::input_parameter< const dvec& >::type offset(r_offset);
	int maxiterPCG = Rf_asInteger(r_maxiterPCG);
	int maxiter = Rf_asInteger(r_maxiter);
	double tolPCG = Rf_asReal(r_tolPCG);
	bool verbose = Rf_asLogical(r_verbose)==TRUE;
	rcpp_result_gen = wrap(get_coeff(y, X, tau, family, alpha0, eta0,
		offset, maxiterPCG, maxiter, tolPCG, verbose));
	return rcpp_result_gen;
END_RCPP
}


// ========================================================================= //

// Modified that (Sigma_iY, Sigma_iX, cov) are input parameters. Previously they are calculated in the function
//      This function needs the function getPCG1ofSigmaAndVector and function getCrossprod and GetTrace
static List get_AI_score(const dvec &Y, const dmat &X, const dvec &w,
	const dvec &tau, const dvec &Sigma_iY, const dmat &Sigma_iX, const dmat &cov,
	int nrun, int maxiterPCG, double tolPCG, double traceCVcutoff)
{
	dmat Sigma_iXt = Sigma_iX.t();
	dvec PY = Sigma_iY - Sigma_iX * (cov * (Sigma_iXt * Y));
	dvec APY;
	get_crossprod_b_grm(PY, APY);
	double YPAPY = dot(PY, APY);

	double Trace = get_trace(Sigma_iX, X, w, tau, cov, nrun, maxiterPCG,
		tolPCG, traceCVcutoff);
	dvec PAPY_1 = get_PCG_diag_sigma(w, tau, APY, maxiterPCG, tolPCG);
	dvec PAPY = PAPY_1 - Sigma_iX * (cov * (Sigma_iXt * PAPY_1));
	double AI = dot(APY, PAPY);

	return List::create(Named("YPAPY") = YPAPY, Named("Trace") = Trace,
		Named("PY") = PY, Named("AI") = AI);
}


RcppExport SEXP saige_get_AI_score(SEXP r_Y, SEXP r_X, SEXP r_w,
	SEXP r_tau, SEXP r_Sigma_iY, SEXP r_Sigma_iX, SEXP r_cov,
	SEXP r_nrun, SEXP r_maxiterPCG, SEXP r_tolPCG, SEXP r_traceCVcutoff)
{
BEGIN_RCPP
	RObject rcpp_result_gen;
	RNGScope rcpp_rngScope_gen;
	traits::input_parameter< const dvec& >::type Y(r_Y);
	traits::input_parameter< const dmat& >::type X(r_X);
	traits::input_parameter< const dvec& >::type w(r_w);
	traits::input_parameter< const dvec& >::type tau(r_tau);
	traits::input_parameter< const dvec& >::type Sigma_iY(r_Sigma_iY);
	traits::input_parameter< const dmat& >::type Sigma_iX(r_Sigma_iX);
	traits::input_parameter< const dmat& >::type cov(r_cov);
	int nrun = Rf_asInteger(r_nrun);
	int maxiterPCG = Rf_asInteger(r_maxiterPCG);
	double tolPCG = Rf_asReal(r_tolPCG);
	double traceCVcutoff = Rf_asReal(r_traceCVcutoff);
	rcpp_result_gen = wrap(get_AI_score(Y, X, w, tau, Sigma_iY, Sigma_iX, cov,
		nrun, maxiterPCG, tolPCG, traceCVcutoff));
	return rcpp_result_gen;
END_RCPP
}


// ========================================================================= //

// Modified that (Sigma_iY, Sigma_iX, cov) are input parameters. Previously they are calculated in the function
// This function needs the function get_PCG_diag_sigma and function get_crossprod, get_AI_score
static List fitglmmaiRPCG(const dvec &Y, const dmat &X, const dvec &w,
	const dvec &in_tau, const dvec &Sigma_iY, const dmat &Sigma_iX,
	const dmat &cov,
	int nrun, int maxiterPCG, double tolPCG, double tol, double traceCVcutoff)
{
	List re = get_AI_score(Y, X, w, in_tau, Sigma_iY, Sigma_iX, cov,
		nrun, maxiterPCG, tolPCG, traceCVcutoff);
  	double YPAPY = re["YPAPY"];
  	double Trace = re["Trace"];
  	double score1 = YPAPY - Trace;
  	double AI1 = re["AI"];
  	double Dtau = score1/AI1;
  	dvec tau = in_tau;
  	dvec tau0 = in_tau;
  	tau(1) = tau0(1) + Dtau;

  	for(int i=0; i<tau.n_elem; ++i)
  	{
		if (tau(i) < tol) tau(i) = 0;
  	}

  	double step = 1.0;
  	while (tau(1) < 0.0)
  	{
		step = step*0.5;
		tau(1) = tau0(1) + step * Dtau;
  	}

  	for(int i=0; i<tau.n_elem; ++i)
  	{
		if (tau(i) < tol) tau(i) = 0;
  	}

  	return List::create(Named("tau") = tau);
}


RcppExport SEXP saige_fitglmmaiRPCG(SEXP r_Y, SEXP r_X, SEXP r_w, SEXP r_tau,
	SEXP r_Sigma_iY, SEXP r_Sigma_iX, SEXP r_cov, SEXP r_nrun,
	SEXP r_maxiterPCG, SEXP r_tolPCG, SEXP r_tol, SEXP r_traceCVcutoff)
{
BEGIN_RCPP
	RObject rcpp_result_gen;
	RNGScope rcpp_rngScope_gen;
	traits::input_parameter< const dvec& >::type Y(r_Y);
	traits::input_parameter< const dmat& >::type X(r_X);
	traits::input_parameter< const dvec& >::type w(r_w);
	traits::input_parameter< const dvec& >::type tau(r_tau);
	traits::input_parameter< const dvec& >::type Sigma_iY(r_Sigma_iY);
	traits::input_parameter< const dmat& >::type Sigma_iX(r_Sigma_iX);
	traits::input_parameter< const dmat& >::type cov(r_cov);
	int nrun = Rf_asInteger(r_nrun);
	int maxiterPCG = Rf_asInteger(r_maxiterPCG);
	double tolPCG = Rf_asReal(r_tolPCG);
	double tol = Rf_asReal(r_tol);
	double traceCVcutoff = Rf_asReal(r_traceCVcutoff);
	rcpp_result_gen = wrap(fitglmmaiRPCG(Y, X, w, tau, Sigma_iY, Sigma_iX,
		cov, nrun, maxiterPCG, tolPCG, tol, traceCVcutoff));
	return rcpp_result_gen;
END_RCPP
}
