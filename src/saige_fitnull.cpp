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
			int b1 = (i << 2) & 0x03;
			int b2 = (i << 4) & 0x03;
			int b3 = (i << 6) & 0x03;
			num_valid[i] = (b0<3 ? 1:0) + (b1<3 ? 1:0) + (b2<3 ? 1:0) + (b3<3 ? 1:0);
			num_sum[i] = (b0<3 ? b0:0) + (b1<3 ? b1:0) + (b2<3 ? b2:0) + (b3<3 ? b3:0);
		}
		lookup_has_init = true;
	}
}


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
		unsigned char *g = (unsigned char*)&RawGeno[0] + Geno_PackedNumSamp*i;
		// calculate allele frequency
		int n_valid=0, sum=0;
		for (size_t j=0; j < Geno_PackedNumSamp; j++)
		{
			n_valid += num_valid[g[j]];
			sum += num_sum[g[j]];
		}
		double af = double(sum) / (2*n_valid);
		double inv = 1 / (af*(1-af));
		double *p = &buf_std_geno[4*i];
		p[0] = (0 - 2*af) * inv; p[1] = (1 - 2*af) * inv;
		p[2] = (2 - 2*af) * inv; p[3] = 0;
	}

	// calculate diag(grm)
	buf_diag_grm = REAL(buf_sigma);
	

END_RCPP
}



// ========================================================================= //

/// Cross-product of standardized genotypes and a numeric vector
/// Input: b (n_samp-length)
/// Output: out_b (n_samp-length)
static void get_crossprod_b_grm(dcolvec &b, dvec &out_b)
{
	out_b.resize(Geno_NumSamp);
	memset(&out_b[0], 0, sizeof(double)*Geno_NumSamp);

	for (size_t i=0; i < Geno_NumVariant; i++)
	{
		unsigned char *g = Geno_PackedRaw + Geno_PackedNumSamp*i;
		double *base = buf_std_geno + 4*i;

		// get dot = sum(std.geno .* b)
		double dot = 0, *pb = &b[0];
		size_t n = Geno_NumSamp;
		for (; n >= 4; n-=4, pb+=4)
		{
			unsigned char gg = *g++;
			dot += base[gg & 0x03] * pb[0] + base[(gg << 2) & 0x03] * pb[1] +
				base[(gg << 4) & 0x03] * pb[2] + base[gg << 6] * pb[3];
		}
		for (unsigned char gg = (n>0 ? *g : 0); n > 0; n--)
		{
			dot += base[gg & 0x03] * (*pb++);
			gg <<= 2;
		}

		// update the output rv += dot .* std.geno
		pb = &out_b[0];
		g = Geno_PackedRaw + Geno_PackedNumSamp*i;
		n = Geno_NumSamp;
		for (; n >= 4; n-=4, pb+=4, pb+=4)
		{
			unsigned char gg = *g++;
			pb[0] += dot * base[gg & 0x03];
			pb[1] += dot * base[(gg << 2) & 0x03];
			pb[2] += dot * base[(gg << 4) & 0x03];
			pb[3] += dot * base[gg << 6];
		}
		for (unsigned char gg = (n>0 ? *g : 0); n > 0; n--)
		{
			(*pb++) += dot * base[gg & 0x03];
			gg <<= 2;
		}
	}

	// normalize
	f64_mul(Geno_NumSamp, 1.0/Geno_NumVariant, &out_b[0]);
}

  
/// Sigma = tau[0] * diag(1/W) + tau[1] * diag(grm)
/// Input: w, tau
/// Output: out_sigma
static void get_diag_sigma(dvec& w, dvec& tau, dvec &out_sigma)
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
static dcolvec get_crossprod(dcolvec& b, dvec& w, dvec& tau)
{
	const double tau0 = tau[0];
	const double tau1 = tau[1];
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
static dvec get_PCG_diag_sigma(dvec& wVec, dvec& tauVec, dvec& bVec,
	int maxiterPCG, float tolPCG)
{
	dvec rVec = bVec;
	dvec r1Vec;
	int Nnomissing = geno.getNnomissing();

	dvec crossProdVec(Nnomissing);
	dvec minvVec = 1/get_diag_sigma(wVec, tauVec);
	float sumr2 = sum(rVec % rVec);

	dvec zVec = minvVec % rVec;
	dvec z1Vec;
	dvec pVec = zVec;

	dvec xVec(Nnomissing);
	xVec.zeros();

	int iter = 0;
	while (sumr2 > tolPCG && iter < maxiterPCG)
	{
		iter = iter + 1;
		dcolvec ApVec = get_crossprod(pVec, wVec, tauVec);
		dvec preA = (rVec.t() * zVec)/(pVec.t() * ApVec);

		float a = preA(0);
		xVec = xVec + a * pVec;
		r1Vec = rVec - a * ApVec;

		z1Vec = minvVec % r1Vec;
		dvec Prebet = (z1Vec.t() * r1Vec)/(zVec.t() * rVec);
		float bet = Prebet(0);
		pVec = z1Vec+ bet*pVec;
		zVec = z1Vec;
		rVec = r1Vec;

		sumr2 = sum(rVec % rVec);
	}

	if (iter >= maxiterPCG)
	{
		cout << "pcg did not converge. You may increase maxiter number." << endl;
	}
	// cout << "iter from get_PCG_diag_sigma " << iter << endl;
	return(xVec);
}


// [[Rcpp::export]]
Rcpp::NumericVector nb(int n) {
  	return(rbinom(n,1,0.5));
}

//This function calculates the coefficients of variation for mean of a vector
// [[Rcpp::export]]
float calCV(dvec& xVec){
  int veclen = xVec.n_elem;
  float vecMean = mean(xVec);
  float vecSd = stddev(xVec);
  float vecCV = (vecSd/vecMean)/veclen;
  return(vecCV);
}



// [[Rcpp::export]]
static double get_trace(dmat Sigma_iX, dmat& Xmat, dvec& wVec, dvec& tauVec,
	dmat& cov1, int nrun, int maxiterPCG, double tolPCG, double traceCVcutoff)
{
	set_seed(200);
	int Nnomissing = geno.getNnomissing();
	dmat Sigma_iXt = Sigma_iX.t();
	dvec Sigma_iu;  
	dcolvec Pu;
	dvec Au;
	dvec uVec;

	int nrunStart = 0;
	int nrunEnd = nrun;
	float traceCV = traceCVcutoff + 0.1;
	dvec tempVec(nrun);
	tempVec.zeros();

	while(traceCV > traceCVcutoff)
	{
		//dvec tempVec(nrun);
		//tempVec.zeros();
		for(int i = nrunStart; i < nrunEnd; i++)
		{
			Rcpp::NumericVector uVec0;
			uVec0 = nb(Nnomissing);
			uVec = as<dvec>(uVec0);
			uVec = uVec*2 - 1;
			Sigma_iu = get_PCG_diag_sigma(wVec, tauVec, uVec, maxiterPCG, tolPCG);
			Pu = Sigma_iu - Sigma_iX * (cov1 *  (Sigma_iXt * uVec));
			Au = get_crossprod_b_grm(uVec);
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
static void get_coefficients(dvec &Y, dmat &X, dvec &w, dvec &tau,
	int maxiterPCG, double tolPCG,
	dvec &Sigma_iY, dmat &Sigma_iX, dmat &cov, dvec &alpha, dvec &eta)
{
	int n_col_X = X.n_cols;
	Sigma_iY = get_PCG_diag_sigma(w, tau, Y, maxiterPCG, tolPCG);
	Sigma_iX.resize(Geno_NumSamp, n_col_X);
	dvec XmatVecTemp;
	for(int i = 0; i < n_col_X; i++)
	{
		XmatVecTemp = Xmat.col(i);
		Sigma_iX.col(i) = get_PCG_diag_sigma(wVec, tauVec, XmatVecTemp,
			maxiterPCG, tolPCG);
	}
	cov = inv_sympd(Xmat.t() * Sigma_iX);
	alpha = cov * (Sigma_iX.t() * Yvec);
	eta = Yvec - tauVec(0) * (Sigma_iY - Sigma_iX * alpha) / wVec;
}



// Modified by SLEE, 04/16/2017
// Modified that (Sigma_iY, Sigma_iX, cov) are input parameters. Previously they are calculated in the function
//      This function needs the function get_PCG_diag_sigma and function get_crossprod and get_trace
static List get_AI_score(dvec& Yvec, dmat& Xmat, dvec& wVec,  dvec& tauVec,
	dvec& Sigma_iY, dmat & Sigma_iX, dmat & cov,
	int nrun, int maxiterPCG, double tolPCG, double traceCVcutoff)
{
	dmat Sigma_iXt = Sigma_iX.t();

	dvec PY1 = Sigma_iY - Sigma_iX * (cov * (Sigma_iXt * Yvec));
	dvec APY = get_crossprod_b_grm(PY1);
	double YPAPY = dot(PY1, APY);

	double Trace = get_trace(Sigma_iX, Xmat, wVec, tauVec, cov, nrun, maxiterPCG, tolPCG, traceCVcutoff);
	dvec PAPY_1 = get_PCG_diag_sigma(wVec, tauVec, APY, maxiterPCG, tolPCG);
	dvec PAPY = PAPY_1 - Sigma_iX * (cov * (Sigma_iXt * PAPY_1));
	double AI = dot(APY, PAPY);

	return List::create(Named("YPAPY") = YPAPY, Named("Trace") = Trace,
		Named("PY") = PY1, Named("AI") = AI);
}


// Modified by SLEE, 04/16/2017
// Modified that (Sigma_iY, Sigma_iX, cov) are input parameters. Previously they are calculated in the function
// This function needs the function get_PCG_diag_sigma and function get_crossprod, get_AI_score
List fitglmmaiRPCG(dvec& Yvec, dmat& Xmat, dvec& wVec,  dvec& tauVec,
	dvec& Sigma_iY, dmat & Sigma_iX, dmat & cov,
	int nrun, int maxiterPCG, double tolPCG, double tol, double traceCVcutoff)
{
  	Rcpp::List re = get_AI_score(Yvec, Xmat,wVec,  tauVec, Sigma_iY, Sigma_iX, cov, nrun, maxiterPCG, tolPCG, traceCVcutoff);
  	double YPAPY = re["YPAPY"];
  	double Trace = re["Trace"];
  	double score1 = YPAPY - Trace;
  	double AI1 = re["AI"];
  	double Dtau = score1/AI1;
  	dvec tau0 = tauVec;
  	tauVec(1) = tau0(1) + Dtau;

  	for(int i=0; i<tauVec.n_elem; ++i)
  	{
		if (tauVec(i) < tol) tauVec(i) = 0;
  	}

  	double step = 1.0;
  	while (tauVec(1) < 0.0)
  	{
		step = step*0.5;
		tauVec(1) = tau0(1) + step * Dtau;
  	}

  	for(int i=0; i<tauVec.n_elem; ++i)
  	{
		if (tauVec(i) < tol) tauVec(i) = 0;
  	}
  	return List::create(Named("tau") = tauVec);
}



/*add for SPA by Wei 04222017*/
// [[Rcpp::export]]
dmat getSigma_X(dvec& wVec, dvec& tauVec,dmat& Xmat, int maxiterPCG, float tolPCG){


  	int Nnomissing = Xmat.n_rows;
  	int colNumX = Xmat.n_cols;

  	cout << colNumX << endl;
  	cout << size(wVec) << endl;
  	cout << size(tauVec) << endl;


  	dmat Sigma_iX1(Nnomissing,colNumX);
  	dvec XmatVecTemp;

  	for(int i = 0; i < colNumX; i++){
    		XmatVecTemp = Xmat.col(i);
    		Sigma_iX1.col(i) = get_PCG_diag_sigma(wVec, tauVec, XmatVecTemp, maxiterPCG, tolPCG);
  	}
  	return(Sigma_iX1);
}


// [[Rcpp::export]]
dvec  getSigma_G(dvec& wVec, dvec& tauVec,dvec& Gvec, int maxiterPCG, float tolPCG){
  	dvec Sigma_iG;
  	Sigma_iG = get_PCG_diag_sigma(wVec, tauVec, Gvec, maxiterPCG, tolPCG);
  	return(Sigma_iG);
}


//This function needs the function get_PCG_diag_sigma and function get_crossprod_b_grm
// [[Rcpp::export]]
dvec GetTrace_q(dmat Sigma_iX, dmat& Xmat, dvec& wVec, dvec& tauVec, dmat& cov1,  int nrun, int maxiterPCG, float tolPCG, float traceCVcutoff){
  	set_seed(200);
  	dmat Sigma_iXt = Sigma_iX.t();
  	int Nnomissing = geno.getNnomissing();
  	dvec tempVec(nrun);
  	tempVec.zeros();
  	dvec tempVec0(nrun);
  	tempVec0.zeros();

        dvec Sigma_iu;
        dcolvec Pu;
        dvec Au;
        dvec uVec;

        int nrunStart = 0;
        int nrunEnd = nrun;
        float traceCV = traceCVcutoff + 0.1;
        float traceCV0 = traceCVcutoff + 0.1;

        while((traceCV > traceCVcutoff) | (traceCV0 > traceCVcutoff)){


  	for(int i = nrunStart; i < nrunEnd; i++){

    		Rcpp::NumericVector uVec0;
    		uVec0 = nb(Nnomissing);
    		uVec = as<dvec>(uVec0);
    		uVec = uVec*2 - 1;
  //  		dvec Sigma_iu;
    		Sigma_iu = get_PCG_diag_sigma(wVec, tauVec, uVec, maxiterPCG, tolPCG);
  //  		dcolvec Pu;
    		Pu = Sigma_iu - Sigma_iX * (cov1 *  (Sigma_iXt * uVec));
  //  		dvec Au;
    		Au = get_crossprod_b_grm(uVec);
    		tempVec(i) = dot(Au, Pu);
    		tempVec0(i) = dot(uVec, Pu);
                Au.clear();
      		Pu.clear();
      		Sigma_iu.clear();
      		uVec.clear();
  	}
	traceCV = calCV(tempVec);
	traceCV0 = calCV(tempVec0);
	
	if((traceCV > traceCVcutoff) | (traceCV0 > traceCVcutoff)){
          nrunStart = nrunEnd;
          nrunEnd = nrunEnd + 10;
          tempVec.resize(nrunEnd);
          tempVec0.resize(nrunEnd);
          cout << "CV for trace random estimator using "<< nrun << " runs is " << traceCV <<  "(> " << traceCVcutoff << endl;
          cout << "try " << nrunEnd << "runs" << endl;
        }

    }   

  	dvec traVec(2);
  	traVec(1) = mean(tempVec);
  	traVec(0) = mean(tempVec0);
	tempVec.clear();
	tempVec0.clear();
  	return(traVec);
}



// ========================================================================= //

// getCoefficients
RcppExport SEXP _SAIGE_getCoefficients(SEXP YvecSEXP, SEXP XmatSEXP, SEXP wVecSEXP, SEXP tauVecSEXP, SEXP maxiterPCGSEXP, SEXP tolPCGSEXP)
{
BEGIN_RCPP
	Rcpp::RObject rcpp_result_gen;
	Rcpp::RNGScope rcpp_rngScope_gen;
	Rcpp::traits::input_parameter< dvec& >::type Yvec(YvecSEXP);
	Rcpp::traits::input_parameter< dmat& >::type Xmat(XmatSEXP);
	Rcpp::traits::input_parameter< dvec& >::type wVec(wVecSEXP);
	Rcpp::traits::input_parameter< dvec& >::type tauVec(tauVecSEXP);
	Rcpp::traits::input_parameter< int >::type maxiterPCG(maxiterPCGSEXP);
	Rcpp::traits::input_parameter< float >::type tolPCG(tolPCGSEXP);
	// rcpp_result_gen = Rcpp::wrap(get_coefficients(Yvec, Xmat, wVec, tauVec, maxiterPCG, tolPCG));
	// return rcpp_result_gen;
END_RCPP
}


// Functon to get working vector and fixed & random coefficients
// Run iterations to get converged alpha and eta
Get_Coef <- function(y, X, tau, family, alpha0, eta0, offset, maxiterPCG,
    tolPCG, maxiter, verbose=FALSE)
{
    tol.coef <- 0.1
    mu <- family$linkinv(eta0)
    mu.eta <- family$mu.eta(eta0)
    Y <- eta0 - offset + (y - mu)/mu.eta
    sqrtW <- mu.eta / sqrt(family$variance(mu))
    W <- sqrtW^2

    for(i in 1:maxiter)
    {
        cat("iGet_Coef: ", i, "\n", sep="")
        re.coef <- getCoefficients(Y, X, W, tau, maxiterPCG, tol=tolPCG)
        alpha <- re.coef$alpha
        eta <- re.coef$eta + offset
        if (verbose)
        {
            cat("Tau: "); print(tau)
            cat("Fixed-effect coefficients: "); print(alpha)
        }
        mu <- family$linkinv(eta)
        mu.eta <- family$mu.eta(eta)
        Y <- eta - offset + (y - mu)/mu.eta
        sqrtW <- mu.eta/sqrt(family$variance(mu))
        W <- sqrtW^2
        if (max(abs(alpha-alpha0)/(abs(alpha)+abs(alpha0)+tol.coef)) < tol.coef)
            break
        alpha0 = alpha
    }

    list(Y=Y, alpha=alpha, eta=eta, W=W, cov=re.coef$cov, sqrtW=sqrtW,
        Sigma_iY=re.coef$Sigma_iY, Sigma_iX=re.coef$Sigma_iX, mu=mu)
}


