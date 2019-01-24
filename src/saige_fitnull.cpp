
#if defined(__clang__)
#pragma clang optimize on
#elif defined(__GNUC__) && ((__GNUC__>4) || (__GNUC__==4 && __GNUC_MINOR__>=4))
#pragma GCC optimize("O3")
#endif


#include <RcppArmadillo.h>
#include <RcppParallel.h>
#include <tbb/parallel_for.h>
#include <algorithm>
#include "vectorization.h"


using namespace std;
using namespace Rcpp;
using namespace arma;
using namespace RcppParallel;


// ========================================================================= //
// Define Intel TBB macro with a thread index starting from 0

#if RCPP_PARALLEL_USE_TBB

#define PARALLEL_HEAD(SIZE)    \
	tbb::parallel_for(tbb::blocked_range<size_t>(0, SIZE,  \
		SIZE/NumThreads + (SIZE % NumThreads ? 1 : 0)),  \
		[&](const tbb::blocked_range<size_t> &r)  \
	{  \
		const int th_idx = tbb::this_task_arena::current_thread_index();  \
		if (th_idx < 0 || th_idx >= NumThreads)  \
			throw "Invalid tbb::this_task_arena::current_thread_index()!";

#define PARALLEL_FOR(i, SIZE)    \
		PARALLEL_HEAD(SIZE)    \
		for (size_t i=r.begin(); i < r.end(); i++)

#define PARALLEL_RANGE(st, ed, SIZE)    \
		PARALLEL_HEAD(SIZE)    \
		const size_t st = r.begin(), ed = r.end();

#define PARALLEL_END    });

#else

#define PARALLEL_FOR(i, SIZE)    \
		const int th_idx = 0;  \
		for (size_t i=0; i < SIZE; i++)
#define PARALLEL_RANGE(st, ed, SIZE)    \
		const int th_idx = 0;  \
		const size_t st = 0, ed = SIZE;
#define PARALLEL_END

#endif



// ========================================================================= //
// R functions for random numbers

inline static void set_seed(unsigned int seed)
{
	Environment base_env("package:base");
	Function set_seed_r = base_env["set.seed"];
	set_seed_r(seed);
}

inline static NumericVector random_binary(int n)
{
	return(rbinom(n, 1, 0.5));
}


// ========================================================================= //
// Store packed SNP genotypes

static int NumThreads = 0;  //< the number of threads

static unsigned char *Geno_PackedRaw = NULL;
static size_t Geno_NumSamp = 0;
static size_t Geno_PackedNumSamp = 0;
static size_t Geno_NumVariant = 0;

static double *buf_std_geno = NULL;   //< a 4-by-n_variant look-up matrix
static double *buf_diag_grm = NULL;   //< n_samp-length, sigma_i = sum_j adj.g[i,j]^2
static double *buf_crossprod = NULL;  //< nThread-by-n_samp matrix

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
RcppExport SEXP saige_store_geno(SEXP rawgeno, SEXP num_samp, SEXP r_buf_geno,
	SEXP r_buf_sigma, SEXP r_buf_crossprod)
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
	buf_std_geno = REAL(r_buf_geno);
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
	buf_diag_grm = REAL(r_buf_sigma);
	memset(buf_diag_grm, 0, sizeof(double)*Geno_NumSamp);
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

	// set the buffer for get_crossprod_b_grm()
	NumericMatrix mat(r_buf_crossprod);
	buf_crossprod = REAL(r_buf_crossprod);
	NumThreads = mat.ncol();
	if (NumThreads > Geno_NumSamp) NumThreads = Geno_NumSamp;
	if (NumThreads > Geno_NumVariant) NumThreads = Geno_NumVariant;
	if (NumThreads < 1) NumThreads = 1;

END_RCPP
}



// ========================================================================= //

/// Cross-product of standardized genotypes and a numeric vector
/// Input: b (n_samp-length)
/// Output: out_b (n_samp-length)
static void get_crossprod_b_grm(const dcolvec &b, dvec &out_b)
{
	// initialize
	memset(buf_crossprod, 0, sizeof(double)*Geno_NumSamp*NumThreads);

	PARALLEL_FOR(i, Geno_NumVariant)
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
		double *pbb = buf_crossprod + Geno_NumSamp * th_idx;
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
	PARALLEL_END

	// initialize out_b
	out_b.resize(Geno_NumSamp);
	PARALLEL_RANGE(st, ed, Geno_NumSamp)
		size_t len = ed - st;
		const double *s = buf_crossprod + st;
		double *p = &out_b[st];
		memset(p, 0, sizeof(double)*len);
		for (int i=0; i < NumThreads; i++)
		{
			f64_add(len, s, p);
			s += Geno_NumSamp;
		}
		f64_mul(len, 1.0/Geno_NumVariant, p);
	PARALLEL_END
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
static dvec get_crossprod(const dcolvec &b, const dvec& w, const dvec& tau)
{
	const double tau0 = tau[0];
	const double tau1 = tau[1];
	if (tau1 == 0)
	{
		return(tau0 * (b % (1/w)));
	} else {
		dvec out_b;
		get_crossprod_b_grm(b, out_b);
		return(tau0 * (b % (1/w)) + tau1 * out_b);
	}
}


/// Sigma = tau[0] * diag(1/W) + tau[1] * grm
/// Input: wVec, tauVec, bVec, maxiterPCG, tolPCG
static dvec get_PCG_diag_sigma(const dvec &w, const dvec &tau, const dvec &b,
	int maxiterPCG, double tolPCG)
{
	dvec r = b, r1, minv;
	get_diag_sigma(w, tau, minv);
	minv = 1 / minv;
	double sumr2 = sum(r % r);

	dvec z = minv % r, z1;
	dvec p = z;
	dvec x(Geno_NumSamp);
	x.zeros();

	int iter = 0;
	while (sumr2 > tolPCG && iter < maxiterPCG)
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

		sumr2 = sum(r % r);
	}

	if (iter >= maxiterPCG)
	{
		cout << "PCG does not converge. You may increase maxiter number." << endl;
	}
	return(x);
}


/// Calculate the coefficient of variation for mean of a vector
static double calcCV(const dvec& x)
{
	double x_mean = mean(x);
	double x_sd = stddev(x);
	double vecCV = (x_sd / x_mean) / int(x.n_elem);
	return(vecCV);
}


// [[export]]
static double get_trace(const dmat &Sigma_iX, const dmat& X, const dvec& w,
	const dvec& tau, const dmat& cov,
	int nrun, int maxiterPCG, double tolPCG, double traceCVcutoff)
{
	set_seed(200);
	dmat Sigma_iXt = Sigma_iX.t();
	dvec Sigma_iu;  
	dcolvec Pu;
	dvec Au, u;

	int nrunStart = 0;
	int nrunEnd = nrun;
	double traceCV = traceCVcutoff + 0.1;
	dvec buf(nrun);
	buf.zeros();

	while(traceCV > traceCVcutoff)
	{
		for(int i = nrunStart; i < nrunEnd; i++)
		{
			u = as<dvec>(random_binary(Geno_NumSamp));
			u = 2*u - 1;
			Sigma_iu = get_PCG_diag_sigma(w, tau, u, maxiterPCG, tolPCG);
			Pu = Sigma_iu - Sigma_iX * (cov *  (Sigma_iXt * u));
			get_crossprod_b_grm(u, Au);
			buf(i) = dot(Au, Pu);
		}
		traceCV = calcCV(buf);
		if(traceCV > traceCVcutoff)
		{
			nrunStart = nrunEnd;
			nrunEnd = nrunEnd + 10;
			buf.resize(nrunEnd);
			cout << "CV for trace random estimator using "<< nrun <<
				" runs is " << traceCV <<  " > " << traceCVcutoff << endl;
			cout << "try " << nrunEnd << " runs" << endl;
		}
	}

	return(mean(buf));
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

	// iteration
	dvec a0 = alpha0;
	for (int i=0; i < maxiter; i++)
	{
		get_coefficients(Y, X, W, tau, maxiterPCG, tolPCG,
			Sigma_iY, Sigma_iX, cov, alpha, eta);

		eta += offset;
		mu = as<dvec>(fc_linkinv(eta));
		mu_eta = as<dvec>(fc_mu_eta(eta));
		Y = eta - offset + (y - mu)/mu_eta;
		W = (mu_eta % mu_eta) / as<dvec>(fc_variance(mu));

		if (max(abs(alpha - a0)/(abs(alpha) + abs(a0) + tol_coef)) < tol_coef)
			break;
		a0 = alpha;
	}
}


// Modified that (Sigma_iY, Sigma_iX, cov) are input parameters. Previously they are calculated in the function
//      This function needs the function getPCG1ofSigmaAndVector and function getCrossprod and GetTrace
static void get_AI_score(const dvec &Y, const dmat &X, const dvec &w,
	const dvec &tau, const dvec &Sigma_iY, const dmat &Sigma_iX, const dmat &cov,
	int nrun, int maxiterPCG, double tolPCG, double traceCVcutoff,
	double &YPAPY, double &Trace, double &AI)
{
	dmat Sigma_iXt = Sigma_iX.t();
	dvec PY = Sigma_iY - Sigma_iX * (cov * (Sigma_iXt * Y));
	dvec APY;
	get_crossprod_b_grm(PY, APY);
	YPAPY = dot(PY, APY);

	Trace = get_trace(Sigma_iX, X, w, tau, cov, nrun, maxiterPCG,
		tolPCG, traceCVcutoff);
	dvec PAPY_1 = get_PCG_diag_sigma(w, tau, APY, maxiterPCG, tolPCG);
	dvec PAPY = PAPY_1 - Sigma_iX * (cov * (Sigma_iXt * PAPY_1));
	AI = dot(APY, PAPY);
}


// Modified that (Sigma_iY, Sigma_iX, cov) are input parameters. Previously they are calculated in the function
// This function needs the function get_PCG_diag_sigma and function get_crossprod, get_AI_score
static dvec fitglmmaiRPCG(const dvec &Y, const dmat &X, const dvec &w,
	const dvec &in_tau, const dvec &Sigma_iY, const dmat &Sigma_iX,
	const dmat &cov,
	int nrun, int maxiterPCG, double tolPCG, double tol, double traceCVcutoff)
{
	double YPAPY, trace, AI;
	get_AI_score(Y, X, w, in_tau, Sigma_iY, Sigma_iX, cov, nrun,
		maxiterPCG, tolPCG, traceCVcutoff,
		YPAPY, trace, AI);
  	double score = YPAPY - trace;
  	double Dtau = score / AI;
  	dvec tau = in_tau;
  	dvec tau0 = in_tau;
  	tau[1] = tau0[1] + Dtau;

  	for(int i=0; i < tau.n_elem; i++)
		if (tau[i] < tol) tau[i] = 0;

  	double step = 1.0;
  	while (tau[1] < 0.0)
  	{
		step *= 0.5;
		tau[1] = tau0[1] + step * Dtau;
  	}

  	for(int i=0; i < tau.n_elem; i++)
		if (tau[i] < tol) tau[i] = 0;

  	return tau;
}


// ========================================================================= //

inline static void print_vec(const char *s, dvec &x)
{
	Rprintf("%s(", s);
	for (int i=0; i < x.n_elem; i++)
	{
		if (i > 0) Rprintf(", ");
		Rprintf("%0.7g", x[i]);
	}
	Rprintf(")\n");
}


RcppExport SEXP saige_fit_AI_PCG_binary(SEXP r_fit0, SEXP r_X, SEXP r_tau,
	SEXP r_param)
{
BEGIN_RCPP

	// parameters for fitting the model
	List param(r_param);
	const double tol = Rf_asReal(param["tol"]);
	const double tol_inv_2 = 1 / (tol*tol);
	const double tolPCG = Rf_asReal(param["tolPCG"]);
	const int maxiter = Rf_asInteger(param["maxiter"]);
	const int maxiterPCG = Rf_asInteger(param["maxiterPCG"]);
	const int nrun = Rf_asInteger(param["nrun"]);
	const double traceCVcutoff = Rf_asReal(param["traceCVcutoff"]);
	const bool verbose = Rf_asLogical(param["verbose"])==TRUE;

	List fit0(r_fit0);
	dvec y = as<dvec>(fit0["y"]);
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

	dvec alpha;
	dmat cov;

	dvec tau = as<dvec>(r_tau);
	dvec tau0 = tau;

	dvec re_Y, re_mu, re_alpha, re_eta, re_W, re_Sigma_iY;
	dmat re_cov, re_Sigma_iX;
	get_coeff(y, X, tau, family, alpha0, eta0, offset, maxiterPCG, maxiter,
		tolPCG, verbose,
		re_Y, re_mu, re_alpha, re_eta, re_W, re_cov, re_Sigma_iY, re_Sigma_iX);

	double YPAPY, Trace, AI;
	get_AI_score(re_Y, X, re_W, tau, re_Sigma_iY, re_Sigma_iX, re_cov, nrun,
		maxiterPCG, tolPCG, traceCVcutoff,
		YPAPY, Trace, AI);

	tau[1] = std::max(0.0, tau0[1] + tau0[1]*tau0[1]*(YPAPY - Trace)/y.size());
	if (verbose)
		print_vec("Variance component estimates: ", tau);

	int iter = 1;
	for (; iter <= maxiter; iter++)
	{
		if (verbose)
		{
			Rprintf("Iteration %d:\n", iter);
			print_vec("    tau: ", tau);
			print_vec("    fixed coeff: ", alpha);
		}

		alpha0 = re_alpha;
		tau0 = tau;
		eta0 = eta;
		get_coeff(y, X, tau, family, alpha0, eta0, offset, maxiterPCG, maxiter,
			tolPCG, verbose,
			re_Y, re_mu, re_alpha, re_eta, re_W, re_cov, re_Sigma_iY, re_Sigma_iX);
		tau = fitglmmaiRPCG(re_Y, X, re_W, tau, re_Sigma_iY,
			re_Sigma_iX, re_cov, nrun, maxiterPCG, tolPCG,
			tol, traceCVcutoff);

		cov = re_cov; alpha = re_alpha; eta = re_eta;
		Y = re_Y; mu = re_mu;

		if (tau[1] == 0) break;
		if (max(abs(tau-tau0)/(abs(tau)+abs(tau0)+tol)) < tol) break;
		if (max(tau) > tol_inv_2)
		{
			Rprintf("Large variance estimate observed in the iterations, model not converged ...");
			iter = maxiter + 1;
			break;
		}
	}

	get_coeff(y, X, tau, family, alpha0, eta0, offset, maxiterPCG, maxiter,
		tolPCG, verbose,
		re_Y, re_mu, re_alpha, re_eta, re_W, re_cov, re_Sigma_iY, re_Sigma_iX);
	cov = re_cov; alpha = re_alpha; eta = re_eta;
	Y = re_Y; mu = re_mu;

	if (verbose)
	{
		print_vec("Final tau: " ,tau);
		print_vec("    fixed coeff: ", alpha);
	}

	return List::create(
		_["theta"] = tau,
		_["coefficients"] = alpha,
		_["linear.predictors"] = eta,
		_["fitted.values"] = mu,
		_["Y"] = Y,
		_["residuals"] = y - mu,
		_["cov"] = cov,
		_["converged"] = bool(iter <= maxiter));

END_RCPP
}
