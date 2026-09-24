// ===========================================================================
//
// saige_surv.cpp: survival (time-to-event) support for SAIGEgds
//
// Cox proportional-hazards frailty model via the Breslow Cox<->Poisson
// equivalence, ported from GATE (https://github.com/weizhou0/GATE, GPL>=2).
//
// The numerical kernels live here:
//   (1) Breslow baseline cumulative hazard  Lambda0(eta, time, status)
//   (2) saddlepoint approximation for the weighted-Poisson score statistic
//   (3) WminusU::build(), the setup for the Cox risk-set covariance operator
//   (4) the Cox model without random effects (independent samples)
//
// Copyright (C) 2026    Xiuwen Zheng
// License: GPL-3
//

#include <RcppArmadillo.h>
#include "saige.h"
#include "saige_surv.h"
#include <vector>
#include <algorithm>
#include <cmath>

using namespace std;
using namespace Rcpp;


namespace SAIGE_SURV
{

// ===========================================================================
// (1) Breslow baseline cumulative hazard
//
// Given linear predictor eta, event time, and event status (1=event, 0=censor),
// returns Lambda0_i = baseline cumulative hazard at subject i's event time, so
// that mu_i = Lambda0_i * exp(eta_i) is the Poisson mean (expected # events).
// Matches GATE's GetLambda0(): ties handled by the Breslow approximation
// (hazard increment d_t / sum_{risk set} exp(eta), risk set = {time >= t}).

void breslow_lambda0(size_t n, const double eta[], const double time[],
	const int status[], double Lambda0[])
{
	// order indices by ascending time (stable, so equal times keep input order)
	std::vector<size_t> ord(n);
	for (size_t i=0; i < n; i++) ord[i] = i;
	std::stable_sort(ord.begin(), ord.end(),
		[&](size_t a, size_t b){ return time[a] < time[b]; });

	// reverse cumulative sum of exp(eta) in sorted order:
	//   S[r] = sum_{r' >= r} exp(eta[ord[r']])  (risk set for time >= time[ord[r]])
	std::vector<double> S(n);
	double acc = 0;
	for (size_t r=n; r-- > 0; )
	{
		acc += std::exp(eta[ord[r]]);
		S[r] = acc;
	}

	// walk distinct-time blocks from smallest to largest, accumulating hazard
	double cumHaz = 0;
	size_t r = 0;
	while (r < n)
	{
		// [r, r2) is the block of identical times
		size_t r2 = r + 1;
		while (r2 < n && time[ord[r2]] == time[ord[r]]) r2++;
		// number of events in this time block
		int d = 0;
		for (size_t k=r; k < r2; k++) d += status[ord[k]];
		// Breslow increment (risk-set sum uses the first row of the block)
		if (d > 0 && S[r] > 0)
			cumHaz += d / S[r];
		// assign cumulative hazard to every subject in this time block
		for (size_t k=r; k < r2; k++)
			Lambda0[ord[k]] = cumHaz;
		r = r2;
	}
}


// ===========================================================================
// (2) Saddlepoint approximation for weighted sum of Poisson variables
//
// Score statistic  S = sum_i (y_i - mu_i) * g_i, where N_i ~ Poisson(mu_i).
// CGF of sum g_i*(N_i - mu_i):
//   K(t)  = sum mu_i (exp(g_i t) - g_i t - 1)         (centered)
//   K1(t) = sum mu_i g_i (exp(g_i t) - 1)
//   K2(t) = sum mu_i g_i^2 exp(g_i t)

static inline double Korg_Poi(double t, size_t n, const double mu[],
	const double g[])
{
	double s = 0;
	for (size_t i=0; i < n; i++)
	{
		double gt = g[i] * t;
		s += mu[i] * (std::exp(gt) - gt - 1.0);
	}
	return s;
}

// K1(t) - q
static inline double K1_adj_Poi(double t, size_t n, const double mu[],
	const double g[], double q)
{
	double s = 0;
	for (size_t i=0; i < n; i++)
		s += mu[i] * g[i] * (std::exp(g[i] * t) - 1.0);
	return s - q;
}

static inline double K2_Poi(double t, size_t n, const double mu[],
	const double g[])
{
	double s = 0;
	for (size_t i=0; i < n; i++)
		s += mu[i] * g[i] * g[i] * std::exp(g[i] * t);
	return s;
}

// Newton root finder for K1_adj_Poi(t) = 0, with the bisection-style safeguard
// from GATE's getroot_K1_Poi(). Returns root; sets converged.
static double getroot_K1_Poi(size_t n, const double mu[], const double g[],
	double q, bool &converged, double tol=1.490116e-08, int maxiter=1000)
{
	double t = 0;
	double K1_eval = K1_adj_Poi(t, n, mu, g, q);
	double prevJump = R_PosInf;
	converged = false;
	for (int rep=1; rep <= maxiter; rep++)
	{
		double K2_eval = K2_Poi(t, n, mu, g);
		double tnew = t - K1_eval / K2_eval;
		if (!R_FINITE(tnew)) { converged = false; break; }
		if (std::fabs(tnew - t) < tol) { t = tnew; converged = true; break; }
		double newK1 = K1_adj_Poi(tnew, n, mu, g, q);
		if (((K1_eval < 0) != (newK1 < 0)))  // sign change
		{
			if (std::fabs(tnew - t) > prevJump - tol)
			{
				tnew = t + ((newK1 > K1_eval) ? 1 : -1) * prevJump / 2;
				newK1 = K1_adj_Poi(tnew, n, mu, g, q);
				prevJump /= 2;
			} else {
				prevJump = std::fabs(tnew - t);
			}
		}
		t = tnew;
		K1_eval = newK1;
	}
	return t;
}

// Lugannani-Rice tail probability at saddlepoint zeta (signed, as in GATE)
static double Get_Saddle_Prob_Poi(double zeta, size_t n, const double mu[],
	const double g[], double q)
{
	double k1 = Korg_Poi(zeta, n, mu, g);
	double k2 = K2_Poi(zeta, n, mu, g);
	double pval;
	if (R_FINITE(k1) && R_FINITE(k2))
	{
		double temp1 = zeta * q - k1;
		double w = ((zeta > 0) ? 1.0 : -1.0) * std::sqrt(2.0 * temp1);
		double v = zeta * std::sqrt(k2);
		double Ztest = w + std::log(v / w) / w;
		if (Ztest > 0)
			pval = ::Rf_pnorm5(Ztest, 0, 1, FALSE, FALSE);  // upper tail
		else
			pval = - ::Rf_pnorm5(Ztest, 0, 1, TRUE, FALSE);  // -lower tail
	} else {
		pval = 0;
	}
	return pval;
}

// Two-sided saddlepoint p-value for the score statistic; falls back to
// pval_noadj if a root fails to converge. Returns p-value, sets converged.
double Saddle_Prob_Poisson(double Score, double pval_noadj, size_t n,
	const double mu[], const double g[], bool &converged)
{
	bool c1, c2;
	double root1 = getroot_K1_Poi(n, mu, g,  Score, c1);
	double root2 = getroot_K1_Poi(n, mu, g, -Score, c2);
	if (c1 && c2)
	{
		double p1 = Get_Saddle_Prob_Poi(root1, n, mu, g,  Score);
		double p2 = Get_Saddle_Prob_Poi(root2, n, mu, g, -Score);
		if (!R_FINITE(p1)) p1 = pval_noadj / 2;
		if (!R_FINITE(p2)) p2 = pval_noadj / 2;
		converged = true;
		return std::fabs(p1) + std::fabs(p2);
	} else {
		converged = false;
		return pval_noadj;
	}
}

// ===========================================================================
// (3) Cox risk-set covariance operator  W - U
//
// Builds the pieces needed by WminusU::apply_inv(): N=exp(eta), Winv=1/mu, the
// distinct event times (Breslow ties), the risk-set sums S_k, the index vector
// RvecIndex, and the m x m Woodbury capacitance inverse Ainv.

void WminusU::build(size_t n_, const double eta[], const double time[],
	const int status[], const double mu[])
{
	n = n_;
	N.set_size(n); Winv.set_size(n);
	for (size_t i=0; i < n; i++)
	{
		N[i] = std::exp(eta[i]);
		Winv[i] = 1.0 / mu[i];
	}
	std::vector<size_t> ord(n);
	for (size_t i=0; i < n; i++) ord[i] = i;
	std::stable_sort(ord.begin(), ord.end(),
		[&](size_t a, size_t b){ return time[a] < time[b]; });
	std::vector<double> Ssuf(n+1); Ssuf[n] = 0;
	for (size_t r=n; r-- > 0; ) Ssuf[r] = Ssuf[r+1] + N[ord[r]];
	std::vector<double> utime, dvec_d, Svec;
	size_t r = 0;
	while (r < n)
	{
		size_t r2 = r + 1;
		while (r2 < n && time[ord[r2]] == time[ord[r]]) r2++;
		int d = 0;
		for (size_t k=r; k < r2; k++) d += status[ord[k]];
		if (d > 0)
		{
			utime.push_back(time[ord[r]]);
			dvec_d.push_back((double)d);
			Svec.push_back(Ssuf[r]);
		}
		r = r2;
	}
	m = utime.size();
	RvecIndex.assign(n, 0);
	for (size_t i=0; i < n; i++)
	{
		size_t lo=0, hi=m;
		while (lo < hi) { size_t mid=(lo+hi)/2;
			if (utime[mid] <= time[i]) lo=mid+1; else hi=mid; }
		RvecIndex[i] = (int)lo;
	}
	arma::vec D(m);
	for (size_t k=0; k < m; k++) D[k] = dvec_d[k] / (Svec[k]*Svec[k]);
	arma::vec M(n);
	for (size_t i=0; i < n; i++) M[i] = N[i]*N[i]*Winv[i];
	arma::vec Tsuf = Rt_times(M);
	arma::mat ACm(m, m);
	for (size_t k=0; k < m; k++)
		for (size_t l=0; l < m; l++)
			ACm(k,l) = Tsuf[std::max(k,l)];
	for (size_t k=0; k < m; k++) ACm(k,k) -= 1.0 / D[k];
	Ainv = arma::pinv(arma::symmatu(ACm));
}


// ===========================================================================
// (4) Cox proportional-hazards model without random effects
//
// Independent samples: the same Cox-via-Poisson iteration as get_coeff() in
// saige_fitnull.cpp, but without the mixed-model (PCG) solver. Given eta, the
// Breslow hazard gives mu = Lambda0*exp(eta); the Poisson working response is
// Y = eta + (y-mu)/mu with weight W = mu, and the coefficients are updated by
// weighted least squares, alpha = (X'WX)^-1 X'WY. The fixed point of the
// iteration is the Cox MLE (Breslow ties). Since the Poisson information X'WX
// ignores the risk-set term U, the iteration converges linearly and is capped
// at 'maxiter' iterations (the caller reports non-convergence).
// Input:  y (0/1 event status), X (no intercept), time, alpha (initial)
// Output: alpha, n_iter, eta, mu, cov = (X'WX)^-1; returns convergence

bool fit_cox_noRE(const arma::vec &y, const arma::mat &X, const double time[],
	arma::vec &alpha, int maxiter, int &n_iter, arma::vec &eta, arma::vec &mu,
	arma::mat &cov)
{
	const size_t n = y.n_elem;
	const double tol_coef = 1e-4;   // as in get_coeff() for survival
	const double w_floor  = 1e-4;   // cf. surv_working_response()
	std::vector<int> status(n);
	for (size_t i=0; i < n; i++) status[i] = (int)std::lround(y[i]);

	arma::vec Lambda0(n), W(n), Y(n);
	eta = X * alpha;
	bool converged = false;
	n_iter = 0;
	for (int iter=1; iter <= maxiter; iter++)
	{
		n_iter = iter;
		breslow_lambda0(n, eta.memptr(), time, &status[0], Lambda0.memptr());
		mu = Lambda0 % arma::exp(eta);
		// subjects with mu ~ 0 (censored before the first event) carry no
		// information: guard the working response (0/0) and the weight
		for (size_t i=0; i < n; i++)
		{
			if (mu[i] > w_floor)
			{
				W[i] = mu[i];
				Y[i] = eta[i] + (y[i] - mu[i]) / mu[i];
			} else {
				W[i] = w_floor;
				Y[i] = eta[i];
			}
		}
		arma::mat XW = X.each_col() % W;                   // W X
		cov = arma::inv_sympd(arma::symmatu(X.t() * XW));  // (X' W X)^-1
		arma::vec a = cov * (XW.t() * Y);
		eta = X * a;
		double dlt = arma::max(arma::abs(a - alpha) /
			(arma::abs(a) + arma::abs(alpha) + tol_coef));
		alpha = a;
		if (dlt < tol_coef) { converged = true; break; }
	}
	// Poisson means at the final linear predictor
	breslow_lambda0(n, eta.memptr(), time, &status[0], Lambda0.memptr());
	mu = Lambda0 % arma::exp(eta);
	return converged;
}

}  // namespace SAIGE_SURV


// ===========================================================================
// R-callable functions

extern "C" {

// fit the Cox model without random effects (see SAIGE_SURV::fit_cox_noRE),
// called from .fit_survival() when neither 'gdsfile' nor 'grm.mat' is given;
// returns the same components as saige_fit_AI_PCG()
RcppExport SEXP saige_surv_fit_noRE(SEXP r_y, SEXP r_X, SEXP r_time,
	SEXP r_coef0, SEXP r_maxiter, SEXP r_verbose)
{
BEGIN_RCPP
	const arma::vec y = as<arma::vec>(r_y);
	const arma::mat X = as<arma::mat>(r_X);
	arma::vec alpha = as<arma::vec>(r_coef0);
	const int maxiter = Rf_asInteger(r_maxiter);
	const bool verbose = Rf_asLogical(r_verbose)==TRUE;
	const size_t n = y.n_elem;
	if (X.n_rows != n)
		throw std::invalid_argument("'X' does not match the sample size.");
	if (!Rf_isReal(r_time) || (size_t)Rf_length(r_time) != n)
		throw std::invalid_argument(
			"'eventTime' should be a numeric vector matching the sample size.");
	if (alpha.n_elem != X.n_cols)
		throw std::invalid_argument("'coef0' does not match the columns of 'X'.");
	if (maxiter < 1)
		throw std::invalid_argument("'maxiter' should be >= 1.");

	int n_iter = 0;
	arma::vec eta, mu;
	arma::mat cov;
	const bool converged = SAIGE_SURV::fit_cox_noRE(y, X, REAL(r_time),
		alpha, maxiter, n_iter, eta, mu, cov);
	if (verbose)
	{
		Rprintf("Cox proportional-hazards model (no random effects):\n");
		Rprintf("    # of iterations: %d\n", n_iter);
		Rprintf("    fixed coeff: (");
		for (size_t i=0; i < alpha.n_elem; i++)
			Rprintf(i ? ", %0.7g" : "%0.7g", alpha[i]);
		Rprintf(")\n");
	}
	arma::vec resid = y - mu;
	return List::create(
		_["coefficients"] = NumericVector(alpha.begin(), alpha.end()),
		_["tau"] = NumericVector::create(1.0, 0.0),
		_["linear.predictors"] = NumericVector(eta.begin(), eta.end()),
		_["fitted.values"] = NumericVector(mu.begin(), mu.end()),
		_["residuals"] = NumericVector(resid.begin(), resid.end()),
		_["cov"] = cov,
		_["converged"] = converged);
END_RCPP
}


// ===========================================================================
// R-callable test wrappers (used for validation against the GATE reference)

RcppExport SEXP saige_surv_lambda0(SEXP eta, SEXP time, SEXP status)
{
BEGIN_RCPP
	const size_t n = Rf_length(eta);
	NumericVector Lambda0(n);
	SAIGE_SURV::breslow_lambda0(n, REAL(eta), REAL(time), INTEGER(status),
		REAL(Lambda0));
	return Lambda0;
END_RCPP
}

RcppExport SEXP saige_surv_spa_poisson(SEXP score, SEXP pval_noadj, SEXP mu,
	SEXP g)
{
BEGIN_RCPP
	const size_t n = Rf_length(mu);
	bool conv = false;
	double p = SAIGE_SURV::Saddle_Prob_Poisson(Rf_asReal(score),
		Rf_asReal(pval_noadj), n, REAL(mu), REAL(g), conv);
	return List::create(_["p.value"]=p, _["converged"]=conv);
END_RCPP
}

// (W-U)^{-1} b for validation against a dense reference
RcppExport SEXP saige_surv_wminusU_inv(SEXP eta, SEXP time, SEXP status,
	SEXP mu, SEXP b)
{
BEGIN_RCPP
	const size_t n = Rf_length(eta);
	SAIGE_SURV::WminusU op;
	op.build(n, REAL(eta), REAL(time), INTEGER(status), REAL(mu));
	arma::vec bv((double*)REAL(b), n, false);
	arma::vec out = op.apply_inv(bv);
	return Rcpp::wrap(out);
END_RCPP
}

}  // extern "C"
