// ===========================================================================
//
// saige_surv.h: survival (time-to-event) support for SAIGEgds
//
// Shared declarations for the Cox-via-Poisson survival path.
// Copyright (C) 2026    Xiuwen Zheng;  License: GPL-3
//

#ifndef _SAIGE_SURV_H_
#define _SAIGE_SURV_H_

#include <RcppArmadillo.h>
#include <vector>
#include <algorithm>
#include <cmath>

namespace SAIGE_SURV
{

// Breslow baseline cumulative hazard (saige_surv.cpp)
void breslow_lambda0(size_t n, const double eta[], const double time[],
	const int status[], double Lambda0[]);

// weighted-Poisson saddlepoint p-value (saige_surv.cpp)
double Saddle_Prob_Poisson(double Score, double pval_noadj, size_t n,
	const double mu[], const double g[], bool &converged);


// Cox risk-set covariance operator  W - U,  with U = N R D R' N.
// (W - U)^{-1} applied via Woodbury; the hot operators are inline here,
// build() lives in saige_surv.cpp.
class WminusU
{
public:
	size_t n, m;
	std::vector<int> RvecIndex;   // n; 1..m (0 if no event time <= t_i)
	arma::vec N, Winv;            // n; N=exp(eta), Winv = 1/mu
	arma::mat Ainv;               // m x m

	// out[i] = sum_{k <= RvecIndex[i]} a[k]   (R %*% a)
	inline arma::vec R_times(const arma::vec &a) const
	{
		arma::vec pref(m+1); pref[0] = 0;
		for (size_t k=0; k < m; k++) pref[k+1] = pref[k] + a[k];
		arma::vec out(n);
		for (size_t i=0; i < n; i++) out[i] = pref[RvecIndex[i]];
		return out;
	}
	// out[k] = sum_{i: RvecIndex[i] >= k+1} v[i]   (t(R) %*% v), k=0..m-1
	inline arma::vec Rt_times(const arma::vec &v) const
	{
		arma::vec bucket(m+1, arma::fill::zeros);
		for (size_t i=0; i < n; i++) bucket[RvecIndex[i]] += v[i];
		arma::vec out(m);
		double suf = 0;
		for (size_t k=m; k >= 1; k--) { suf += bucket[k]; out[k-1] = suf; }
		return out;
	}

	// (W - U)^{-1} b
	inline arma::vec apply_inv(const arma::vec &b) const
	{
		arma::vec NWinvb = N % Winv % b;
		arma::vec rm = Rt_times(NWinvb);
		arma::vec am = Ainv * rm;
		arma::vec corr = N % Winv % R_times(am);
		return Winv % b - corr;
	}

	// build from eta, event time, status (1=event) and mu
	// (saige_surv.cpp)
	void build(size_t n_, const double eta[], const double time[],
		const int status[], const double mu[]);
};

}  // namespace SAIGE_SURV

#endif // _SAIGE_SURV_H_
