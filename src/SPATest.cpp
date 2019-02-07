// ===========================================================
//
// SPATest.cpp: C implementation of part of the R SPATest package
//
// Copyright (C) 2019    Xiuwen Zheng
//
// This file is part of SAIGEgds. It was created based on the R codes in the
// SPATest package.
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

#include "vectorization.h"
#include <Rdefines.h>
#include <R.h>
#include <Rmath.h>
#include <float.h>


inline static double sq(double v) { return v*v; }

// inline static int sign(double v) { return (v>0) ? 1 : ((v<0) ? -1 : 0); }

static double Korg(double t, size_t n_g, const double mu[], const double g[])
{
	double sum = 0;
	for (size_t i=0; i < n_g; i++)
	{
		double temp = log(1 - mu[i] + mu[i] * exp(g[i] * t));
		if (R_FINITE(temp)) sum += temp;
	}
	return sum;
}


static COREARRAY_TARGET_CLONES
	double K1_adj(double t, size_t n_g, const double mu[], const double g[], double q)
{
	double sum = 0;
	for (size_t i=0; i < n_g; i++)
	{
		double temp1 = (1 - mu[i])* exp(-g[i] * t) + mu[i];
		double temp2 = mu[i] * g[i];
		double temp = temp2 / temp1;
		if (R_FINITE(temp)) sum += temp;
	}
	return sum - q;
}


static COREARRAY_TARGET_CLONES
	double K2(double t, size_t n_g, const double mu[], const double g[])
{
	double sum = 0;
	for (size_t i=0; i < n_g; i++)
	{
		double temp1 = sq((1 - mu[i])* exp(-g[i]*t) + mu[i]);
		double temp2 = (1 - mu[i]) * mu[i] * sq(g[i]) * exp(-g[i]*t);
		double temp = temp2 / temp1;
		if (R_FINITE(temp)) sum += temp;
	}
	return sum;
}


// .Machine$double.eps^0.25
static const double root_tol = sqrt(sqrt(DBL_EPSILON));

static void COREARRAY_TARGET_CLONES
	getroot_K1(double &root, int &n_iter, bool &converged,
		double init, size_t n_g, const double mu[],
		const double g[], double q, double tol=root_tol, int maxiter=1000)
{
	double g_pos = 0;
	for (size_t i=0; i < n_g; i++) g_pos += (g[i] > 0) ? g[i] : 0;
	double g_neg = 0;
	for (size_t i=0; i < n_g; i++) g_neg += (g[i] < 0) ? g[i] : 0;

	if (q>=g_pos || q<=g_neg)
	{
		root = R_PosInf; n_iter = 0;
		converged = true;

	} else {
		double t = root = init;
		double K1_eval = K1_adj(t, n_g, mu, g, q);
		double prevJump = R_PosInf;
		converged = false;
		for (int n_iter=1; n_iter <= maxiter; n_iter++)
		{
			double K2_eval = K2(t, n_g, mu, g);
			double tnew = t - K1_eval/K2_eval;
			if (!R_FINITE(tnew)) break;
			if (abs(tnew - t) < tol)
			{
				converged = true;
				break;
			}
			double newK1 = K1_adj(tnew, n_g, mu, g, q);
			if (sign(K1_eval) != sign(newK1))
			{
				if (abs(tnew - t) > prevJump-tol)
				{
					tnew = t + sign(newK1 - K1_eval)*prevJump/2;
					newK1 = K1_adj(tnew, n_g, mu, g, q);
					prevJump = prevJump / 2;
				} else {
					prevJump = abs(tnew - t);
				}
			}
			root = t = tnew;
			K1_eval = newK1;
		}
	}
}


static double Get_Saddle_Prob(double zeta, size_t n_g, const double mu[],
	const double g[], double q)
{
	double k1 = Korg(zeta, n_g, mu, g);
	double k2 = K2(zeta, n_g, mu, g);
	double pval = 0;
	if (R_FINITE(k1) && R_FINITE(k1))
	{
		double temp1 = zeta * q - k1;
		double w = sign(zeta) * sqrt(2 * temp1);
		double v = zeta * sqrt(k2);
		double Z_test = w + 1/w * log(v/w);
		if (Z_test > 0)
			pval = ::Rf_pnorm5(Z_test, 0, 1, FALSE, FALSE);
		else
			pval = - ::Rf_pnorm5(Z_test, 0, 1, TRUE, FALSE);
	}
	return pval;
}




// m1<-sum(mu * g),  var1<-sum(mu * (1-mu) * g^2)
extern "C" double Saddle_Prob(double q, double m1, double var1, size_t n_g,
	const double mu[], const double g[], double cutoff, bool &converged)
{
	double s = q - m1;
	double qinv = -s + m1;
	double pval_noadj = Rf_pchisq(s*s/var1, 1, FALSE, FALSE);
	double pval;

	while (true)
	{
		converged = true;
		if (cutoff < 0.1) cutoff = 0.1;

		if (abs(q - m1)/sqrt(var1) < cutoff)
		{
			pval = pval_noadj;
		} else {
			double uni1_root, uni2_root;
			int n_iter1, n_iter2;
			bool conv1, conv2;
			getroot_K1(uni1_root, n_iter1, conv1, 0, n_g, mu, g, q);
			getroot_K1(uni2_root, n_iter2, conv2, 0, n_g, mu, g, qinv);
			if (conv1 && conv2)
			{
				double p1 = Get_Saddle_Prob(uni1_root, n_g, mu, g, q);
				double p2 = Get_Saddle_Prob(uni2_root, n_g, mu, g, qinv);
				pval = abs(p1) + abs(p2);
			} else {
				pval = pval_noadj;
				converged = false;
			}				
		}

		if (pval!=0 && pval_noadj/pval>1000)
			cutoff *= 2;
		else
			break;
	}
	return pval;
}
