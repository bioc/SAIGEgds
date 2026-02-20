// ===========================================================
//
// saige_permu.cpp: Miscellaneous functions linked to the SKAT package
//

#include <RcppArmadillo.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cmath>
#include <Rmath.h>


// Compute the hypergeometric probability
typedef void (*Type_GetProb)(int &k, int &ngroup, int &ncase, int group[],
	double weight[], double prob[]);
// Compute the p-value from exact tests
typedef void (*Type_SKATExact)(int *resarray, int &nres, int *nres_k,
	double *Z0, double *Z1, int &k, int &m, int &total, int *total_k,
	double *prob_k, double *odds, double *p1, int *IsExact, double *pval,
	double *pval_same, double &minP, int &test_type, double &epsilon);

// Package name
static const char *PKG_SKAT = "SKAT";
// Whether has loaded the C functions in SKAT or not
static bool Have_SKAT_Func = false;

static Type_GetProb skat_GetProb = NULL;
static Type_SKATExact skat_SKATExact = NULL;

#define LOAD_SKAT(name, fc)    { \
	DL_FUNC f = R_FindSymbol(name, PKG_SKAT, NULL); \
	if (!f) Rf_error("Not found '%s' in the SKAT package!", name); \
	std::memcpy(&fc, &f, sizeof(f)); \
}

// Get the functions in the SKAT package
static void Init_SKAT_Pkg_Functions()
{
	if (!Have_SKAT_Func)
	{
		// load SKAT pacakge
		Rcpp::Environment base = Rcpp::Environment::base_env();
		Rcpp::Function requireNamespace = base["requireNamespace"];
		requireNamespace("SKAT");
		// load C functions
		LOAD_SKAT("RGetProb", skat_GetProb);
		LOAD_SKAT("RSKATExact", skat_SKATExact);
		Have_SKAT_Func = true;
	}
}

RcppExport SEXP saige_init_skat_pkg()
{
	Init_SKAT_Pkg_Functions();
	return R_NilValue;
}


// Compute the hypergeometric probability
static void GetHyperGeoProb(int k, int ngroup, int ncase, int group[],
	double weight[], double prob[])
{
	Init_SKAT_Pkg_Functions();
	// call RGetProb in the SKAT package
	(*skat_GetProb)(k, ngroup, ncase, group, weight, prob);
}

// Computes the binomial coefficient "x choose y"
static inline void Get_Total_K(int k, std::vector<int> &n_total_k)
{
	for(int i=0; i <= k; i++)
		n_total_k[i] = (int)Rf_choose(k, i);
}


void SKATExactBin_ComputeProb_Group(arma::uvec &idx, arma::uvec &idxCompVec,
	arma::vec &pi1, uint32_t n, uint32_t ncase, int type_group,
	std::vector<double> &prob)
{
	const int k = idx.n_elem;
	// use the default value as ER will only be used for variants with MAC <= 10
	int ngroup1 = 10;
	arma::vec p1 = pi1(idx);
	arma::vec p2 = pi1(idxCompVec);
	arma::uvec id_temp = arma::find(p1 >= 1);

	if (id_temp.n_elem > 0)
	{
		for(unsigned int j = 0; j < id_temp.n_elem; j++)
		{
			unsigned int id_temp_j = id_temp(j);
			p1(id_temp_j) = 0.999;
		}
	}

	std::vector<double> weight;
	std::vector<int> group;
	arma::uvec a1Vec, a2Vec, IDX;
	double a1, a2, p1temp, p2temp, oddtemp, p2oddtemp;
	arma::vec p1tempVec;

	for(int i=0; i < ngroup1; i++)
	{
		a1 = double(i) / ngroup1;
		a2 = double(i + 1) / ngroup1;
		if ((i+1) < ngroup1)
		{
			a1Vec = arma::find(p1 >= a1);				
			a2Vec = arma::find(p1 < a2);				
		} else {	
			a1Vec = arma::find(p1 >= a1);				
			a2Vec = arma::find(p1 <= a2);					
		}
		IDX = arma::intersect(a1Vec, a2Vec);
		
		if (IDX.n_elem > 0)
		{
			p1tempVec = p1(IDX);
			p1temp = arma::mean(p1tempVec);
			oddtemp = p1temp/(1-p1temp);
			weight.push_back(oddtemp);
			group.push_back(IDX.n_elem);
		}	
    }
	p2temp = arma::mean(p2);
	p2oddtemp = p2temp / (1-p2temp);
	weight.push_back(p2oddtemp);

	for(int i=0; i < weight.size(); i++)
		weight[i] = weight[i] / p2oddtemp;
	group.push_back(n-k);
	//std::vector<double> prob_k(k+1, 0.0);
	int ngroup = group.size();
	int ncasei = int(ncase);
	GetHyperGeoProb(k, ngroup, ncasei, &group[0], &weight[0], &prob[0]);
}


void SKATExactBin_ComputProb_New(arma::uvec &idx, arma::uvec &idxCompVec,
	arma::vec &pi1, uint32_t n, uint32_t ncase, int NResampling, int ExactMax,
	int test_type, int type_group, std::vector<double> &prob,
	std::vector<int> &IsExactVec, std::vector<int> &n_total_k, int &n_total,
	bool & Is_ExactP)
{
	const int k = idx.n_elem;
	SKATExactBin_ComputeProb_Group(idx, idxCompVec, pi1, n, ncase,
		type_group, prob);

	Get_Total_K(k, n_total_k);
	n_total = std::accumulate(n_total_k.begin(), n_total_k.end(), 0);
	Is_ExactP = true;
	if (n_total > NResampling)
	{
		for(int i=0; i <= k; i++)
		{
			if (n_total_k[i] > ExactMax)
			{
				n_total_k[i] = int(ceil(NResampling * prob[i]));
				IsExactVec[i] = 0;
			}	
		}
		Is_ExactP = false;
	}	
	n_total = std::accumulate(n_total_k.begin(), n_total_k.end(), 0);
}	


void Get_Res_Arrays(arma::mat &res_out, arma::uvec &idx,
	std::vector<int> &resarray, int &nres, std::vector<int> &nres_k)
{
	arma::vec res_out_colvec;
	for (int i=0; i < nres; i++)
	{
		res_out_colvec = res_out.col(i);
		arma::uvec res_out_i = arma::find( res_out_colvec > 0);
		arma::uvec res_out_i_s = arma::sort(res_out_i);
		int res_out_i_s_k = res_out_i_s.n_elem;
		nres_k[i] = res_out_i_s.n_elem;
		for(int k=0; k < res_out_i_s_k; k++)
			resarray.push_back(res_out_i_s(k));
	}
}


double SKATExactBin_Work(arma::mat &Z, arma::vec &res, arma::vec &pi1,
	uint32_t ncase, arma::uvec &idx, arma::uvec &idxCompVec,
	arma::mat &res_out, int NResampling, int ExactMax, double epsilon,
	int test_type)
{
	uint32_t n = res.n_elem;
	arma::vec p1 = pi1(idx);
	arma::vec p2 = pi1(idxCompVec);

	arma::mat Z_1 = Z.rows(idx);
	arma::mat Z1temp = (Z_1 % (-p1)).t();
	arma::vec Z0 = arma::vectorise(Z1temp);

	arma::mat Z1temp2 = (Z_1 % (1-p1)).t();
	arma::vec Z1 =  arma::vectorise(Z1temp2);

	int m = Z_1.n_cols;
	int k = idx.n_elem;
	std::vector<int> n_total_k(k+1, 0);

	std::vector<double> prob(k+1, 0.0);
	std::vector<int> IsExactVec(k+1, 1);

	int n_total = 0;
	bool Is_ExactP = false;
	SKATExactBin_ComputProb_New(idx, idxCompVec, pi1, n, ncase,
		NResampling, ExactMax, test_type, 2, prob, IsExactVec,
		n_total_k, n_total, Is_ExactP);

	double p1mean = arma::mean(p1);
	arma::vec p1_adj = p1 / p1mean;
	arma::vec odds = p1 / (1 - p1);
	int test_type_new=1;

	if(res_out.is_empty())
	{
		res_out = res(idx);
	} else {
		arma::vec res_out2 = arma::join_cols(res(idx), res_out(idx));
		res_out.resize(res_out2.n_elem);
		res_out = res_out2;
	}	

	std::vector<int> resarray;
	int nres = res_out.n_cols;
	std::vector<int> nres_k(nres, 0);
	
	Get_Res_Arrays(res_out, idx, resarray, nres, nres_k);
	std::vector<double> pval(nres, 0.0); 
	std::vector<double> pval1(nres, 0.0); 
	double minP = 100;
	
	typedef std::vector<double> stdvec;
	stdvec Z1std = arma::conv_to< stdvec >::from(Z1);
	stdvec Z0std = arma::conv_to< stdvec >::from(Z0);
	stdvec oddsstd = arma::conv_to< stdvec >::from(odds);
	stdvec p1_adjstd = arma::conv_to< stdvec >::from(p1_adj);

	(*skat_SKATExact)(&resarray[0], nres, &nres_k[0], &Z0std[0], &Z1std[0],
		k, m, n_total, &n_total_k[0], &prob[0], &oddsstd[0], &p1_adjstd[0],
		&IsExactVec[0], &pval[0], &pval1[0], minP, test_type_new, epsilon);

	// output p-value
	return pval[0] - pval1[0]/2;
}
