// ===========================================================
//
// saige.h: the header file for the SAIGEgds package
//
// Copyright (C) 2022    Xiuwen Zheng / AbbVie-ComputationalGenomics
//
// This file is part of SAIGEgds. It was created based on the original SAIGE
// C++ and R codes in the SAIGE package. Compared with the original SAIGE,
// all single-precision floating-point numbers are changed to double precision,
// and a more efficient algorithm with packed 2-bit and sparse genotypes is
// implemented to calculate the cross product of genetic relationship matrix.
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
// with SAIGEgds.
// If not, see <http://www.gnu.org/licenses/>.

#include <Rdefines.h>

// to avoid the conflict with C++
#ifdef length
#   undef length
#endif

#ifdef isNull
#   undef isNull
#endif

#include <RcppParallel.h>
#include <tbb/parallel_for.h>



// ========================================================================= //

/// define BYTE type
#define BYTE    unsigned char


// ========================================================================= //
// Define Intel TBB macro with a thread index starting from 0
// requiring C++11

#if defined(RCPP_PARALLEL_USE_TBB) && RCPP_PARALLEL_USE_TBB

#define PARALLEL_HEAD(SIZE, balancing)    \
	tbb::parallel_for(  \
		balancing ? tbb::blocked_range<size_t>(0, SIZE) :  \
		tbb::blocked_range<size_t>(0, SIZE, SIZE/SAIGE_NumThread + (SIZE % SAIGE_NumThread ? 1:0)),  \
		[&](const tbb::blocked_range<size_t> &r)  \
	{  \
		const int th_idx = tbb::this_task_arena::current_thread_index();  \
		if (th_idx < 0 || th_idx >= SAIGE_NumThread)  \
			throw std::invalid_argument( \
				"Invalid tbb::this_task_arena::current_thread_index()!");

#define PARALLEL_FOR(i, SIZE, balancing)    \
		PARALLEL_HEAD(SIZE, balancing)    \
		for (size_t i=r.begin(); i < r.end(); i++)

#define PARALLEL_RANGE(st, ed, SIZE, balancing)    \
		PARALLEL_HEAD(SIZE, balancing)    \
		const size_t st = r.begin(), ed = r.end();

#define PARALLEL_END    });

#define PARALLEL_THREAD_BLOCK    \
		tbb::task_arena arena(SAIGE_NumThread); \
		arena.execute([&]{
#define PARALLEL_THREAD_BLOCK_END    });

#else

#define PARALLEL_FOR(i, SIZE, balancing)    \
	{  \
		const int th_idx = 0;  \
		for (size_t i=0; i < (size_t)SIZE; i++)
#define PARALLEL_RANGE(st, ed, SIZE, balancing)    \
	{  \
		const int th_idx = 0;  \
		const size_t st = 0, ed = SIZE;
#define PARALLEL_END    }

#define PARALLEL_THREAD_BLOCK
#define PARALLEL_THREAD_BLOCK_END

#endif


// ========================================================================= //

namespace SAIGE
{
	/// Type of outcome variables
	enum TTrait
	{
		Unknown = 0,  //< used for the initial value
		Quant   = 1,  //< continuous outcomes
		Binary  = 2   //< binary outcomes, e.g., case-control study
	};

	/// the number of threads used in the SAIGEgds package
	extern int SAIGE_NumThread;


	// =================================================================

	/// const wrapper for dense matrix
	class Type_Matrix
	{
	public:
		const double *val;  //< pointer to the data

		// constructor
		Type_Matrix(SEXP mat=NULL) { reset(mat); }
		// destructor
		~Type_Matrix() { reset(NULL); }

		void reset(SEXP mat);
		inline int nrow() const { return m_nrow; }
		inline int ncol() const { return m_ncol; }
		inline bool empty() const { return !val; }
	protected:
		int m_nrow, m_ncol;
	};


	/// const wrapper for Matrix::dgCMatrix
	class Type_dgCMatrix
	{
	public:
		const int *i;     //< Matrix::dgCMatrix@i
		const int *p;     //< Matrix::dgCMatrix@p
		const double *x;  //< Matrix::dgCMatrix@x

		// constructor
		Type_dgCMatrix(SEXP mat=NULL) { reset(mat); }
		// destructor
		~Type_dgCMatrix() { reset(NULL); }

		void reset(SEXP mat);
		inline int nrow() const { return m_nrow; }
		inline int ncol() const { return m_ncol; }
		inline bool empty() const { return !i || !p || !x; }
	protected:
		int m_nrow, m_ncol;
	};


	/// const wrapper for Matrix::dsCMatrix
	class Type_dsCMatrix: public Type_dgCMatrix
	{
	public:
		// constructor
		Type_dsCMatrix(SEXP mat=NULL) { reset(mat); }
		// destructor
		~Type_dsCMatrix() { reset(NULL); }

		void reset(SEXP mat);
	};


	// =================================================================

	/// square
	static inline double sq(double v) { return v*v; }

}


// ========================================================================= //
/// SPAtest

extern "C" double Saddle_Prob(double q, double m1, double var1, size_t n_g,
	const double mu[], const double g[], double cutoff, bool &converged, double *p_noadj);

extern "C" double Saddle_Prob_Fast(double q, double m1, double var1, size_t n_g,
	const double mu[], const double g[], size_t n_nonzero, const int nonzero_idx[],
	double cutoff, bool &converged, double buf_spa[], double *p_noadj);

