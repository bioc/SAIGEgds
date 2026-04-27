// ===========================================================
//
// gpu_opencl.h: Optional OpenCL GPU acceleration
//
// Copyright (C) 2026    Xiuwen Zheng / AbbVie-ComputationalGenomics
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
// with SAIGEgds.
// If not, see <http://www.gnu.org/licenses/>.

#ifndef GPU_OPENCL_H
#define GPU_OPENCL_H

#include <cstddef>
#include <cstdint>
#include <vector>

#include <Rdefines.h>

// to avoid conflicts with R macros
#ifdef length
#   undef length
#endif

// ========================================================================= //
// Public GPU interface

/// Initialize the GPU context. Returns true if a usable GPU was found.
/// If OpenCL library is not available or no GPU device exists, returns false.
bool gpu_opencl_init(bool verbose);

/// Release all GPU resources (context, queue, buffers, program).
void gpu_opencl_cleanup();

/// Query whether GPU acceleration is currently available and initialized.
bool gpu_opencl_available();

/// Upload dense 2-bit packed genotypes to GPU device memory.
/// packed_geno: pointer to 2-bit packed genotype bytes (packed_nsamp * nvar)
/// std_geno: pointer to 4-by-nvar lookup table of standardized genotype values
/// nsamp, nvar, packed_nsamp: dimensions
/// Returns true on success.
bool gpu_upload_2b_geno(const uint8_t *packed_geno, const double *std_geno,
	int nsamp, int nvar, int packed_nsamp);

/// Upload a user-defined dense GRM matrix to GPU device memory.
/// mat: column-major nsamp x nsamp dense matrix
/// nsamp: matrix dimension
/// Returns true on success.
bool gpu_upload_dense_grm(const double *mat, int nsamp);

/// Perform GPU-accelerated GRM cross-product: out_b = (G'G / nvar) * b
/// b: input vector (nsamp-length)
/// out_b: output vector (nsamp-length), caller must allocate
/// Returns true if GPU computation succeeded; false to fall back to CPU.
bool gpu_grm_crossprod_2b(const double *b, double *out_b, int nsamp,
	int nvar, int packed_nsamp);

/// Perform GPU-accelerated dense GRM matrix-vector product: out_b = mat * b
/// b: input vector (nsamp-length)
/// out_b: output vector (nsamp-length), caller must allocate
/// Returns true if GPU computation succeeded; false to fall back to CPU.
bool gpu_grm_crossprod_dense(const double *b, double *out_b, int nsamp);

/// Free GPU genotype buffers (called when GRM storage is reset).
void gpu_free_buffers();


// ========================================================================= //
// GRM construction API

/// Progress callback type: called after each block with user-supplied data
typedef void (*gpu_progress_fn)(void *data);

/// Compute dense GRM on GPU from 2-bit packed genotypes.
/// g_pack: [n_byte x nsamp] raw matrix (col-major, sample-major packing).
/// grm_out: caller-allocated nsamp x nsamp double array (col-major).
/// Returns true on success; false to fall back to CPU.
bool gpu_grm_dense_calc(const uint8_t *g_pack, int nsamp, int nvar,
	int n_byte, bool use_f64, int bs, double *grm_out,
	gpu_progress_fn prog_cb = nullptr, void *prog_data = nullptr);

/// Sparse GRM step 1: scan all sample pairs in blocks, return (i,j) pairs
/// with GRM value >= rel_cutoff.
/// g_pack: [n_byte x nsamp], out_i/out_j: output 0-based sample index pairs.
/// Returns true on success.
bool gpu_grm_sparse_scan(const uint8_t *g_pack, int nsamp, int nvar,
	int n_byte, double rel_cutoff, int bs, std::vector<int> &out_i,
	std::vector<int> &out_j,
	gpu_progress_fn prog_cb = nullptr, void *prog_data = nullptr);

/// Sparse GRM step 2: refine GRM values for specific (i,j) pairs using full
/// variant set on GPU.
/// pair_i/pair_j: 0-based sample indices, npairs: count.
/// out_x: caller-allocated output array of length npairs.
/// Returns true on success.
bool gpu_grm_sparse_refine(const uint8_t *g_pack, int nsamp, int nvar,
	int n_byte, const int *pair_i, const int *pair_j, int npairs,
	double *out_x);


// R-callable entry points
extern "C" {
	SEXP saige_gpu_init(SEXP verbose);
	SEXP saige_gpu_cleanup();
}

#endif  // GPU_OPENCL_H
