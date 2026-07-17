// ===========================================================
//
// gpu_opencl.cpp: Optional OpenCL GPU acceleration (runtime-loaded)
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

#include "gpu_opencl.h"

#include <Rdefines.h>
#include <R_ext/Print.h>

#include <cmath>
#include <cstring>
#include <cstdlib>
#include <string>
#include <vector>
#include <algorithm>


// ========================================================================= //
// Platform-specific dynamic loading

#ifdef _WIN32
#   define WIN32_LEAN_AND_MEAN
#   include <windows.h>
    typedef HMODULE lib_handle_t;
#   define LIB_OPEN(name)    LoadLibraryA(name)
#   define LIB_SYM(h, name)  ((void*)GetProcAddress(h, name))
#   define LIB_CLOSE(h)      FreeLibrary(h)
    static const char *opencl_lib_names[] = { "OpenCL.dll", NULL };
#else
#   include <dlfcn.h>
    typedef void* lib_handle_t;
#   define LIB_OPEN(name)    dlopen(name, RTLD_LAZY)
#   define LIB_SYM(h, name)  dlsym(h, name)
#   define LIB_CLOSE(h)      dlclose(h)
#   ifdef __APPLE__
    static const char *opencl_lib_names[] = {
        "/System/Library/Frameworks/OpenCL.framework/OpenCL",
        "libOpenCL.dylib", NULL };
#   else
    static const char *opencl_lib_names[] = {
        "libOpenCL.so.1", "libOpenCL.so", NULL };
#   endif
#endif


// ========================================================================= //
// Minimal OpenCL type definitions (so we don't need CL headers at compile time)

typedef int32_t  cl_int;
typedef uint32_t cl_uint;
typedef uint64_t cl_ulong;
typedef cl_uint  cl_bool;
typedef cl_ulong cl_bitfield;
typedef cl_bitfield cl_device_type;
typedef cl_uint  cl_platform_info;
typedef cl_uint  cl_device_info;
typedef cl_uint  cl_context_info;
typedef cl_uint  cl_command_queue_properties;  // for OpenCL 1.x
typedef cl_bitfield cl_mem_flags;
typedef cl_uint  cl_program_build_info;

typedef struct _cl_platform_id *    cl_platform_id;
typedef struct _cl_device_id *      cl_device_id;
typedef struct _cl_context *        cl_context;
typedef struct _cl_command_queue *  cl_command_queue;
typedef struct _cl_mem *            cl_mem;
typedef struct _cl_program *        cl_program;
typedef struct _cl_kernel *         cl_kernel;
typedef struct _cl_event *          cl_event;

#define CL_SUCCESS                  0
#define CL_DEVICE_TYPE_GPU          (1 << 2)
#define CL_DEVICE_TYPE_ALL          0xFFFFFFFF
#define CL_PLATFORM_NAME            0x0902
#define CL_DEVICE_NAME              0x102B
#define CL_DEVICE_GLOBAL_MEM_SIZE   0x101F
#define CL_DEVICE_MAX_WORK_GROUP_SIZE  0x1004
#define CL_DEVICE_MAX_COMPUTE_UNITS    0x1002
#define CL_MEM_READ_ONLY            (1 << 2)
#define CL_MEM_WRITE_ONLY           (1 << 1)
#define CL_MEM_READ_WRITE           (1 << 0)
#define CL_MEM_COPY_HOST_PTR        (1 << 5)
#define CL_PROGRAM_BUILD_LOG        0x1183
#define CL_TRUE                     1
#define CL_FALSE                    0


// ========================================================================= //
// OpenCL function pointer types

typedef cl_int (*pfn_clGetPlatformIDs)(cl_uint, cl_platform_id *, cl_uint *);
typedef cl_int (*pfn_clGetPlatformInfo)(cl_platform_id, cl_platform_info, size_t, void *, size_t *);
typedef cl_int (*pfn_clGetDeviceIDs)(cl_platform_id, cl_device_type, cl_uint, cl_device_id *, cl_uint *);
typedef cl_int (*pfn_clGetDeviceInfo)(cl_device_id, cl_device_info, size_t, void *, size_t *);
typedef cl_context (*pfn_clCreateContext)(const void *, cl_uint, const cl_device_id *, void *, void *, cl_int *);
typedef cl_command_queue (*pfn_clCreateCommandQueue)(cl_context, cl_device_id, cl_command_queue_properties, cl_int *);
typedef cl_mem (*pfn_clCreateBuffer)(cl_context, cl_mem_flags, size_t, void *, cl_int *);
typedef cl_program (*pfn_clCreateProgramWithSource)(cl_context, cl_uint, const char **, const size_t *, cl_int *);
typedef cl_int (*pfn_clBuildProgram)(cl_program, cl_uint, const cl_device_id *, const char *, void *, void *);
typedef cl_int (*pfn_clGetProgramBuildInfo)(cl_program, cl_device_id, cl_program_build_info, size_t, void *, size_t *);
typedef cl_kernel (*pfn_clCreateKernel)(cl_program, const char *, cl_int *);
typedef cl_int (*pfn_clSetKernelArg)(cl_kernel, cl_uint, size_t, const void *);
typedef cl_int (*pfn_clEnqueueNDRangeKernel)(cl_command_queue, cl_kernel, cl_uint, const size_t *, const size_t *, const size_t *, cl_uint, const cl_event *, cl_event *);
typedef cl_int (*pfn_clEnqueueWriteBuffer)(cl_command_queue, cl_mem, cl_bool, size_t, size_t, const void *, cl_uint, const cl_event *, cl_event *);
typedef cl_int (*pfn_clEnqueueReadBuffer)(cl_command_queue, cl_mem, cl_bool, size_t, size_t, void *, cl_uint, const cl_event *, cl_event *);
typedef cl_int (*pfn_clFinish)(cl_command_queue);
typedef cl_int (*pfn_clReleaseMemObject)(cl_mem);
typedef cl_int (*pfn_clReleaseKernel)(cl_kernel);
typedef cl_int (*pfn_clReleaseProgram)(cl_program);
typedef cl_int (*pfn_clReleaseCommandQueue)(cl_command_queue);
typedef cl_int (*pfn_clReleaseContext)(cl_context);
typedef cl_int (*pfn_clEnqueueFillBuffer)(cl_command_queue, cl_mem, const void *, size_t, size_t, size_t, cl_uint, const cl_event *, cl_event *);


// ========================================================================= //
// Resolved function pointers (all initially NULL)

static struct {
	pfn_clGetPlatformIDs           clGetPlatformIDs;
	pfn_clGetPlatformInfo          clGetPlatformInfo;
	pfn_clGetDeviceIDs             clGetDeviceIDs;
	pfn_clGetDeviceInfo            clGetDeviceInfo;
	pfn_clCreateContext            clCreateContext;
	pfn_clCreateCommandQueue       clCreateCommandQueue;
	pfn_clCreateBuffer             clCreateBuffer;
	pfn_clCreateProgramWithSource  clCreateProgramWithSource;
	pfn_clBuildProgram             clBuildProgram;
	pfn_clGetProgramBuildInfo      clGetProgramBuildInfo;
	pfn_clCreateKernel             clCreateKernel;
	pfn_clSetKernelArg             clSetKernelArg;
	pfn_clEnqueueNDRangeKernel     clEnqueueNDRangeKernel;
	pfn_clEnqueueWriteBuffer       clEnqueueWriteBuffer;
	pfn_clEnqueueReadBuffer        clEnqueueReadBuffer;
	pfn_clFinish                   clFinish;
	pfn_clReleaseMemObject         clReleaseMemObject;
	pfn_clReleaseKernel            clReleaseKernel;
	pfn_clReleaseProgram           clReleaseProgram;
	pfn_clReleaseCommandQueue      clReleaseCommandQueue;
	pfn_clReleaseContext           clReleaseContext;
	pfn_clEnqueueFillBuffer        clEnqueueFillBuffer;
} CL = {};


// ========================================================================= //
// GPU state

static lib_handle_t cl_lib  = NULL;
static bool gpu_initialized = false;

static cl_context        gpu_context  = NULL;
static cl_command_queue  gpu_queue    = NULL;
static cl_device_id      gpu_device   = NULL;
static cl_program        gpu_program  = NULL;
static size_t            gpu_max_wg   = 0;
static size_t            gpu_local_wg = 0;  // preferred local work-group size

// Kernels
static cl_kernel kern_dot_2b     = NULL;  // pass 1: dot product per variant
static cl_kernel kern_accum_2b   = NULL;  // pass 2: accumulate per sample
static cl_kernel kern_gemv       = NULL;  // dense GRM GEMV
static cl_kernel kern_grm_pair   = NULL;  // GRM construction: pairwise block
static cl_kernel kern_grm_ijx    = NULL;  // GRM construction: refine (i,j) list

// Device buffers for 2-bit packed genotypes
static cl_mem buf_d_packed_geno  = NULL;  // packed genotype bytes
static cl_mem buf_d_std_geno     = NULL;  // 4-by-nvar standardized geno lookup
static cl_mem buf_d_b            = NULL;  // input vector b
static cl_mem buf_d_out          = NULL;  // output vector
static cl_mem buf_d_dots         = NULL;  // per-variant dot products

// Device buffers for dense GRM
static cl_mem buf_d_grm_mat      = NULL;  // dense GRM matrix

// Stored dimensions
static int gpu_nsamp = 0, gpu_nvar = 0, gpu_packed_nsamp = 0;
static bool gpu_has_2b = false, gpu_has_dense_grm = false;

// Persistent host-side float buffers (avoids per-call allocation + zero-init)
static float *host_b_f   = NULL;  // for double->float b conversion
static float *host_out_f = NULL;  // for float->double output conversion
static int    host_f_len = 0;     // allocated length of host_b_f / host_out_f


// ========================================================================= //
// OpenCL kernel source (embedded)

static const char *kernel_source = R"CL(

// Pass 1: compute dot product of standardized genotype column v with vector b
// Each work item handles one variant.
// packed_geno: 2-bit packed genotypes [packed_nsamp x nvar], column-major
// std_geno: lookup table [4 x nvar] for standardized genotype values
// b: input vector [nsamp]
// dots: output dot products [nvar]
__kernel void dot_2b(
	__global const uchar *packed_geno,
	__global const float *std_geno,
	__global const float *b,
	__global float *dots,
	const int nsamp,
	const int packed_nsamp,
	const int nvar)
{
	const int v = get_global_id(0);  // variant index
	if (v >= nvar) return;  // guard for rounded-up global size
	__global const uchar *g = packed_geno + (long)packed_nsamp * v;
	__global const float *base = std_geno + (v << 2);
	const float lut[4] = { base[0], base[1], base[2], base[3] };

	float dot = 0.0f;
	int j = 0;
	// Branch-free loop over full bytes (4 samples each)
	const int full_bytes = nsamp >> 2;
	for (int p = 0; p < full_bytes; p++)
	{
		uchar gg = g[p];
		dot += lut[gg & 0x03]        * b[j];
		dot += lut[(gg >> 2) & 0x03] * b[j+1];
		dot += lut[(gg >> 4) & 0x03] * b[j+2];
		dot += lut[gg >> 6]          * b[j+3];
		j += 4;
	}
	// Remaining 0-3 samples
	if (j < nsamp)
	{
		uchar gg = g[full_bytes];
		dot += lut[gg & 0x03] * b[j]; j++;
		if (j < nsamp) { dot += lut[(gg >> 2) & 0x03] * b[j]; j++; }
		if (j < nsamp) { dot += lut[(gg >> 4) & 0x03] * b[j]; }
	}
	dots[v] = dot;
}

// Pass 2: accumulate dot .* std_geno into output vector (one work item per sample)
// For each sample j, iterate over all variants and accumulate:
//   out[j] += dots[v] * std_geno[geno(j,v)] / nvar
__kernel void accum_2b(
	__global const uchar *packed_geno,
	__global const float *std_geno,
	__global const float *dots,
	__global float *out,
	const int nsamp,
	const int nvar,
	const int packed_nsamp,
	const float inv_nvar)
{
	const int j = get_global_id(0);  // sample index
	if (j >= nsamp) return;

	// position within the packed byte: sample j is at byte (j/4), shift (j%4)*2
	const int byte_off = j >> 2;
	const int bit_shift = (j & 3) << 1;

	float acc = 0.0f;
	for (int v = 0; v < nvar; v++)
	{
		uchar gg = packed_geno[(long)packed_nsamp * v + byte_off];
		int c = (gg >> bit_shift) & 0x03;
		acc += dots[v] * std_geno[(v << 2) + c];
	}
	out[j] = acc * inv_nvar;
}

// Dense GRM matrix-vector product: out = mat * b (float4 vectorized)
// One work item per row (sample)
__kernel void gemv_dense(
	__global const float *mat,
	__global const float *b,
	__global float *out,
	const int nsamp)
{
	const int i = get_global_id(0);
	if (i >= nsamp) return;

	__global const float *row = mat + (long)nsamp * i;
	float4 sum4 = (float4)(0.0f);
	const int nvec = nsamp >> 2;
	for (int k = 0; k < nvec; k++)
		sum4 += vload4(k, row) * vload4(k, b);
	float sum = sum4.x + sum4.y + sum4.z + sum4.w;
	for (int j = nvec << 2; j < nsamp; j++)
		sum += row[j] * b[j];
	out[i] = sum;
}


// ---- GRM construction kernels ----

// Compute GRM(i,j) for a block of sample pairs from 2-bit packed genotypes.
// std_geno_lut layout: 8 floats per group of 4 variants (byte), indexed by
//   combined genotype pair code. For each byte position b:
//   std_geno_lut[b*32 + code] = product of normalized genotypes.
// However, the R-side lookup is structured differently (8 x nrow*4):
//   For each byte b and 4 SNPs within, there are 8 entries encoding all
//   pair products. We use a simpler approach:
// std_geno: float[4 * n_byte * 4], for each variant v, std_geno[v*4 + c]
//   gives the standardized genotype for code c (0,1,2,3=missing->0).
// grm_out: float[block_i_n * block_j_n], row-major within the block.
// Each work item computes one (i,j) pair.
__kernel void grm_pairwise_2b(
	__global const uchar *packed_geno,
	__global const float *std_geno,
	__global float *grm_out,
	const int nsamp, const int nvar, const int n_byte,
	const int block_i_st, const int block_i_n,
	const int block_j_st, const int block_j_n)
{
	const int li = get_global_id(0);  // local index within block_i
	const int lj = get_global_id(1);  // local index within block_j
	if (li >= block_i_n || lj >= block_j_n) return;

	const int gi = block_i_st + li;  // global sample i
	const int gj = block_j_st + lj;  // global sample j
	if (gj < gi) return;  // upper triangle only

	// Precompute column base pointers (avoids per-variant address calc)
	__global const uchar *gi_col = packed_geno + (long)gi * n_byte;
	__global const uchar *gj_col = packed_geno + (long)gj * n_byte;
	float sum = 0.0f;
	int valid_n = 0;

	// Process 4 variants per byte (unrolled), one byte load per sample
	const int n_full_bytes = nvar >> 2;
	for (int bp = 0; bp < n_full_bytes; bp++)
	{
		uchar gi_byte = gi_col[bp];
		uchar gj_byte = gj_col[bp];
		// std_geno base for the 4 variants in this byte: bp*4 variants * 4 entries
		__global const float *sg = std_geno + (bp << 4);
		int ci, cj;
		// Variant 0 (bits 0-1)
		ci = gi_byte & 0x03;  cj = gj_byte & 0x03;
		sum += sg[ci] * sg[cj];  // std_geno[3]=0, so missing contributes 0
		valid_n += (ci < 3) & (cj < 3);
		// Variant 1 (bits 2-3)
		ci = (gi_byte >> 2) & 0x03;  cj = (gj_byte >> 2) & 0x03;
		sum += sg[4 + ci] * sg[4 + cj];
		valid_n += (ci < 3) & (cj < 3);
		// Variant 2 (bits 4-5)
		ci = (gi_byte >> 4) & 0x03;  cj = (gj_byte >> 4) & 0x03;
		sum += sg[8 + ci] * sg[8 + cj];
		valid_n += (ci < 3) & (cj < 3);
		// Variant 3 (bits 6-7)
		ci = gi_byte >> 6;  cj = gj_byte >> 6;
		sum += sg[12 + ci] * sg[12 + cj];
		valid_n += (ci < 3) & (cj < 3);
	}
	// Remaining 0-3 variants
	const int rem = nvar & 3;
	if (rem > 0)
	{
		uchar gi_byte = gi_col[n_full_bytes];
		uchar gj_byte = gj_col[n_full_bytes];
		__global const float *sg = std_geno + (n_full_bytes << 4);
		for (int r = 0; r < rem; r++)
		{
			int ci = (gi_byte >> (r << 1)) & 0x03;
			int cj = (gj_byte >> (r << 1)) & 0x03;
			sum += sg[(r << 2) + ci] * sg[(r << 2) + cj];
			valid_n += (ci < 3) & (cj < 3);
		}
	}

	float val = (valid_n > 0) ? (sum / valid_n) : 0.0f;
	grm_out[li * block_j_n + lj] = val;
}


// Compute GRM values for an explicit list of (i,j) pairs.
// pair_i, pair_j: arrays of sample indices (0-based)
// grm_x: output array of GRM values
// Each work item handles one pair.
__kernel void grm_refine_ijx(
	__global const uchar *packed_geno,
	__global const float *std_geno,
	__global const int *pair_i, __global const int *pair_j,
	__global float *grm_x,
	const int nvar, const int n_byte, const int npairs)
{
	const int k = get_global_id(0);
	if (k >= npairs) return;

	const int gi = pair_i[k];
	const int gj = pair_j[k];

	// Precompute column base pointers
	__global const uchar *gi_col = packed_geno + (long)gi * n_byte;
	__global const uchar *gj_col = packed_geno + (long)gj * n_byte;
	float sum = 0.0f;
	int valid_n = 0;

	// Process 4 variants per byte (unrolled)
	const int n_full_bytes = nvar >> 2;
	for (int bp = 0; bp < n_full_bytes; bp++)
	{
		uchar gi_byte = gi_col[bp];
		uchar gj_byte = gj_col[bp];
		__global const float *sg = std_geno + (bp << 4);
		int ci, cj;
		ci = gi_byte & 0x03;  cj = gj_byte & 0x03;
		sum += sg[ci] * sg[cj];  // std_geno[3]=0, so missing contributes 0
		valid_n += (ci < 3) & (cj < 3);
		ci = (gi_byte >> 2) & 0x03;  cj = (gj_byte >> 2) & 0x03;
		sum += sg[4 + ci] * sg[4 + cj];
		valid_n += (ci < 3) & (cj < 3);
		ci = (gi_byte >> 4) & 0x03;  cj = (gj_byte >> 4) & 0x03;
		sum += sg[8 + ci] * sg[8 + cj];
		valid_n += (ci < 3) & (cj < 3);
		ci = gi_byte >> 6;  cj = gj_byte >> 6;
		sum += sg[12 + ci] * sg[12 + cj];
		valid_n += (ci < 3) & (cj < 3);
	}
	const int rem = nvar & 3;
	if (rem > 0)
	{
		uchar gi_byte = gi_col[n_full_bytes];
		uchar gj_byte = gj_col[n_full_bytes];
		__global const float *sg = std_geno + (n_full_bytes << 4);
		for (int r = 0; r < rem; r++)
		{
			int ci = (gi_byte >> (r << 1)) & 0x03;
			int cj = (gj_byte >> (r << 1)) & 0x03;
			sum += sg[(r << 2) + ci] * sg[(r << 2) + cj];
			valid_n += (ci < 3) & (cj < 3);
		}
	}

	grm_x[k] = (valid_n > 0) ? (sum / valid_n) : 0.0f;
}

)CL";


// ========================================================================= //
// Internal helpers

/// Try to load the OpenCL shared library at runtime
static bool load_opencl_library()
{
	if (cl_lib) return true;
	for (int i = 0; opencl_lib_names[i]; i++)
	{
		cl_lib = LIB_OPEN(opencl_lib_names[i]);
		if (cl_lib) break;
	}
	if (!cl_lib) return false;

	// Resolve all required function pointers
	#define RESOLVE(fn) \
		CL.fn = (pfn_##fn)LIB_SYM(cl_lib, #fn); \
		if (!CL.fn) { LIB_CLOSE(cl_lib); cl_lib = NULL; return false; }

	RESOLVE(clGetPlatformIDs)
	RESOLVE(clGetPlatformInfo)
	RESOLVE(clGetDeviceIDs)
	RESOLVE(clGetDeviceInfo)
	RESOLVE(clCreateContext)
	RESOLVE(clCreateCommandQueue)
	RESOLVE(clCreateBuffer)
	RESOLVE(clCreateProgramWithSource)
	RESOLVE(clBuildProgram)
	RESOLVE(clGetProgramBuildInfo)
	RESOLVE(clCreateKernel)
	RESOLVE(clSetKernelArg)
	RESOLVE(clEnqueueNDRangeKernel)
	RESOLVE(clEnqueueWriteBuffer)
	RESOLVE(clEnqueueReadBuffer)
	RESOLVE(clFinish)
	RESOLVE(clReleaseMemObject)
	RESOLVE(clReleaseKernel)
	RESOLVE(clReleaseProgram)
	RESOLVE(clReleaseCommandQueue)
	RESOLVE(clReleaseContext)

	#undef RESOLVE

	// clEnqueueFillBuffer is optional (OpenCL 1.2+), not fatal if missing
	CL.clEnqueueFillBuffer = (pfn_clEnqueueFillBuffer)LIB_SYM(cl_lib, "clEnqueueFillBuffer");

	return true;
}


/// Round up n to the next multiple of wg
static inline size_t round_up_wg(size_t n, size_t wg)
{
	return ((n + wg - 1) / wg) * wg;
}

/// Release a single cl_mem buffer and set pointer to NULL
static void release_mem(cl_mem &m)
{
	if (m) { CL.clReleaseMemObject(m); m = NULL; }
}


// ========================================================================= //
// Public API implementation

/// Initialize the GPU context. Returns true if a usable GPU was found.
/// If OpenCL library is not available or no GPU device exists, returns false
bool gpu_opencl_init(bool verbose)
{
	if (gpu_initialized) return true;

	// Step 1: load the OpenCL library at runtime
	if (!load_opencl_library())
	{
		if (verbose) Rprintf("GPU: OpenCL library not found, using CPU.\n");
		return false;
	}

	// Step 2: find a GPU device
	cl_uint num_platforms = 0;
	if (CL.clGetPlatformIDs(0, NULL, &num_platforms) != CL_SUCCESS || num_platforms == 0)
	{
		if (verbose) Rprintf("GPU: no OpenCL platforms found, using CPU.\n");
		return false;
	}

	std::string chosen_platform_name;
	cl_platform_id *platforms = new cl_platform_id[num_platforms];
	CL.clGetPlatformIDs(num_platforms, platforms, NULL);

	bool found = false;
	for (cl_uint pi = 0; pi < num_platforms && !found; pi++)
	{
		cl_uint num_devs = 0;
		// Try GPU first, then fall back to any accelerator
		if (CL.clGetDeviceIDs(platforms[pi], CL_DEVICE_TYPE_GPU, 1, &gpu_device, &num_devs) == CL_SUCCESS
			&& num_devs > 0)
		{
			// Get platform name
			char pname[256] = {};
			CL.clGetPlatformInfo(platforms[pi], CL_PLATFORM_NAME, sizeof(pname), pname, NULL);
			chosen_platform_name = pname;

			cl_int err;
			gpu_context = CL.clCreateContext(NULL, 1, &gpu_device, NULL, NULL, &err);
			if (err != CL_SUCCESS) continue;
			gpu_queue = CL.clCreateCommandQueue(gpu_context, gpu_device, 0, &err);
			if (err != CL_SUCCESS) { CL.clReleaseContext(gpu_context); gpu_context = NULL; continue; }
			found = true;
		}
	}
	delete[] platforms;

	if (!found)
	{
		if (verbose) Rprintf("GPU: no OpenCL GPU device found, using CPU.\n");
		return false;
	}

	// Step 3: query device info
	char dev_name[256] = {};
	CL.clGetDeviceInfo(gpu_device, CL_DEVICE_NAME, sizeof(dev_name), dev_name, NULL);
	CL.clGetDeviceInfo(gpu_device, CL_DEVICE_MAX_WORK_GROUP_SIZE, sizeof(gpu_max_wg), &gpu_max_wg, NULL);

	cl_ulong dev_mem = 0;
	CL.clGetDeviceInfo(gpu_device, CL_DEVICE_GLOBAL_MEM_SIZE, sizeof(dev_mem), &dev_mem, NULL);

	cl_uint dev_cu = 0;
	CL.clGetDeviceInfo(gpu_device, CL_DEVICE_MAX_COMPUTE_UNITS, sizeof(dev_cu), &dev_cu, NULL);

	// Step 4: compile kernels
	cl_int err;
	const size_t src_len = strlen(kernel_source);
	gpu_program = CL.clCreateProgramWithSource(gpu_context, 1, &kernel_source, &src_len, &err);
	if (err != CL_SUCCESS)
	{
		if (verbose) Rprintf("GPU: failed to create OpenCL program (err=%d).\n", err);
		gpu_opencl_cleanup();
		return false;
	}

	// Build kernels (single-precision, no fp64 required)
	err = CL.clBuildProgram(gpu_program, 1, &gpu_device, "-cl-std=CL1.2", NULL, NULL);
	if (err != CL_SUCCESS)
	{
		if (verbose)
		{
			Rprintf("GPU: OpenCL kernel build failed (err=%d).\n", err);
			// Print build log
			size_t log_size = 0;
			CL.clGetProgramBuildInfo(gpu_program, gpu_device, CL_PROGRAM_BUILD_LOG, 0, NULL, &log_size);
			if (log_size > 0 && log_size < 8192)
			{
				char *log_buf = new char[log_size + 1];
				CL.clGetProgramBuildInfo(gpu_program, gpu_device, CL_PROGRAM_BUILD_LOG, log_size, log_buf, NULL);
				log_buf[log_size] = '\0';
				Rprintf("GPU build log:\n%s\n", log_buf);
				delete[] log_buf;
			}
		}
		gpu_opencl_cleanup();
		return false;
	}

	// Create kernels
	kern_dot_2b   = CL.clCreateKernel(gpu_program, "dot_2b", &err);
	if (err != CL_SUCCESS) { gpu_opencl_cleanup(); return false; }
	kern_accum_2b = CL.clCreateKernel(gpu_program, "accum_2b", &err);
	if (err != CL_SUCCESS) { gpu_opencl_cleanup(); return false; }
	kern_gemv     = CL.clCreateKernel(gpu_program, "gemv_dense", &err);
	if (err != CL_SUCCESS) { gpu_opencl_cleanup(); return false; }
	kern_grm_pair = CL.clCreateKernel(gpu_program, "grm_pairwise_2b", &err);
	if (err != CL_SUCCESS) { gpu_opencl_cleanup(); return false; }
	kern_grm_ijx  = CL.clCreateKernel(gpu_program, "grm_refine_ijx", &err);
	if (err != CL_SUCCESS) { gpu_opencl_cleanup(); return false; }

	gpu_initialized = true;

	// Choose local work-group size: min of device max and 256
	gpu_local_wg = (gpu_max_wg < 256) ? gpu_max_wg : 256;

	if (verbose)
	{
		Rprintf("GPU: OpenCL initialized\n");
		Rprintf("    Platform: %s\n", chosen_platform_name.c_str());
		Rprintf("    Device: %s\n", dev_name);
		Rprintf("    Compute units: %u, Global memory: %.1f GB\n",
			dev_cu, dev_mem / (1024.0*1024.0*1024.0));
		Rprintf("    Max/local work group size: %zu/%zu\n", gpu_max_wg, gpu_local_wg);
	}

	return true;
}


/// Release all GPU resources (context, queue, buffers, program)
void gpu_opencl_cleanup()
{
	gpu_free_buffers();
	if (kern_dot_2b)   { CL.clReleaseKernel(kern_dot_2b);   kern_dot_2b = NULL; }
	if (kern_accum_2b) { CL.clReleaseKernel(kern_accum_2b); kern_accum_2b = NULL; }
	if (kern_gemv)     { CL.clReleaseKernel(kern_gemv);     kern_gemv = NULL; }
	if (kern_grm_pair) { CL.clReleaseKernel(kern_grm_pair); kern_grm_pair = NULL; }
	if (kern_grm_ijx)  { CL.clReleaseKernel(kern_grm_ijx);  kern_grm_ijx = NULL; }
	if (gpu_program) { CL.clReleaseProgram(gpu_program);    gpu_program = NULL; }
	if (gpu_queue)   { CL.clReleaseCommandQueue(gpu_queue); gpu_queue = NULL; }
	if (gpu_context) { CL.clReleaseContext(gpu_context);    gpu_context = NULL; }
	gpu_device = NULL;
	gpu_initialized = false;
	gpu_has_2b = false;
	gpu_has_dense_grm = false;
}


/// Query whether GPU acceleration is currently available and initialized
bool gpu_opencl_available()
{
	return gpu_initialized;
}


/// Free GPU genotype buffers (called when GRM storage is reset)
void gpu_free_buffers()
{
	release_mem(buf_d_packed_geno);
	release_mem(buf_d_std_geno);
	release_mem(buf_d_b);
	release_mem(buf_d_out);
	release_mem(buf_d_dots);
	release_mem(buf_d_grm_mat);
	delete[] host_b_f;   host_b_f = NULL;
	delete[] host_out_f; host_out_f = NULL;
	host_f_len = 0;
	gpu_has_2b = false;
	gpu_has_dense_grm = false;
	gpu_nsamp = gpu_nvar = gpu_packed_nsamp = 0;
}


/// Upload dense 2-bit packed genotypes to GPU device memory
bool gpu_upload_2b_geno(const uint8_t *packed_geno, const double *std_geno,
	int nsamp, int nvar, int packed_nsamp)
{
	if (!gpu_initialized) return false;

	// Check memory requirements vs device memory (80% threshold)
	cl_ulong dev_mem = 0;
	CL.clGetDeviceInfo(gpu_device, CL_DEVICE_GLOBAL_MEM_SIZE, sizeof(dev_mem), &dev_mem, NULL);
	size_t needed = (size_t)packed_nsamp * nvar          // packed genotypes
		+ (size_t)4 * nvar * sizeof(float)               // std_geno lookup
		+ (size_t)nsamp * sizeof(float) * 2              // b + out vectors
		+ (size_t)nvar * sizeof(float);                  // dots vector
	if (needed > (size_t)(dev_mem * 0.8))
	{
		Rprintf("GPU: data size (%.1f GB) exceeds 80%% of device memory (%.1f GB), using CPU.\n",
			needed / (1024.0*1024.0*1024.0), dev_mem / (1024.0*1024.0*1024.0));
		return false;
	}

	// Release any previous buffers
	gpu_free_buffers();

	// Convert std_geno from double to float for GPU
	size_t std_geno_n = (size_t)4 * nvar;
	float *std_geno_f = new float[std_geno_n];
	for (size_t i = 0; i < std_geno_n; i++)
		std_geno_f[i] = (float)std_geno[i];

	cl_int err;
	buf_d_packed_geno = CL.clCreateBuffer(gpu_context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
		(size_t)packed_nsamp * nvar, (void*)packed_geno, &err);
	if (err != CL_SUCCESS) { delete[] std_geno_f; gpu_free_buffers(); return false; }

	buf_d_std_geno = CL.clCreateBuffer(gpu_context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
		std_geno_n * sizeof(float), (void*)std_geno_f, &err);
	delete[] std_geno_f;
	if (err != CL_SUCCESS) { gpu_free_buffers(); return false; }

	buf_d_b = CL.clCreateBuffer(gpu_context, CL_MEM_READ_ONLY,
		(size_t)nsamp * sizeof(float), NULL, &err);
	if (err != CL_SUCCESS) { gpu_free_buffers(); return false; }

	buf_d_out = CL.clCreateBuffer(gpu_context, CL_MEM_WRITE_ONLY,
		(size_t)nsamp * sizeof(float), NULL, &err);
	if (err != CL_SUCCESS) { gpu_free_buffers(); return false; }

	buf_d_dots = CL.clCreateBuffer(gpu_context, CL_MEM_READ_WRITE,
		(size_t)nvar * sizeof(float), NULL, &err);
	if (err != CL_SUCCESS) { gpu_free_buffers(); return false; }

	gpu_nsamp = nsamp;
	gpu_nvar = nvar;
	gpu_packed_nsamp = packed_nsamp;
	gpu_has_2b = true;

	// Allocate persistent host float buffers for compute calls
	if (host_f_len < nsamp)
	{
		delete[] host_b_f;   delete[] host_out_f;
		host_b_f   = new float[nsamp];
		host_out_f = new float[nsamp];
		host_f_len = nsamp;
	}

	return true;
}


/// Upload a user-defined dense GRM matrix to GPU device memory
bool gpu_upload_dense_grm(const double *mat, int nsamp)
{
	if (!gpu_initialized) return false;

	// Check memory
	cl_ulong dev_mem = 0;
	CL.clGetDeviceInfo(gpu_device, CL_DEVICE_GLOBAL_MEM_SIZE, sizeof(dev_mem), &dev_mem, NULL);
	size_t needed = (size_t)nsamp * nsamp * sizeof(float)    // GRM matrix
		+ (size_t)nsamp * sizeof(float) * 2;                 // b + out vectors
	if (needed > (size_t)(dev_mem * 0.8))
	{
		Rprintf("GPU: dense GRM size (%.1f GB) exceeds 80%% of device memory, using CPU.\n",
			needed / (1024.0*1024.0*1024.0));
		return false;
	}

	// Release previous dense GRM buffer (keep 2b buffers if any)
	release_mem(buf_d_grm_mat);

	// Convert dense GRM from double to float for GPU
	size_t mat_n = (size_t)nsamp * nsamp;
	float *mat_f = new float[mat_n];
	for (size_t i = 0; i < mat_n; i++)
		mat_f[i] = (float)mat[i];

	cl_int err;
	buf_d_grm_mat = CL.clCreateBuffer(gpu_context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
		mat_n * sizeof(float), (void*)mat_f, &err);
	delete[] mat_f;
	if (err != CL_SUCCESS) { release_mem(buf_d_grm_mat); return false; }

	// Ensure b and out buffers are allocated for this nsamp
	if (!buf_d_b || gpu_nsamp < nsamp)
	{
		release_mem(buf_d_b);
		release_mem(buf_d_out);
		buf_d_b = CL.clCreateBuffer(gpu_context, CL_MEM_READ_ONLY,
			(size_t)nsamp * sizeof(float), NULL, &err);
		if (err != CL_SUCCESS) { release_mem(buf_d_grm_mat); return false; }
		buf_d_out = CL.clCreateBuffer(gpu_context, CL_MEM_WRITE_ONLY,
			(size_t)nsamp * sizeof(float), NULL, &err);
		if (err != CL_SUCCESS) { release_mem(buf_d_grm_mat); release_mem(buf_d_b); return false; }
	}

	gpu_nsamp = nsamp;
	gpu_has_dense_grm = true;

	// Allocate persistent host float buffers for compute calls
	if (host_f_len < nsamp)
	{
		delete[] host_b_f;   delete[] host_out_f;
		host_b_f   = new float[nsamp];
		host_out_f = new float[nsamp];
		host_f_len = nsamp;
	}

	return true;
}


/// Perform GPU-accelerated GRM cross-product: out_b = G'G * b / nvar
/// b: input vector (nsamp-length)
/// out_b: output vector (nsamp-length), caller must allocate
/// Returns true if GPU computation succeeded; false to fall back to CPU
bool gpu_grm_crossprod_2b(const double *b, double *out_b, int nsamp,
	int nvar, int packed_nsamp)
{
	if (!gpu_initialized || !gpu_has_2b) return false;
	if (nsamp != gpu_nsamp || nvar != gpu_nvar || packed_nsamp != gpu_packed_nsamp)
		return false;

	// Convert b from double to float using persistent buffer
	for (int i = 0; i < nsamp; i++) host_b_f[i] = (float)b[i];

	cl_int err;

	// Upload b vector to device (float)
	err = CL.clEnqueueWriteBuffer(gpu_queue, buf_d_b, CL_TRUE, 0,
		(size_t)nsamp * sizeof(float), host_b_f, 0, NULL, NULL);
	if (err != CL_SUCCESS) return false;

	// Pass 1: compute dot products (one work item per variant)
	CL.clSetKernelArg(kern_dot_2b, 0, sizeof(cl_mem), &buf_d_packed_geno);
	CL.clSetKernelArg(kern_dot_2b, 1, sizeof(cl_mem), &buf_d_std_geno);
	CL.clSetKernelArg(kern_dot_2b, 2, sizeof(cl_mem), &buf_d_b);
	CL.clSetKernelArg(kern_dot_2b, 3, sizeof(cl_mem), &buf_d_dots);
	CL.clSetKernelArg(kern_dot_2b, 4, sizeof(int), &nsamp);
	CL.clSetKernelArg(kern_dot_2b, 5, sizeof(int), &packed_nsamp);
	CL.clSetKernelArg(kern_dot_2b, 6, sizeof(int), &nvar);

	size_t global_nvar = round_up_wg((size_t)nvar, gpu_local_wg);
	err = CL.clEnqueueNDRangeKernel(gpu_queue, kern_dot_2b, 1, NULL,
		&global_nvar, &gpu_local_wg, 0, NULL, NULL);
	if (err != CL_SUCCESS) return false;

	// Pass 2: accumulate into output (one work item per sample)
	float inv_nvar_f = 1.0f / nvar;
	CL.clSetKernelArg(kern_accum_2b, 0, sizeof(cl_mem), &buf_d_packed_geno);
	CL.clSetKernelArg(kern_accum_2b, 1, sizeof(cl_mem), &buf_d_std_geno);
	CL.clSetKernelArg(kern_accum_2b, 2, sizeof(cl_mem), &buf_d_dots);
	CL.clSetKernelArg(kern_accum_2b, 3, sizeof(cl_mem), &buf_d_out);
	CL.clSetKernelArg(kern_accum_2b, 4, sizeof(int), &nsamp);
	CL.clSetKernelArg(kern_accum_2b, 5, sizeof(int), &nvar);
	CL.clSetKernelArg(kern_accum_2b, 6, sizeof(int), &packed_nsamp);
	CL.clSetKernelArg(kern_accum_2b, 7, sizeof(float), &inv_nvar_f);

	size_t global_nsamp = round_up_wg((size_t)nsamp, gpu_local_wg);
	err = CL.clEnqueueNDRangeKernel(gpu_queue, kern_accum_2b, 1, NULL,
		&global_nsamp, &gpu_local_wg, 0, NULL, NULL);
	if (err != CL_SUCCESS) return false;

	// Read back result (float) and convert to double
	err = CL.clEnqueueReadBuffer(gpu_queue, buf_d_out, CL_TRUE, 0,
		(size_t)nsamp * sizeof(float), host_out_f, 0, NULL, NULL);
	if (err != CL_SUCCESS) return false;

	for (int i = 0; i < nsamp; i++) out_b[i] = (double)host_out_f[i];

	return true;
}


/// Perform GPU-accelerated dense GRM matrix-vector product: out_b = mat * b
/// b: input vector (nsamp-length)
/// out_b: output vector (nsamp-length), caller must allocate
/// Returns true if GPU computation succeeded; false to fall back to CPU
bool gpu_grm_crossprod_dense(const double *b, double *out_b, int nsamp)
{
	if (!gpu_initialized || !gpu_has_dense_grm) return false;
	if (nsamp != gpu_nsamp) return false;

	// Convert b from double to float using persistent buffer
	for (int i = 0; i < nsamp; i++) host_b_f[i] = (float)b[i];

	cl_int err;

	// Upload b (float)
	err = CL.clEnqueueWriteBuffer(gpu_queue, buf_d_b, CL_TRUE, 0,
		(size_t)nsamp * sizeof(float), host_b_f, 0, NULL, NULL);
	if (err != CL_SUCCESS) return false;

	// Launch GEMV kernel
	CL.clSetKernelArg(kern_gemv, 0, sizeof(cl_mem), &buf_d_grm_mat);
	CL.clSetKernelArg(kern_gemv, 1, sizeof(cl_mem), &buf_d_b);
	CL.clSetKernelArg(kern_gemv, 2, sizeof(cl_mem), &buf_d_out);
	CL.clSetKernelArg(kern_gemv, 3, sizeof(int), &nsamp);

	size_t global_nsamp = round_up_wg((size_t)nsamp, gpu_local_wg);
	err = CL.clEnqueueNDRangeKernel(gpu_queue, kern_gemv, 1, NULL,
		&global_nsamp, &gpu_local_wg, 0, NULL, NULL);
	if (err != CL_SUCCESS) return false;

	// Read back (float) and convert to double
	err = CL.clEnqueueReadBuffer(gpu_queue, buf_d_out, CL_TRUE, 0,
		(size_t)nsamp * sizeof(float), host_out_f, 0, NULL, NULL);
	if (err != CL_SUCCESS) return false;

	for (int i = 0; i < nsamp; i++) out_b[i] = (double)host_out_f[i];

	return true;
}


// ========================================================================= //
// GRM construction API

/// Helper: convert the R-side g_lookup (8 x nrow*4 double matrix with pairwise
/// products) into a flat float[nvar*4] array of per-variant standardized
/// genotype values (for codes 0,1,2,3).
/// The R-side lookup stores: for each group of 4 SNPs (one byte position),
/// 8 entries encode pairwise products. We need individual standardized values.
/// From the initialization in grm_sp_init_lookup (saige_misc.cpp):
///   F[0]=g0*g0, F[1]=g0*g1, F[2]=0, F[3]=g1*g1,
///   F[4]=g0*g2, F[5]=g1*g2, F[6]=g2*g2, F[7]=0
/// So: g0 = sqrt(F[0]), g1 = F[1]/g0 (if g0!=0), g2 = F[4]/g0 (if g0!=0)
/// But we actually need to recover individual values, not products.
/// Instead, we compute them directly from the packed genotypes.
///
/// Simpler approach: compute std_geno directly from g_pack on the host side.
static float *compute_std_geno_f32(const uint8_t *g_pack, int n_byte, int nsamp,
	int nvar)
{
	float *std_geno = new float[(size_t)nvar * 4];
	for (int v = 0; v < nvar; v++)
	{
		int bp = v >> 2;      // byte position
		int bs = (v & 3) << 1; // bit shift
		int n = 0, s = 0;
		for (int j = 0; j < nsamp; j++)
		{
			uint8_t byte_val = g_pack[bp + (size_t)j * n_byte];
			int c = (byte_val >> bs) & 0x03;
			if (c < 3) { n++; s += c; }
		}
		double f = (n > 0) ? (double)s / (2 * n) : 0.0;
		double d = 1.0 / sqrt(2.0 * f * (1.0 - f));
		if (!R_FINITE(d)) { d = 0; f = 0; }
		double f2 = 2.0 * f;
		float g0 = (float)((0 - f2) * d);
		float g1 = (float)((1 - f2) * d);
		float g2 = (float)((2 - f2) * d);
		std_geno[v * 4 + 0] = g0;
		std_geno[v * 4 + 1] = g1;
		std_geno[v * 4 + 2] = g2;
		std_geno[v * 4 + 3] = 0.0f;  // missing -> 0
	}
	return std_geno;
}


/// Compute dense GRM on GPU. g_pack: [n_byte x nsamp] raw matrix (col-major).
/// grm_out: caller-allocated nsamp x nsamp double array (col-major).
/// Returns true on success.
bool gpu_grm_dense_calc(const uint8_t *g_pack, int nsamp, int nvar,
	int n_byte, bool use_f64, int bs, double *grm_out,
	gpu_progress_fn prog_cb, void *prog_data)
{
	if (!gpu_initialized) return false;

	// Check GPU memory: need packed geno + std_geno + output block buffer
	cl_ulong dev_mem = 0;
	CL.clGetDeviceInfo(gpu_device, CL_DEVICE_GLOBAL_MEM_SIZE, sizeof(dev_mem), &dev_mem, NULL);
	size_t geno_size = (size_t)n_byte * nsamp;
	size_t std_size = (size_t)nvar * 4 * sizeof(float);
	size_t block_out_size = (size_t)bs * bs * sizeof(float);
	size_t needed = geno_size + std_size + block_out_size;
	if (needed > (size_t)(dev_mem * 0.8))
	{
		Rprintf("GPU: dense GRM data (%.1f GB) exceeds 80%% of device memory, using CPU.\n",
			needed / (1024.0*1024.0*1024.0));
		return false;
	}

	// Compute standardized genotype lookup on host
	float *std_geno_f = compute_std_geno_f32(g_pack, n_byte, nsamp, nvar);

	cl_int err;
	// Upload packed genotypes
	cl_mem d_geno = CL.clCreateBuffer(gpu_context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
		geno_size, (void*)g_pack, &err);
	if (err != CL_SUCCESS) { delete[] std_geno_f; return false; }

	// Upload std_geno
	cl_mem d_std = CL.clCreateBuffer(gpu_context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
		std_size, (void*)std_geno_f, &err);
	delete[] std_geno_f;
	if (err != CL_SUCCESS) { CL.clReleaseMemObject(d_geno); return false; }

	// Allocate output block buffer on device
	cl_mem d_out = CL.clCreateBuffer(gpu_context, CL_MEM_WRITE_ONLY,
		(size_t)bs * bs * sizeof(float), NULL, &err);
	if (err != CL_SUCCESS)
	{
		CL.clReleaseMemObject(d_geno);
		CL.clReleaseMemObject(d_std);
		return false;
	}

	// Host buffer for reading back block results
	float *h_block = new float[(size_t)bs * bs];

	// Process blocks over upper triangle
	int n_block = (nsamp + bs - 1) / bs;
	bool ok = true;

	// Set loop-invariant kernel args once before the block loop
	CL.clSetKernelArg(kern_grm_pair, 0, sizeof(cl_mem), &d_geno);
	CL.clSetKernelArg(kern_grm_pair, 1, sizeof(cl_mem), &d_std);
	CL.clSetKernelArg(kern_grm_pair, 2, sizeof(cl_mem), &d_out);
	CL.clSetKernelArg(kern_grm_pair, 3, sizeof(int), &nsamp);
	CL.clSetKernelArg(kern_grm_pair, 4, sizeof(int), &nvar);
	CL.clSetKernelArg(kern_grm_pair, 5, sizeof(int), &n_byte);

	for (int bi = 0; bi < n_block && ok; bi++)
	{
		int i_st = bi * bs;
		int i_n = std::min(i_st + bs, nsamp) - i_st;
		for (int bj = bi; bj < n_block && ok; bj++)
		{
			int j_st = bj * bs;
			int j_n = std::min(j_st + bs, nsamp) - j_st;

			// Only per-block args need updating each iteration
			CL.clSetKernelArg(kern_grm_pair, 6, sizeof(int), &i_st);
			CL.clSetKernelArg(kern_grm_pair, 7, sizeof(int), &i_n);
			CL.clSetKernelArg(kern_grm_pair, 8, sizeof(int), &j_st);
			CL.clSetKernelArg(kern_grm_pair, 9, sizeof(int), &j_n);

			// 2D launch
			size_t local_wg = 16;
			size_t global[2] = { round_up_wg(i_n, local_wg), round_up_wg(j_n, local_wg) };
			size_t local[2] = { local_wg, local_wg };
			err = CL.clEnqueueNDRangeKernel(gpu_queue, kern_grm_pair, 2, NULL,
				global, local, 0, NULL, NULL);
			if (err != CL_SUCCESS) { ok = false; break; }

			// Read back block
			err = CL.clEnqueueReadBuffer(gpu_queue, d_out, CL_TRUE, 0,
				(size_t)i_n * j_n * sizeof(float), h_block, 0, NULL, NULL);
			if (err != CL_SUCCESS) { ok = false; break; }

			// Copy to output matrix (col-major)
			for (int li = 0; li < i_n; li++)
			{
				int gi = i_st + li;
				for (int lj = 0; lj < j_n; lj++)
				{
					int gj = j_st + lj;
					if (gj >= gi)
					{
						double v = (double)h_block[li * j_n + lj];
						grm_out[gi + (size_t)gj * nsamp] = v;
						grm_out[gj + (size_t)gi * nsamp] = v;
					}
				}
			}
			// update progress
			if (prog_cb) prog_cb(prog_data);
		}
	}

	delete[] h_block;
	CL.clReleaseMemObject(d_out);
	CL.clReleaseMemObject(d_std);
	CL.clReleaseMemObject(d_geno);

	return ok;
}


/// Sparse GRM step 1: scan all pairs in blocks, return (i,j) pairs >= rel_cutoff.
/// g_pack: [n_byte x nsamp], std_geno computed internally.
/// out_i, out_j: output vectors of 0-based sample indices.
/// Returns true on success.
bool gpu_grm_sparse_scan(const uint8_t *g_pack, int nsamp, int nvar,
	int n_byte, double rel_cutoff, int bs, std::vector<int> &out_i,
	std::vector<int> &out_j,
	gpu_progress_fn prog_cb, void *prog_data)
{
	if (!gpu_initialized) return false;

	cl_ulong dev_mem = 0;
	CL.clGetDeviceInfo(gpu_device, CL_DEVICE_GLOBAL_MEM_SIZE, sizeof(dev_mem), &dev_mem, NULL);
	size_t geno_size = (size_t)n_byte * nsamp;
	size_t std_size = (size_t)nvar * 4 * sizeof(float);
	size_t block_out_size = (size_t)bs * bs * sizeof(float);
	size_t needed = geno_size + std_size + block_out_size;
	if (needed > (size_t)(dev_mem * 0.8))
		return false;

	float *std_geno_f = compute_std_geno_f32(g_pack, n_byte, nsamp, nvar);

	cl_int err;
	cl_mem d_geno = CL.clCreateBuffer(gpu_context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
		geno_size, (void*)g_pack, &err);
	if (err != CL_SUCCESS) { delete[] std_geno_f; return false; }

	cl_mem d_std = CL.clCreateBuffer(gpu_context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
		std_size, (void*)std_geno_f, &err);
	delete[] std_geno_f;
	if (err != CL_SUCCESS) { CL.clReleaseMemObject(d_geno); return false; }

	cl_mem d_out = CL.clCreateBuffer(gpu_context, CL_MEM_WRITE_ONLY,
		(size_t)bs * bs * sizeof(float), NULL, &err);
	if (err != CL_SUCCESS)
	{
		CL.clReleaseMemObject(d_geno);
		CL.clReleaseMemObject(d_std);
		return false;
	}

	float *h_block = new float[(size_t)bs * bs];
	// Use slightly relaxed cutoff to avoid float32 false negatives
	float rel_cutoff_f = (float)(rel_cutoff * 0.99);

	out_i.clear();
	out_j.clear();

	int n_block = (nsamp + bs - 1) / bs;
	bool ok = true;

	// Set loop-invariant kernel args once before the block loop
	CL.clSetKernelArg(kern_grm_pair, 0, sizeof(cl_mem), &d_geno);
	CL.clSetKernelArg(kern_grm_pair, 1, sizeof(cl_mem), &d_std);
	CL.clSetKernelArg(kern_grm_pair, 2, sizeof(cl_mem), &d_out);
	CL.clSetKernelArg(kern_grm_pair, 3, sizeof(int), &nsamp);
	CL.clSetKernelArg(kern_grm_pair, 4, sizeof(int), &nvar);
	CL.clSetKernelArg(kern_grm_pair, 5, sizeof(int), &n_byte);

	for (int bi = 0; bi < n_block && ok; bi++)
	{
		int i_st = bi * bs;
		int i_n = std::min(i_st + bs, nsamp) - i_st;
		for (int bj = bi; bj < n_block && ok; bj++)
		{
			int j_st = bj * bs;
			int j_n = std::min(j_st + bs, nsamp) - j_st;

			// Only per-block args need updating each iteration
			CL.clSetKernelArg(kern_grm_pair, 6, sizeof(int), &i_st);
			CL.clSetKernelArg(kern_grm_pair, 7, sizeof(int), &i_n);
			CL.clSetKernelArg(kern_grm_pair, 8, sizeof(int), &j_st);
			CL.clSetKernelArg(kern_grm_pair, 9, sizeof(int), &j_n);

			size_t local_wg = 16;
			size_t global[2] = { round_up_wg(i_n, local_wg), round_up_wg(j_n, local_wg) };
			size_t local[2] = { local_wg, local_wg };
			err = CL.clEnqueueNDRangeKernel(gpu_queue, kern_grm_pair, 2, NULL,
				global, local, 0, NULL, NULL);
			if (err != CL_SUCCESS) { ok = false; break; }

			err = CL.clEnqueueReadBuffer(gpu_queue, d_out, CL_TRUE, 0,
				(size_t)i_n * j_n * sizeof(float), h_block, 0, NULL, NULL);
			if (err != CL_SUCCESS) { ok = false; break; }

			// Filter on host
			for (int li = 0; li < i_n; li++)
			{
				int gi = i_st + li;
				for (int lj = 0; lj < j_n; lj++)
				{
					int gj = j_st + lj;
					if (gj >= gi && h_block[li * j_n + lj] >= rel_cutoff_f)
					{
						out_i.push_back(gi);
						out_j.push_back(gj);
					}
				}
			}
			// update progress
			if (prog_cb) prog_cb(prog_data);
		}
	}

	delete[] h_block;
	CL.clReleaseMemObject(d_out);
	CL.clReleaseMemObject(d_std);
	CL.clReleaseMemObject(d_geno);

	return ok;
}


/// Sparse GRM step 2: compute GRM values for specific (i,j) pairs using full
/// variant set on GPU.
/// g_pack: [n_byte x nsamp], pair_i/pair_j: 0-based indices, npairs: count.
/// out_x: caller-allocated output array of length npairs.
/// Returns true on success.
bool gpu_grm_sparse_refine(const uint8_t *g_pack, int nsamp, int nvar,
	int n_byte, const int *pair_i, const int *pair_j, int npairs,
	double *out_x)
{
	if (!gpu_initialized || npairs <= 0) return false;

	cl_ulong dev_mem = 0;
	CL.clGetDeviceInfo(gpu_device, CL_DEVICE_GLOBAL_MEM_SIZE, sizeof(dev_mem), &dev_mem, NULL);
	size_t geno_size = (size_t)n_byte * nsamp;
	size_t std_size = (size_t)nvar * 4 * sizeof(float);
	size_t pair_size = (size_t)npairs * sizeof(int) * 2 + (size_t)npairs * sizeof(float);
	size_t needed = geno_size + std_size + pair_size;
	if (needed > (size_t)(dev_mem * 0.8))
		return false;

	float *std_geno_f = compute_std_geno_f32(g_pack, n_byte, nsamp, nvar);

	cl_int err;
	cl_mem d_geno = CL.clCreateBuffer(gpu_context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
		geno_size, (void*)g_pack, &err);
	if (err != CL_SUCCESS) { delete[] std_geno_f; return false; }

	cl_mem d_std = CL.clCreateBuffer(gpu_context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
		std_size, (void*)std_geno_f, &err);
	delete[] std_geno_f;
	if (err != CL_SUCCESS) { CL.clReleaseMemObject(d_geno); return false; }

	cl_mem d_pi = CL.clCreateBuffer(gpu_context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
		(size_t)npairs * sizeof(int), (void*)pair_i, &err);
	if (err != CL_SUCCESS)
	{
		CL.clReleaseMemObject(d_geno); CL.clReleaseMemObject(d_std);
		return false;
	}

	cl_mem d_pj = CL.clCreateBuffer(gpu_context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
		(size_t)npairs * sizeof(int), (void*)pair_j, &err);
	if (err != CL_SUCCESS)
	{
		CL.clReleaseMemObject(d_geno); CL.clReleaseMemObject(d_std);
		CL.clReleaseMemObject(d_pi);
		return false;
	}

	cl_mem d_x = CL.clCreateBuffer(gpu_context, CL_MEM_WRITE_ONLY,
		(size_t)npairs * sizeof(float), NULL, &err);
	if (err != CL_SUCCESS)
	{
		CL.clReleaseMemObject(d_geno); CL.clReleaseMemObject(d_std);
		CL.clReleaseMemObject(d_pi); CL.clReleaseMemObject(d_pj);
		return false;
	}

	// Set kernel args
	CL.clSetKernelArg(kern_grm_ijx, 0, sizeof(cl_mem), &d_geno);
	CL.clSetKernelArg(kern_grm_ijx, 1, sizeof(cl_mem), &d_std);
	CL.clSetKernelArg(kern_grm_ijx, 2, sizeof(cl_mem), &d_pi);
	CL.clSetKernelArg(kern_grm_ijx, 3, sizeof(cl_mem), &d_pj);
	CL.clSetKernelArg(kern_grm_ijx, 4, sizeof(cl_mem), &d_x);
	CL.clSetKernelArg(kern_grm_ijx, 5, sizeof(int), &nvar);
	CL.clSetKernelArg(kern_grm_ijx, 6, sizeof(int), &n_byte);
	CL.clSetKernelArg(kern_grm_ijx, 7, sizeof(int), &npairs);

	size_t global_n = round_up_wg((size_t)npairs, gpu_local_wg);
	err = CL.clEnqueueNDRangeKernel(gpu_queue, kern_grm_ijx, 1, NULL,
		&global_n, &gpu_local_wg, 0, NULL, NULL);
	if (err != CL_SUCCESS)
	{
		CL.clReleaseMemObject(d_x); CL.clReleaseMemObject(d_pj);
		CL.clReleaseMemObject(d_pi); CL.clReleaseMemObject(d_std);
		CL.clReleaseMemObject(d_geno);
		return false;
	}

	// Read back results
	float *h_x = new float[npairs];
	err = CL.clEnqueueReadBuffer(gpu_queue, d_x, CL_TRUE, 0,
		(size_t)npairs * sizeof(float), h_x, 0, NULL, NULL);

	if (err == CL_SUCCESS)
	{
		for (int i = 0; i < npairs; i++)
			out_x[i] = (double)h_x[i];
	}

	delete[] h_x;
	CL.clReleaseMemObject(d_x);
	CL.clReleaseMemObject(d_pj);
	CL.clReleaseMemObject(d_pi);
	CL.clReleaseMemObject(d_std);
	CL.clReleaseMemObject(d_geno);

	return (err == CL_SUCCESS);
}


// ========================================================================= //
// R-callable entry points

extern "C" SEXP saige_gpu_init(SEXP verbose)
{
	bool v = (Rf_asLogical(verbose) == TRUE);
	bool ok = gpu_opencl_init(v);
	return Rf_ScalarLogical(ok ? TRUE : FALSE);
}

extern "C" SEXP saige_gpu_cleanup()
{
	gpu_opencl_cleanup();
	return R_NilValue;
}
