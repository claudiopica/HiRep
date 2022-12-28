/***************************************************************************\
 * Copyright (c) 2012, Ari Hietanen and Ulrik Søndergaard                   *
 * All rights reserved.                                                     *
\***************************************************************************/
/* Copyright (c) 2022, NVIDIA CORPORATION. All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *  * Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *  * Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 *  * Neither the name of NVIDIA CORPORATION nor the names of its
 *    contributors may be used to endorse or promote products derived
 *    from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS ``AS IS'' AND ANY
 * EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
 * OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */


#ifdef WITH_GPU
//This file should not be compiled if !WITH_GPU

#include "inverters.h"
#include "libhr_core.h"

#define GSUM_BLOCK_SIZE 256     // No more than 1024 on Tesla
#define BLOCK_SIZE_REM 64

/************************************************************************************************/
/************************************************************************************************/
/************************************************************************************************/

// Utility class used to avoid linker errors with extern
// unsized shared memory arrays with templated type
template <class T>
struct SharedMemory {
  __device__ inline operator T *() {
    extern __shared__ int __smem[];
    return (T *)__smem;
  }

  __device__ inline operator const T *() const {
    extern __shared__ int __smem[];
    return (T *)__smem;
  }
};

// specialize for double to avoid unaligned memory
// access compile errors
template <>
struct SharedMemory<double> {
  __device__ inline operator double *() {
    extern __shared__ double __smem_d[];
    return (double *)__smem_d;
  }

  __device__ inline operator const double *() const {
    extern __shared__ double __smem_d[];
    return (double *)__smem_d;
  }
};
// For hr_complex
template <>
struct SharedMemory<hr_complex> {
  __device__ inline operator hr_complex *() {
    extern __shared__ double __smem_d[];
    return (hr_complex *)__smem_d;
  }

  __device__ inline operator const hr_complex *() const {
    extern __shared__ double __smem_d[];
    return (hr_complex *)__smem_d;
  }
};

template <class T>
__device__ __forceinline__ T warpReduceSum(unsigned int mask, T mySum) {
  for (int offset = warpSize / 2; offset > 0; offset /= 2) {
    mySum += __shfl_down_sync(mask, mySum, offset);
  }
  return mySum;
}

template <>
__device__ __forceinline__ hr_complex warpReduceSum(unsigned int mask, hr_complex mySum) {
  for (int offset = warpSize / 2; offset > 0; offset /= 2) {
    mySum.re += __shfl_down_sync(mask, creal(mySum), offset);
    mySum.im += __shfl_down_sync(mask, cimag(mySum), offset);
  }
  return mySum;
}

template <>
__device__ __forceinline__ hr_complex_flt warpReduceSum(unsigned int mask, hr_complex_flt mySum) {
  for (int offset = warpSize / 2; offset > 0; offset /= 2) {
    mySum.re += __shfl_down_sync(mask, creal(mySum), offset);
    mySum.im += __shfl_down_sync(mask, cimag(mySum), offset);
  }
  return mySum;
}

#if __CUDA_ARCH__ >= 800
// Specialize warpReduceFunc for int inputs to use __reduce_add_sync intrinsic
// when on SM 8.0 or higher
template <>
__device__ __forceinline__ int warpReduceSum<int>(unsigned int mask, int mySum) {
  mySum = __reduce_add_sync(mask, mySum);
  return mySum;
}
#endif

/*
    Parallel sum reduction using shared memory
    - takes log(n) steps for n input elements
    - uses n threads
    - only works for power-of-2 arrays

    This version adds multiple elements per thread sequentially.  This reduces
   the overall cost of the algorithm while keeping the work complexity O(n) and
   the step complexity O(log n). (Brent's Theorem optimization)
    Note, this kernel needs a minimum of 64*sizeof(T) bytes of shared memory.
    In other words if blockSize <= 32, allocate 64*sizeof(T) bytes.
    If blockSize > 32, allocate blockSize*sizeof(T) bytes.
*/
template <typename T, unsigned int blockSize, bool nIsPow2>
__global__ void reduce7(const T *__restrict__ g_idata, T *__restrict__ g_odata, unsigned int n) {
  T *sdata = SharedMemory<T>();

  // perform first level of reduction,
  // reading from global memory, writing to shared memory
  unsigned int tid = threadIdx.x;
  unsigned int gridSize = blockSize * gridDim.x;
  unsigned int maskLength = (blockSize & 31);  // 31 = warpSize-1
  maskLength = (maskLength > 0) ? (32 - maskLength) : maskLength;
  const unsigned int mask = (0xffffffff) >> maskLength;

  T mySum = 0;

  // we reduce multiple elements per thread.  The number is determined by the
  // number of active thread blocks (via gridDim).  More blocks will result
  // in a larger gridSize and therefore fewer elements per thread
  if (nIsPow2) {
    unsigned int i = blockIdx.x * blockSize * 2 + threadIdx.x;
    gridSize = gridSize << 1;

    while (i < n) {
      mySum += g_idata[i];
      // ensure we don't read out of bounds -- this is optimized away for
      // powerOf2 sized arrays
      if ((i + blockSize) < n) {
        mySum += g_idata[i + blockSize];
      }
      i += gridSize;
    }
  } else {
    unsigned int i = blockIdx.x * blockSize + threadIdx.x;
    while (i < n) {
      mySum += g_idata[i];
      i += gridSize;
    }
  }

  // Reduce within warp using shuffle or reduce_add if T==int & CUDA_ARCH ==
  // SM 8.0
  mySum = warpReduceSum<T>(mask, mySum);

  // each thread puts its local sum into shared memory
  if ((tid % warpSize) == 0) {
    sdata[tid / warpSize] = mySum;
  }

  __syncthreads();

  const unsigned int shmem_extent =
      (blockSize / warpSize) > 0 ? (blockSize / warpSize) : 1;
  const unsigned int ballot_result = __ballot_sync(mask, tid < shmem_extent);
  if (tid < shmem_extent) {
    mySum = sdata[tid];
    // Reduce final warp using shuffle or reduce_add if T==int & CUDA_ARCH ==
    // SM 8.0
    mySum = warpReduceSum<T>(ballot_result, mySum);
  }

  // write result for this block to global mem
  if (tid == 0) {
    g_odata[blockIdx.x] = mySum;
  }
}

static bool isPow2(unsigned int x) { return ((x & (x - 1)) == 0); }

////////////////////////////////////////////////////////////////////////////////
// Wrapper function for kernel launch
////////////////////////////////////////////////////////////////////////////////
template <class T>
void reduce(int size, int threads, int blocks, T *d_idata, T *d_odata) {
  dim3 dimBlock(threads, 1, 1);
  dim3 dimGrid(blocks, 1, 1);

  // when there is only one warp per block, we need to allocate two warps
  // worth of shared memory so that we don't index shared memory out of bounds
  int smemSize = (threads <= 32) ? 2 * threads * sizeof(T) : threads * sizeof(T);

  // For reduce7 kernel we require only blockSize/warpSize
  // number of elements in shared memory
  smemSize = ((threads / 32) + 1) * sizeof(T);
  if (isPow2(size)) {
    switch (threads) {
      case 1024:
        reduce7<T, 1024, true><<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata, size);
        break;
      case 512:
        reduce7<T, 512, true><<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata, size);
        break;
      case 256:
        reduce7<T, 256, true><<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata, size);
        break;
      case 128:
        reduce7<T, 128, true><<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata, size);
        break;
      case 64:
        reduce7<T, 64, true><<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata, size);
        break;
      case 32:
        reduce7<T, 32, true><<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata, size);
        break;
      case 16:
        reduce7<T, 16, true><<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata, size);
        break;
      case 8:
        reduce7<T, 8, true><<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata, size);
        break;
      case 4:
        reduce7<T, 4, true><<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata, size);
        break;
      case 2:
        reduce7<T, 2, true><<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata, size);
        break;
      case 1:
        reduce7<T, 1, true><<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata, size);
        break;
      }
  } else {
    switch (threads) {
      case 1024:
        reduce7<T, 1024, false><<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata, size);
        break;
      case 512:
        reduce7<T, 512, false><<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata, size);
        break;
      case 256:
        reduce7<T, 256, false><<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata, size);
        break;
      case 128:
        reduce7<T, 128, false><<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata, size);
        break;
      case 64:
        reduce7<T, 64, false><<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata, size);
        break;
      case 32:
        reduce7<T, 32, false><<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata, size);
        break;
      case 16:
        reduce7<T, 16, false><<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata, size);
        break;
      case 8:
        reduce7<T, 8, false><<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata, size);
        break;
      case 4:
        reduce7<T, 4, false><<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata, size);
        break;
      case 2:
        reduce7<T, 2, false><<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata, size);
        break;
      case 1:
        reduce7<T, 1, false><<<dimGrid, dimBlock, smemSize>>>(d_idata, d_odata, size);
        break;
    }
  }
}

template <class T>
T global_sum_gpu(T *vector, int size) {
  int threads = 1024; while (size<threads) { threads /= 2; }
  // Blocks calculation from Nvidia reduction sample, see subroutine getNumBlocksAndThreads()
  const int blocks = (size + (threads * 2 - 1)) / (threads * 2);
  T res = 0;
  T *vector_host = (T *)malloc(blocks * sizeof(T));

  // Reduction over blocks
  T *vector_out = NULL;
  cudaMalloc((void **)&vector_out, blocks * sizeof(T));

  reduce<T>(size, threads, blocks, vector, vector_out);
  cudaMemcpy(vector_host, vector_out, blocks * sizeof(T), cudaMemcpyDeviceToHost);

  for (int i = 0; i < blocks; i++) {
    res += vector_host[i];
  }

  // free and return
  free(vector_host);
  cudaFree(vector_out);
  return res;
}

template int global_sum_gpu<int>(int* vector, int size);
template float global_sum_gpu<float>(float* vector, int size);
template double global_sum_gpu<double>(double* vector, int size);
template hr_complex_flt global_sum_gpu<hr_complex_flt>(hr_complex_flt* vector, int size);
template hr_complex global_sum_gpu<hr_complex>(hr_complex* vector, int size);


#ifdef __cplusplus
  extern "C" {
#endif
// The following function is to expose the global sum to C code
int global_sum_gpu_int(int* vector, int size){
  int res;
  int* vector_d;
  cudaMalloc((void **)&vector_d, size*sizeof(int));
  cudaMemcpy(vector_d, vector, size*sizeof(int), cudaMemcpyHostToDevice);
  res = global_sum_gpu<int>(vector_d, size);
  cudaFree(vector_d);
  return res;
}

// The following function is to expose the global sum to C code
float global_sum_gpu_float(float* vector, int size){
  float res;
  float* vector_d;
  cudaMalloc((void **)&vector_d, size*sizeof(float));
  cudaMemcpy(vector_d, vector, size*sizeof(float), cudaMemcpyHostToDevice);
  res = global_sum_gpu<float>(vector_d, size);
  cudaFree(vector_d);
  return res;
}

// The following function is to expose the global sum to C code
double global_sum_gpu_double(double* vector, int size){
  double res;
  double* vector_d;
  cudaMalloc((void **)&vector_d, size*sizeof(double));
  cudaMemcpy(vector_d, vector, size*sizeof(double), cudaMemcpyHostToDevice);
  res = global_sum_gpu<double>(vector_d, size);
  cudaFree(vector_d);
  return res;
}

// The following function is to expose the global sum to C code
hr_complex_flt global_sum_gpu_complex_flt(hr_complex_flt* vector, int size){
  hr_complex_flt res;
  hr_complex_flt* vector_d;
  cudaMalloc((void **)&vector_d, size*sizeof(hr_complex_flt));
  cudaMemcpy(vector_d, vector, size*sizeof(hr_complex_flt), cudaMemcpyHostToDevice);
  res = global_sum_gpu<hr_complex_flt>(vector_d, size);
  cudaFree(vector_d);
  return res;
}

// The following function is to expose the global sum to C code
hr_complex global_sum_gpu_complex(hr_complex* vector, int size){
  hr_complex res;
  hr_complex* vector_d;
  cudaMalloc((void **)&vector_d, size*sizeof(hr_complex));
  cudaMemcpy(vector_d, vector, size*sizeof(hr_complex), cudaMemcpyHostToDevice);
  res = global_sum_gpu<hr_complex>(vector_d, size);
  cudaFree(vector_d);
  return res;
}

#ifdef __cplusplus
  }
#endif

#endif
