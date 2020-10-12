/***************************************************************************\
 * Copyright (c) 2012, Ari Hietanen and Ulrik Søndergaard                   *   
 * All rights reserved.                                                     * 
\***************************************************************************/


#ifndef GLOBAL_SUM_GPU_C
#define GLOBAL_SUM_GPU_C
#include "hr_complex.h"
#include "global.h"
 #include "gpu.h"

#define GSUM_BLOCK_SIZE 256     // No more than 1024 on Tesla
#define BLOCK_SIZE_REM 64

template<typename REAL>
__global__ void copy_with_zero_padding(REAL* out, REAL* in, int N, int new_n){
  int i = blockIdx.x*BLOCK_SIZE + threadIdx.x;
  i = min(i,new_n);
  if (i<N){
    out[i]=in[i];
  }
  else {
    out[i]=0.0;
  }
}

template<typename COMPLEX>
__global__ void copy_with_zero_padding_complex(COMPLEX* out, COMPLEX* in, int N, int new_n){
  int i = blockIdx.x*BLOCK_SIZE + threadIdx.x;
  i = min(i,new_n);
  if (i<N){
    out[i]=in[i];
  }
  else {
    out[i].re=0.0;
    out[i].im=0.0;
  }
}

unsigned int next_pow2( unsigned int n )
{
  register unsigned int val = n;
  --val;
  val |= (val >> 1);
  val |= (val >> 2);
  val |= (val >> 4);
  val |= (val >> 8);
  val |= (val >> 16);
  return ++val; // Val is now the next highest power of 2.
} 



////////////////////// Optimized Global Sum Kernel for REALS //////////////////////////
// Requires multiple kernel calls in order to produce global synchronization 
// (You might think that this gives some overhead, but you are wrong.)
//Reduction of the last elements 

template <typename REAL,int blockSize>
__forceinline__ __device__ void warp_reduce(volatile REAL* sdata, int tid){
  if (blockSize>=64) sdata[tid] += sdata[tid+32];
  if (blockSize>=32) sdata[tid] += sdata[tid+16];
  if (blockSize>=16) sdata[tid] += sdata[tid+8];
  if (blockSize>=8) sdata[tid] += sdata[tid+4];
  if (blockSize>=4) sdata[tid] += sdata[tid+2];
  if (blockSize>=2) sdata[tid] += sdata[tid+1];
}

template <typename REAL,int blockSize>
__global__ void sum_reduce(REAL *vec,int ns){
  extern __shared__ REAL sdata[];
  unsigned int i = ns*blockIdx.x*blockSize + threadIdx.x;
  unsigned int tid = threadIdx.x;
  sdata[tid] = 0.0;
  for (int j=0;j<ns;++j,i+=blockSize){
    sdata[tid] += vec[i];
  }
  __syncthreads();

  if (blockSize >= 1024) {
    if (tid < 512) sdata[tid]+=sdata[tid+512];
    __syncthreads();
  }

  if (blockSize >= 512) {
    if (tid < 256) sdata[tid]+=sdata[tid+256];
    __syncthreads();
  }

  if (blockSize >= 256) {
    if (tid < 128) sdata[tid]+=sdata[tid+128];
    __syncthreads();
  }

  if (blockSize >= 128) {
    if (tid < 64) sdata[tid]+=sdata[tid+64];
    __syncthreads();
  }

  //Unroll the last warp
  if (tid < 32){
    warp_reduce<double,blockSize>(sdata,tid);
  }

  if (tid == 0) vec[blockIdx.x] = sdata[0];
}


//Kernel calling function

void global_reduction_sum(double* resField, unsigned int Npow2){
  unsigned int max_grid = (grid_size_max_gpu >> 1u)+1u; 
  int ns = Npow2/GSUM_BLOCK_SIZE/max_grid;
  int n;
  ns = (4>ns) ? 4 : ns;
  for (n=Npow2;n>GSUM_BLOCK_SIZE*ns;n/=GSUM_BLOCK_SIZE*ns){
    int grid = n/GSUM_BLOCK_SIZE/ns;
    sum_reduce<double,GSUM_BLOCK_SIZE><<<grid,BLOCK_SIZE,(size_t) BLOCK_SIZE*sizeof(double)>>>(resField,ns);
  }
  int nst = ns;
  while (nst>n) nst/=2;
  n/=nst;
  switch (n){
  case 1024:
    sum_reduce<double,1024><<<1,n,n*sizeof(double)>>>(resField,nst);
    break;
  case 512:
    sum_reduce<double,512><<<1,n,n*sizeof(double)>>>(resField,nst);
    break;
  case 256: 
    sum_reduce<double,256><<<1,n,n*sizeof(double)>>>(resField,nst);
    break;
  case 128: 
    sum_reduce<double,128><<<1,n,n*sizeof(double)>>>(resField,nst);
    break;
  case 64: 
    sum_reduce<double,64><<<1,n,n*sizeof(double)>>>(resField,nst);
    break;
  case 32: 
    sum_reduce<double,32><<<1,n,48*sizeof(double)>>>(resField,nst);
    break;
  case 16: 
    sum_reduce<double,16><<<1,n,24*sizeof(double)>>>(resField,nst);
    break;
  case 8: 
    sum_reduce<double,8><<<1,n,12*sizeof(double)>>>(resField,nst);
    break;
  case 4: 
    sum_reduce<double,4><<<1,n,6*sizeof(double)>>>(resField,nst);
    break;
  case 2: 
    sum_reduce<double,2><<<1,n,3*sizeof(double)>>>(resField,nst);
    break;
  case 1: 
    sum_reduce<double,1><<<1,n,2*sizeof(double)>>>(resField,nst);
    break;
  }
}


////////////////////////////////////////////////////// Optimized Global Sum Kernel
// Requires multiple kernel calls in order to produce global synchronization for global sum
// (You might think that this gives some overhead, but you are wrong.)


template <typename COMPLEX,int blockSize>
__forceinline__ __device__ void warp_reduce_complex(volatile COMPLEX* sdata, int tid){
  if (blockSize>=64) {_complex_add(sdata[tid],sdata[tid],sdata[tid+32]);}
  if (blockSize>=32) {_complex_add(sdata[tid],sdata[tid],sdata[tid+16]);}
  if (blockSize>=16) {_complex_add(sdata[tid],sdata[tid],sdata[tid+8]);}
  if (blockSize>=8) {_complex_add(sdata[tid],sdata[tid],sdata[tid+4]);}
  if (blockSize>=4) {_complex_add(sdata[tid],sdata[tid],sdata[tid+2]);}
  if (blockSize>=2) {_complex_add(sdata[tid],sdata[tid],sdata[tid+1]);}}


template <typename COMPLEX,int blockSize>
__global__ void sum_reduce_complex(COMPLEX *vec, int ns){
  extern __shared__ COMPLEX sdata_cmplx[];
  unsigned int i = ns*blockIdx.x*blockSize + threadIdx.x;
  unsigned int tid = threadIdx.x;
  sdata_cmplx[tid].re = 0.0;
  sdata_cmplx[tid].im = 0.0;
  for (int j=0;j<ns;++j,i+=blockSize){
    _complex_add(sdata_cmplx[tid],sdata_cmplx[tid],vec[i]);
  }
  __syncthreads();

  if (blockSize >= 1024) {
    if (tid < 512) {_complex_add(sdata_cmplx[tid],sdata_cmplx[tid],sdata_cmplx[tid+512]);}
    __syncthreads();
  }

  if (blockSize >= 512) {
    if (tid < 256) {_complex_add(sdata_cmplx[tid],sdata_cmplx[tid],sdata_cmplx[tid+256]);}
    __syncthreads();
  }

  if (blockSize >= 256) {
    if (tid < 128) {_complex_add(sdata_cmplx[tid],sdata_cmplx[tid],sdata_cmplx[tid+128]);}
    __syncthreads();
  }

  if (blockSize >= 128) {
    if (tid < 64) {_complex_add(sdata_cmplx[tid],sdata_cmplx[tid],sdata_cmplx[tid+64]);}
    __syncthreads();
  }

  //Unroll the last warp
  if (tid < 32){
    warp_reduce_complex<hr_complex,blockSize>(sdata_cmplx,tid);
  }

  if (tid == 0) vec[blockIdx.x] = sdata_cmplx[0];
}


void global_reduction_complex_sum(hr_complex* resField, unsigned int Npow2){
  unsigned int max_grid = (grid_size_max_gpu >> 1u)+1u; 
  int ns = Npow2/GSUM_BLOCK_SIZE/max_grid;
  int n;
  ns = (4>ns) ? 4 : ns;
  for (n=Npow2;n>GSUM_BLOCK_SIZE*ns;n/=GSUM_BLOCK_SIZE*ns){
    int grid = n/GSUM_BLOCK_SIZE/ns;
    sum_reduce_complex<hr_complex,GSUM_BLOCK_SIZE><<<grid,BLOCK_SIZE,(size_t) BLOCK_SIZE*sizeof(hr_complex)>>>(resField,ns);
  }
  int nst = ns;
  while (nst>n) nst/=2;
  n/=nst;
  switch (n){
  case 1024:
    sum_reduce_complex<hr_complex,1024><<<1,n,n*sizeof(hr_complex)>>>(resField,nst);
    break;
  case 512:
    sum_reduce_complex<hr_complex,512><<<1,n,n*sizeof(hr_complex)>>>(resField,nst);
    break;
  case 256: 
    sum_reduce_complex<hr_complex,256><<<1,n,n*sizeof(hr_complex)>>>(resField,nst);
    break;
  case 128: 
    sum_reduce_complex<hr_complex,128><<<1,n,n*sizeof(hr_complex)>>>(resField,nst);
    break;
  case 64: 
    sum_reduce_complex<hr_complex,64><<<1,n,n*sizeof(hr_complex)>>>(resField,nst);
    break;
  case 32: 
    sum_reduce_complex<hr_complex,32><<<1,n,48*sizeof(hr_complex)>>>(resField,nst);
    break;
  case 16: 
    sum_reduce_complex<hr_complex,16><<<1,n,24*sizeof(hr_complex)>>>(resField,nst);
    break;
  case 8: 
    sum_reduce_complex<hr_complex,8><<<1,n,12*sizeof(hr_complex)>>>(resField,nst);
    break;
  case 4: 
    sum_reduce_complex<hr_complex,4><<<1,n,6*sizeof(hr_complex)>>>(resField,nst);
    break;
  case 2: 
    sum_reduce_complex<hr_complex,2><<<1,n,3*sizeof(hr_complex)>>>(resField,nst);
    break;
  case 1: 
    sum_reduce_complex<hr_complex,1><<<1,n,2*sizeof(hr_complex)>>>(resField,nst);
    break;
  }
}

double global_sum_gpu(double* vector, int n){
  unsigned int new_n = next_pow2(n);
  int grid_size = new_n / BLOCK_SIZE;
  double* padded_vector;
  double res;
  padded_vector=alloc_double_sum_field(new_n);
  copy_with_zero_padding<<<grid_size,BLOCK_SIZE>>>(padded_vector, vector, n, new_n);
  global_reduction_sum(padded_vector,new_n);
  cudaMemcpy(&res,padded_vector,sizeof(res),cudaMemcpyDeviceToHost);
  return res;
}

hr_complex global_sum_gpu(hr_complex* vector, int n){
  unsigned int new_n = next_pow2(n);
  int grid_size = new_n / BLOCK_SIZE;
  hr_complex* padded_vector;
  hr_complex res;
  padded_vector=alloc_complex_sum_field(new_n);
  copy_with_zero_padding_complex<<<grid_size,BLOCK_SIZE>>>(padded_vector, vector, n, new_n);
  global_reduction_complex_sum(padded_vector,new_n);
  cudaMemcpy(&res,padded_vector,sizeof(res),cudaMemcpyDeviceToHost);
  return res;
}

#endif
