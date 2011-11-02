#ifndef LINEAR_ALGEBRA_GPU_CU
#define LINEAR_ALGEBRA_GPU_CU

#include <cuda.h>
#include "global.h"


template<typename REAL>
__device__ void warpReduce(volatile REAL* sdata, int tid){
  sdata[tid]+=sdata[tid+32];
  sdata[tid]+=sdata[tid+16];
  sdata[tid]+=sdata[tid+8];
  sdata[tid]+=sdata[tid+4];
  sdata[tid]+=sdata[tid+2];
  sdata[tid]+=sdata[tid+1];
}

//Global sum calculation
template<typename REAL>
__global__ void global_sum_gpu(REAL* in, REAL* out, int N){
  __shared__ REAL sdata[BLOCK_SIZE];
  int tid = threadIdx.x;
  int i;
  //  REAL res;
  sdata[tid]=in[tid];
  for (i=tid+BLOCK_SIZE;i<N;i+= BLOCK_SIZE){ // Sum over all blocks 
    sdata[tid]+=in[i];
  }
  __syncthreads();
  
    for (i = BLOCK_SIZE/2;i>32;i/=2){
    if (tid < i) sdata[tid]+=sdata[tid+i];
    __syncthreads();
  }
  if (tid<32){ //Unroll all the threads in a WARP
    warpReduce(sdata,tid);
  }
  __syncthreads();
  if (tid == 0)  out[0] = sdata[0];
  /*  if (tid == 0) {
    res =0;
    for (i=0;i<BLOCK_SIZE;++i){
      res+=sdata[i];
    }
    out[0] = res;
    }*/
}


/* Im <s1,s2> */
template<typename SPINOR_TYPE, typename REAL>
  __global__ void spinor_field_prod_im_gpu(SPINOR_TYPE* s1, SPINOR_TYPE* s2, REAL* resField,int N){
  int i = blockIdx.x*BLOCK_SIZE + threadIdx.x;
  if (i<N){
    _spinor_prod_im_f(resField[i],s1[i],s2[i]);
  }
}

/* s1+=r*s2 r real */
template< typename SPINOR_TYPE, typename REAL >
__global__ void spinor_field_mul_add_assign_gpu(SPINOR_TYPE *s1, REAL r, SPINOR_TYPE *s2,int N){
  int i = blockIdx.x*BLOCK_SIZE + threadIdx.x;
  if (i<N) {
    _spinor_mul_add_assign_f(s1[i],r,s2[i]);
  }
}


#endif
