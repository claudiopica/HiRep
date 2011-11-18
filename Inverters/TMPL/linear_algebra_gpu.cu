#ifndef LINEAR_ALGEBRA_GPU_CU
#define LINEAR_ALGEBRA_GPU_CU

#include "global.h"
#include "gpu.h"



//reduction of the last elements
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
}

//Reduction of the last elements
template<typename COMPLEX>
__device__ void warpReduce_complex(volatile COMPLEX* sdata, int tid){
  _complex_add(sdata[tid],sdata[tid],sdata[tid+32]);
  _complex_add(sdata[tid],sdata[tid],sdata[tid+16]);
  _complex_add(sdata[tid],sdata[tid],sdata[tid+8]);
  _complex_add(sdata[tid],sdata[tid],sdata[tid+4]);
  _complex_add(sdata[tid],sdata[tid],sdata[tid+2]);
  _complex_add(sdata[tid],sdata[tid],sdata[tid+1]);
}

//Global sum calculation
template<typename COMPLEX>
__global__ void global_sum_complex_gpu(COMPLEX* in, COMPLEX* out, int N){
  __shared__ COMPLEX sdata[BLOCK_SIZE];
  int tid = threadIdx.x;
  int i;
  //  REAL res;
  sdata[tid]=in[tid];
  for (i=tid+BLOCK_SIZE;i<N;i+= BLOCK_SIZE){ // Sum over all blocks 
    _complex_add(sdata[tid],sdata[tid],in[i]);    
  }
  __syncthreads();
  
  for (i = BLOCK_SIZE/2;i>32;i/=2){
    if (tid < i){
      _complex_add(sdata[tid],sdata[tid],sdata[tid+i]);
    }
    __syncthreads();
  }
  if (tid<32){ //Unroll all the threads in a WARP
    warpReduce_complex(sdata,tid);
  }
  __syncthreads();
  if (tid == 0)  out[0] = sdata[0];
}

/* Re <s1,s2> */
template<typename COMPLEX, typename REAL>
  __global__ void spinor_field_prod_re_gpu(COMPLEX* s1, COMPLEX* s2, REAL* resField,int N){
  int i = blockIdx.x*BLOCK_SIZE + threadIdx.x;
  i=min(i,N-1);
  resField[i] = _complex_prod_re(s1[i],s2[i]);
}

/* Im <s1,s2> */
template<typename COMPLEX, typename REAL>
  __global__ void spinor_field_prod_im_gpu(COMPLEX* s1, COMPLEX* s2, REAL* resField,int N){
  int i = blockIdx.x*BLOCK_SIZE + threadIdx.x;
  i=min(i,N-1);
  resField[i]=_complex_prod_im(s1[i],s2[i]);
}

/* <s1,s2> */
template< typename COMPLEX>
  __global__ void spinor_field_prod_gpu(COMPLEX* s1, COMPLEX* s2, COMPLEX* resField,int N){
  int i = blockIdx.x*BLOCK_SIZE + threadIdx.x;
  i=min(i,N-1);
  _complex_prod(resField[i],s1[i],s2[i]);
}

/* Re <g5*s1,s2> */
template<typename COMPLEX, typename REAL>
  __global__ void spinor_field_g5_prod_re_gpu(COMPLEX* s1, COMPLEX* s2, REAL* resField,int N){
  int i = blockIdx.x*BLOCK_SIZE + threadIdx.x;
  i=min(i,N-1);
  resField[i]=_complex_prod_re(s1[i],s2[i]);
  if (i>((N>>1)-1)){
  	resField[i]=-resField[i];
  }
}

/* Im <g5*s1,s2> */
template<typename COMPLEX, typename REAL>
  __global__ void spinor_field_g5_prod_im_gpu(COMPLEX* s1, COMPLEX* s2, REAL* resField,int N){
  int i = blockIdx.x*BLOCK_SIZE + threadIdx.x;
  i=min(i,N-1);
  resField[i]=_complex_prod_im(s1[i],s2[i]);
  if (i>((N>>1)-1)){
  	resField[i]=-resField[i];
  }
}

/* Re <s1,s1> */ 
template<typename COMPLEX, typename REAL>
  __global__ void spinor_field_sqnorm_gpu(COMPLEX* s1, REAL* resField,int N){
  int i = blockIdx.x*BLOCK_SIZE + threadIdx.x;
  i=min(i,N-1);
  resField[i] = _complex_prod_re(s1[i],s1[i]);
}

/* s1+=r*s2 r real */
template< typename COMPLEX, typename REAL >
__global__ void spinor_field_mul_add_assign_gpu(COMPLEX *s1, REAL r, COMPLEX *s2,int N){
  int i = blockIdx.x*BLOCK_SIZE + threadIdx.x;
  i=min(i,N-1);
  _complex_mulr_assign(s1[i],r,s2[i]);
}

/* s1=r*s2 r real */
template< typename COMPLEX , typename REAL >
__global__ void spinor_field_mul_gpu(COMPLEX *s1, REAL r, COMPLEX *s2,int N){
  int i = blockIdx.x*BLOCK_SIZE + threadIdx.x;
  i=min(i,N-1);
  _complex_mulr(s1[i],r,s2[i]);

}

/* s1+=c*s2 c complex */
template< typename COMPLEX >
__global__ void spinor_field_mulc_add_assign_gpu(COMPLEX *s1, COMPLEX c, COMPLEX *s2,int N){
  int i = blockIdx.x*BLOCK_SIZE + threadIdx.x;
  i=min(i,N-1);
  _complex_mul_assign(s1[i],c,s2[i]);
}

/* s1=c*s2 c complex */
template< typename COMPLEX >
__global__ void spinor_field_mulc_gpu(COMPLEX *s1, COMPLEX c, COMPLEX *s2,int N){
  int i = blockIdx.x*BLOCK_SIZE + threadIdx.x;
  i=min(i,N-1);
  _complex_mul(s1[i],c,s2[i]);
}

/* r=s1+s2 */
template< typename COMPLEX>
__global__ void spinor_field_add_gpu(COMPLEX *r, COMPLEX *s1, COMPLEX *s2,int N){
  int i = blockIdx.x*BLOCK_SIZE + threadIdx.x;
  i=min(i,N-1);
  _complex_add(r[i],s1[i],s2[i]);

}

/* r=s1-s2 */
template< typename COMPLEX >
__global__ void spinor_field_sub_gpu(COMPLEX *r, COMPLEX * s1, COMPLEX *s2,int N){
  int i = blockIdx.x*BLOCK_SIZE + threadIdx.x;
  i=min(i,N-1);
  _complex_sub(r[i],s1[i],s2[i]);
}

/* s1+=s2 */
template< typename COMPLEX>
__global__ void spinor_field_add_assign_gpu(COMPLEX *s1, COMPLEX *s2,int N){
  int i = blockIdx.x*BLOCK_SIZE + threadIdx.x;
  i=min(i,N-1);
  _complex_add_assign(s1[i],s2[i]);
}

/* s1-=s2 */
template< typename COMPLEX >
__global__ void spinor_field_sub_assign_gpu(COMPLEX* s1, COMPLEX *s2,int N){
  int i = blockIdx.x*BLOCK_SIZE + threadIdx.x;
  i=min(i,N-1);
  _complex_sub_assign(s1[i],s2[i]);
}

/* s1=0 */
template< typename COMPLEX>
__global__ void spinor_field_zero_gpu(COMPLEX *s1,int N){
  int i = blockIdx.x*BLOCK_SIZE + threadIdx.x;
  i=min(i,N-1);
  _complex_0(s1[i]);
}

/* s1=-s2 */
template< typename COMPLEX >
__global__ void spinor_field_minus_gpu(COMPLEX* s1, COMPLEX *s2,int N){
  int i = blockIdx.x*BLOCK_SIZE + threadIdx.x;
  i=min(i,N-1);
  _complex_minus(s1[i],s2[i]);
}

/* s1=r1*s2+r2*s3 */
template< typename COMPLEX , typename REAL >
__global__ void spinor_field_lc_gpu(COMPLEX *s1, REAL r1, COMPLEX *s2, REAL r2, COMPLEX *s3, int N){
  int i = blockIdx.x*BLOCK_SIZE + threadIdx.x;
  i=min(i,N-1);
  _complex_rlc(s1[i],r1,s2[i],r2,s3[i]);
}

/* s1+=r*s2 r real */
template< typename COMPLEX, typename REAL >
__global__ void spinor_field_lc_add_assign_gpu(COMPLEX *s1, REAL r1, COMPLEX *s2, REAL r2, COMPLEX *s3, int N){
  int i = blockIdx.x*BLOCK_SIZE + threadIdx.x;
  i=min(i,N-1);
  _complex_rlc_assign(s1[i],r1,s2[i],r2,s3[i]);
}


/* s1=cd1*s2+cd2*s3 cd1, cd2 complex*/
template< typename COMPLEX >
__global__ void spinor_field_clc_gpu(COMPLEX *s1, COMPLEX c1, COMPLEX *s2, COMPLEX c2, COMPLEX *s3, int N){
  int i = blockIdx.x*BLOCK_SIZE + threadIdx.x;
  i=min(i,N-1);
  _complex_clc(s1[i],c1,s2[i],c2,s3[i]);
}

/* s1+=r*s2 r real */
template< typename COMPLEX >
__global__ void spinor_field_clc_add_assign_gpu(COMPLEX *s1, COMPLEX c1, COMPLEX *s2, COMPLEX c2, COMPLEX *s3, int N){
  int i = blockIdx.x*BLOCK_SIZE + threadIdx.x;
  i=min(i,N-1);
  _complex_clc_assign(s1[i],c1,s2[i],c2,s3[i]);
}

/* s1=g5*s2  */
template< typename COMPLEX>
__global__ void spinor_field_g5_gpu(COMPLEX *s1, COMPLEX *s2,int N){
  int i = blockIdx.x*BLOCK_SIZE + threadIdx.x;
  i=min(i,N-1);
  if ( i < (N>>1) ) {
    s1[i]=s2[i];
  }
  else{
    _complex_minus(s1[i],s2[i]);
  }
}

/* s1=g5*s1 */
template< typename COMPLEX >
__global__ void spinor_field_g5_assign_gpu(COMPLEX* s1,int N){
  int i = blockIdx.x*BLOCK_SIZE + threadIdx.x;
  i=min(i,N-1);
  if (i>((N>>1)-1)){
    _complex_minus(s1[i],s1[i]);
  }
}

/* tools per eva.c  */
template< typename COMPLEX , typename REAL >
__global__ void spinor_field_lc1_gpu(REAL r, COMPLEX *s1, COMPLEX *s2, int N){
  int i = blockIdx.x*BLOCK_SIZE + threadIdx.x;
  i=min(i,N-1);
  _complex_mulr_assign(s1[i],r,s2[i]);
}


template< typename COMPLEX, typename REAL >
__global__ void spinor_field_lc2_gpu(REAL r1, REAL r2, COMPLEX *s1, COMPLEX *s2, int N){
  int i = blockIdx.x*BLOCK_SIZE + threadIdx.x;
  i=min(i,N-1);
  _complex_rlc(s1[i],r1,s1[i],r2,s2[i]);
}

template< typename COMPLEX , typename REAL >
__global__ void spinor_field_lc3_gpu(REAL r1,REAL r2, COMPLEX *s1, COMPLEX *s2, COMPLEX *s3, int N){
  int i = blockIdx.x*BLOCK_SIZE + threadIdx.x;
  i=min(i,N-1);
  _complex_rlc_assign(s3[i],r1,s1[i],r2,s2[i]);
  _complex_minus(s3[i],s3[i]);
}

/* c1=0 */
/*
template< typename COMPLEX>
__global__ void complex_field_zero_gpu(COMPLEX *c1,int N){
  int i = blockIdx.x*BLOCK_SIZE + threadIdx.x;
  if (i<N) {
    c1[i].re=0;
    c1[i].im=0;
  }
}
*/


#endif
