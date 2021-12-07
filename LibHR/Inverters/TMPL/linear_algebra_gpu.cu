#ifndef LINEAR_ALGEBRA_GPU_CU
#define LINEAR_ALGEBRA_GPU_CU

extern "C" {
#include "global.h"
#include "gpu.h"
}

/* Re <s1,s2>  FOR OPTIMIZED GLOBAL SUM*/
template<typename COMPLEX, typename REAL>
  __global__ void spinor_field_prod_re_padded_gpu(COMPLEX* s1, COMPLEX* s2, REAL* resField,unsigned int N){
  int i = blockIdx.x*BLOCK_SIZE + threadIdx.x;
  COMPLEX c1 = s1[i];
  COMPLEX c2 = s2[i];
  int i2 = i+gridDim.x*BLOCK_SIZE;
  REAL tmp1,tmp2;
  if (i2<N){
    COMPLEX c3 = s1[i2];
    COMPLEX c4 = s2[i2];
    tmp1 = _complex_prod_re(c1,c2);
    tmp2 = _complex_prod_re(c3,c4);
  }
  else{
    tmp1 = _complex_prod_re(c1,c2);
    tmp2 = 0.0;
  }
  resField[i] = tmp1+tmp2;
}

/* Im <s1,s2> */
template<typename COMPLEX, typename REAL>
  __global__ void spinor_field_prod_im_padded_gpu(COMPLEX* s1, COMPLEX* s2, REAL* resField,unsigned int N){
  int i = blockIdx.x*BLOCK_SIZE + threadIdx.x;
  COMPLEX c1 = s1[i];
  COMPLEX c2 = s2[i];
  int i2 = i+gridDim.x*BLOCK_SIZE;
  int i3 = min(i2,N-1);
  COMPLEX c3 = s1[i3];
  COMPLEX c4 = s2[i3];
  REAL tmp1 = _complex_prod_im(c1,c2);
  REAL tmp2 = (i2<N) ? _complex_prod_im(c3,c4) : 0;
  resField[i] = tmp1+tmp2;
}

/* <s1,s2> */
template< typename COMPLEX>
__global__ void spinor_field_prod_padded_gpu(COMPLEX* s1, COMPLEX* s2, hr_complex* resField,unsigned int N){
  int i = blockIdx.x*BLOCK_SIZE + threadIdx.x;
  COMPLEX c1 = s1[i];
  COMPLEX c2 = s2[i];
  int i2 = i+gridDim.x*BLOCK_SIZE;
  COMPLEX tmp1,tmp2;
  if (i2<N) {
    COMPLEX c3 = s1[i2];
    COMPLEX c4 = s2[i2];
    _complex_prod(tmp1,c1,c2);
    _complex_prod(tmp2,c3,c4);
  }
  else{
    _complex_prod(tmp1,c1,c2);
    tmp2.re=tmp2.im=0.;
  }
  resField[i].re=tmp1.re+tmp2.re;
  resField[i].im=tmp1.im+tmp2.im;
}

/* Re <g5*s1,s2> */
/*template<typename COMPLEX, typename REAL>
  __global__ void spinor_field_g5_prod_re_padded_gpu(COMPLEX* s1, COMPLEX* s2, REAL* resField,unsigned int N){
  int i = blockIdx.x*BLOCK_SIZE + threadIdx.x;
  COMPLEX c1 = s1[i];
  COMPLEX c2 = s2[i];
  int i2 = i+gridDim.x*BLOCK_SIZE;
  REAL tmp1,tmp2;
  if (i2<N){
    COMPLEX c3 = s1[i2];
    COMPLEX c4 = s2[i2];
    tmp1 = _complex_prod_re(c1,c2);
    //    tmp1 = (i<(N>>1)) ? _complex_prod_re(c1,c2) : -_complex_prod_re(c1,c2);
    tmp2 = -_complex_prod_re(c3,c4);
  }
  else{
    tmp1 = (i<(N>>1)) ? _complex_prod_re(c1,c2) : -_complex_prod_re(c1,c2);
    tmp2 = 0.0;
  }
  resField[i]= tmp1+tmp2;
  }*/

template<typename COMPLEX, typename REAL>
  __global__ void spinor_field_minus_prod_re_padded_gpu(COMPLEX* s1, COMPLEX* s2, REAL* resField,unsigned int N){
  int i = blockIdx.x*BLOCK_SIZE + threadIdx.x;
  COMPLEX c1 = s1[i];
  COMPLEX c2 = s2[i];
  int i2 = i+gridDim.x*BLOCK_SIZE;
  REAL tmp1,tmp2;
  if (i2<N){
    COMPLEX c3 = s1[i2];
    COMPLEX c4 = s2[i2];
    tmp1 = _complex_prod_re(c1,c2);
    tmp2 = _complex_prod_re(c3,c4);
  }
  else{
    tmp1 = _complex_prod_re(c1,c2);
    tmp2 = 0.0;
  }
  resField[i] = -(tmp1+tmp2);
}


/* Im <g5*s1,s2> */
/*template<typename COMPLEX, typename REAL>
  __global__ void spinor_field_g5_prod_im_padded_gpu(COMPLEX* s1, COMPLEX* s2, REAL* resField,unsigned int N){
  int i = blockIdx.x*BLOCK_SIZE + threadIdx.x;
  COMPLEX c1 = s1[i];
  COMPLEX c2 = s2[i];
  int i2 = i+gridDim.x*BLOCK_SIZE;
  REAL tmp1,tmp2;
  if (i2<N){
    COMPLEX c3 = s1[i2];
    COMPLEX c4 = s2[i2];
    tmp1 = _complex_prod_im(c1,c2);
    //    tmp1 = (i<(N>>1)) ? _complex_prod_re(c1,c2) : -_complex_prod_re(c1,c2);
    tmp2 = -_complex_prod_im(c3,c4);
  }
  else{
    tmp1 = (i<(N>>1)) ? _complex_prod_im(c1,c2) : -_complex_prod_im(c1,c2);
    tmp2 = 0.0;
  }
  resField[i]= tmp1+tmp2;
  }*/

template<typename COMPLEX, typename REAL>
  __global__ void spinor_field_minus_prod_im_padded_gpu(COMPLEX* s1, COMPLEX* s2, REAL* resField,unsigned int N){
  int i = blockIdx.x*BLOCK_SIZE + threadIdx.x;
  COMPLEX c1 = s1[i];
  COMPLEX c2 = s2[i];
  int i2 = i+gridDim.x*BLOCK_SIZE;
  REAL tmp1,tmp2;
  if (i2<N){
    COMPLEX c3 = s1[i2];
    COMPLEX c4 = s2[i2];
    tmp1 = _complex_prod_im(c1,c2);
    tmp2 = _complex_prod_im(c3,c4);
  }
  else{
    tmp1 = _complex_prod_im(c1,c2);
    tmp2 = 0.0;
  }
  resField[i] = -(tmp1+tmp2);
}

/* Re <s1,s1> */ 
template<typename COMPLEX, typename REAL>
  __global__ void spinor_field_sqnorm_padded_gpu(COMPLEX* s1, REAL* resField,unsigned int N){
  int i = blockIdx.x*BLOCK_SIZE + threadIdx.x;
  COMPLEX c1 = s1[i];
  int i2 = i+gridDim.x*BLOCK_SIZE;
  REAL tmp1,tmp2;
  if (i2<N){
    COMPLEX c2 = s1[i2];
    tmp1 = _complex_prod_re(c1,c1);
    tmp2 = _complex_prod_re(c2,c2);
  }
  else{
    tmp1 = _complex_prod_re(c1,c1);
    tmp2 = 0.0;
  }
  resField[i] = tmp1+tmp2;
}

/////////// *********************



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

/* s1=-s1 */
template< typename COMPLEX >
__global__ void spinor_field_minus_assign_gpu(COMPLEX* s1,int N){
  int i = blockIdx.x*BLOCK_SIZE + threadIdx.x;
  if  (i<N){
    s1[i].re = -s1[i].re;
    s1[i].im = -s1[i].im;
  }

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

template< typename COMPLEX>
__global__ void spinor_field_copy_gpu_to_gpu_gpu(COMPLEX* dst, COMPLEX* src, int N){
  int i = blockIdx.x*BLOCK_SIZE + threadIdx.x;
  i=min(i,N-1);
  dst[i]=src[i];
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
