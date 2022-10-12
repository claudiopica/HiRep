#ifndef LINEAR_ALGEBRA_GPU_CU
#define LINEAR_ALGEBRA_GPU_CU

extern "C" {
#include "global.h"
#include "gpu.h"
}

/* Re <s1,s2> */
template<typename COMPLEX, typename REAL>
__global__ void spinor_field_prod_re_gpu(COMPLEX* s1, COMPLEX* s2, REAL* resField, int N)
{
  for (int i = blockIdx.x * blockDim.x + threadIdx.x;
       i < N;
       i += blockDim.x * gridDim.x)
    {
       resField[i] = _complex_prod_re(s1[i],s2[i]);
    }
}

/* Im <s1,s2> */
template<typename COMPLEX, typename REAL>
__global__ void spinor_field_prod_im_gpu(COMPLEX* s1, COMPLEX* s2, REAL* resField, int N)
{
  for (int i = blockIdx.x * blockDim.x + threadIdx.x;
       i < N;
       i += blockDim.x * gridDim.x)
    {
       resField[i]=_complex_prod_im(s1[i],s2[i]);
    }
}

/* <s1,s2> */
template< typename COMPLEX>
__global__ void spinor_field_prod_gpu(COMPLEX* s1, COMPLEX* s2, hr_complex* resField, int N)
{
  for (int i = blockIdx.x * blockDim.x + threadIdx.x;
       i < N;
       i += blockDim.x * gridDim.x)
    {
       resField[i] = _complex_prod(s1[i],s2[i]);
    }
}

/* Re <g5*s1,s2> */
template<typename COMPLEX, typename REAL>
__global__ void spinor_field_g5_prod_re_gpu(COMPLEX* s1, COMPLEX* s2, REAL* resField, int N)
{
  for (int i = blockIdx.x * blockDim.x + threadIdx.x;
       i < N;
       i += blockDim.x * gridDim.x)
    {
       resField[i] = - _complex_prod_re(s1[i],s2[i]);
    }
}

/* Im <g5*s1,s2> */
template<typename COMPLEX, typename REAL>
 __global__ void spinor_field_g5_prod_im_gpu(COMPLEX* s1, COMPLEX* s2, REAL* resField, int N)
{
  for (int i = blockIdx.x * blockDim.x + threadIdx.x;
       i < N;
       i += blockDim.x * gridDim.x)
    {
       resField[i] = - _complex_prod_im(s1[i],s2[i]);
    }
}

/* Re <s1,s1> */
template<typename COMPLEX, typename REAL>
__global__ void spinor_field_sqnorm_gpu(COMPLEX* s1, REAL* resField, int N)
{
  for (int i = blockIdx.x * blockDim.x + threadIdx.x;
       i < N;
       i += blockDim.x * gridDim.x)
    {
       resField[i] = _complex_prod_re(s1[i],s1[i]);
    }
}

/* s1+=r*s2 r real */
template< typename COMPLEX, typename REAL >
__global__ void spinor_field_mul_add_assign_gpu(COMPLEX *s1, REAL r, COMPLEX *s2, int N)
{
  for (int i = blockIdx.x * blockDim.x + threadIdx.x;
       i < N;
       i += blockDim.x * gridDim.x)
    {
       _complex_mulr_assign(s1[i],r,s2[i]);
    }
}

/* s1=r*s2 r real */
template< typename COMPLEX , typename REAL >
__global__ void spinor_field_mul_gpu(COMPLEX *s1, REAL r, COMPLEX *s2, int N)
{
  for (int i = blockIdx.x * blockDim.x + threadIdx.x;
       i < N;
       i += blockDim.x * gridDim.x)
    {
       _complex_mulr(s1[i],r,s2[i]);
    }
}

/* s1+=c*s2 c complex */
template< typename COMPLEX >
__global__ void spinor_field_mulc_add_assign_gpu(COMPLEX *s1, COMPLEX c, COMPLEX *s2, int N)
{
  for (int i = blockIdx.x * blockDim.x + threadIdx.x;
       i < N;
       i += blockDim.x * gridDim.x)
    {
       _complex_mul_assign(s1[i],c,s2[i]);
    }
}

/* s1=c*s2 c complex */
template< typename COMPLEX >
__global__ void spinor_field_mulc_gpu(COMPLEX *s1, COMPLEX c, COMPLEX *s2, int N)
{
  for (int i = blockIdx.x * blockDim.x + threadIdx.x;
       i < N;
       i += blockDim.x * gridDim.x)
    {
       _complex_mul(s1[i],c,s2[i]);
    }
}

/* r=s1+s2 */
template< typename COMPLEX>
__global__ void spinor_field_add_gpu(COMPLEX *r, COMPLEX *s1, COMPLEX *s2, int N)
{
  for (int i = blockIdx.x * blockDim.x + threadIdx.x;
       i < N;
       i += blockDim.x * gridDim.x)
    {
       _complex_add(r[i],s1[i],s2[i]);
    }
}

/* r=s1-s2 */
template< typename COMPLEX >
__global__ void spinor_field_sub_gpu(COMPLEX *r, COMPLEX * s1, COMPLEX *s2, int N)
{
  for (int i = blockIdx.x * blockDim.x + threadIdx.x;
       i < N;
       i += blockDim.x * gridDim.x)
    {
       _complex_sub(r[i],s1[i],s2[i]);
    }
}

/* s1+=s2 */
template< typename COMPLEX>
__global__ void spinor_field_add_assign_gpu(COMPLEX *s1, COMPLEX *s2, int N)
{
  for (int i = blockIdx.x * blockDim.x + threadIdx.x;
       i < N;
       i += blockDim.x * gridDim.x)
    {
       _complex_add_assign(s1[i],s2[i]);
    }
}

/* s1-=s2 */
template< typename COMPLEX >
__global__ void spinor_field_sub_assign_gpu(COMPLEX* s1, COMPLEX *s2, int N)
{
  for (int i = blockIdx.x * blockDim.x + threadIdx.x;
       i < N;
       i += blockDim.x * gridDim.x)
    {
       _complex_sub_assign(s1[i],s2[i]);
    }
}

/* s1=0 */
template< typename COMPLEX>
__global__ void spinor_field_zero_gpu(COMPLEX *s1, int N)
{
  for (int i = blockIdx.x * blockDim.x + threadIdx.x;
       i < N;
       i += blockDim.x * gridDim.x)
    {
       _complex_0(s1[i]);
    }
}

/* s1=-s2 */
template< typename COMPLEX >
__global__ void spinor_field_minus_gpu(COMPLEX * s1, COMPLEX * s2, int N)
{
  for (int i = blockIdx.x * blockDim.x + threadIdx.x;
       i < N;
       i += blockDim.x * gridDim.x)
    {
       _complex_minus(s1[i],s2[i]);
    }
}

/* s1=-s1 */
template< typename COMPLEX >
__global__ void spinor_field_minus_assign_gpu(COMPLEX* s1, int N)
{
  for (int i = blockIdx.x * blockDim.x + threadIdx.x;
       i < N;
       i += blockDim.x * gridDim.x)
    {
       s1[i] = -s1[i];
    }
}

/* s1=r1*s2+r2*s3 r1,r2 real */
template< typename COMPLEX , typename REAL >
__global__ void spinor_field_lc_gpu(COMPLEX *s1, REAL r1, COMPLEX *s2, REAL r2, COMPLEX *s3, int N)
{
  for (int i = blockIdx.x * blockDim.x + threadIdx.x;
       i < N;
       i += blockDim.x * gridDim.x)
    {
       _complex_rlc(s1[i],r1,s2[i],r2,s3[i]);
    }
}

/* s1+=r1*s2+r2*s3 r1,r2 real */
template< typename COMPLEX, typename REAL >
__global__ void spinor_field_lc_add_assign_gpu(COMPLEX *s1, REAL r1, COMPLEX *s2, REAL r2, COMPLEX *s3, int N)
{
  for (int i = blockIdx.x * blockDim.x + threadIdx.x;
       i < N;
       i += blockDim.x * gridDim.x)
    {
       _complex_rlc_assign(s1[i],r1,s2[i],r2,s3[i]);
    }
}

/* s1=cd1*s2+cd2*s3 cd1, cd2 complex */
template< typename COMPLEX >
__global__ void spinor_field_clc_gpu(COMPLEX *s1, COMPLEX c1, COMPLEX *s2, COMPLEX c2, COMPLEX *s3, int N)
{
  for (int i = blockIdx.x * blockDim.x + threadIdx.x;
       i < N;
       i += blockDim.x * gridDim.x)
    {
       _complex_clc(s1[i],c1,s2[i],c2,s3[i]);
    }
}

/* s1+=cd1*s2+cd2*s3 cd1,cd2 complex */
template< typename COMPLEX >
__global__ void spinor_field_clc_add_assign_gpu(COMPLEX *s1, COMPLEX c1, COMPLEX *s2, COMPLEX c2, COMPLEX *s3, int N)
{
  for (int i = blockIdx.x * blockDim.x + threadIdx.x;
       i < N;
       i += blockDim.x * gridDim.x)
    {
       _complex_clc_assign(s1[i],c1,s2[i],c2,s3[i]);
    }
}

/* s1=g5*s2  */
template< typename COMPLEX>
__global__ void spinor_field_g5_gpu(COMPLEX *s1, COMPLEX *s2,int N)
{
  for (int i = blockIdx.x * blockDim.x + threadIdx.x;
       i < N;
       i += blockDim.x * gridDim.x)
    {
       _complex_minus(s1[i],s2[i]);
    }
}

/* s1=g5*s1 */
template< typename COMPLEX >
__global__ void spinor_field_g5_assign_gpu(COMPLEX* s1,int N)
{
  for (int i = blockIdx.x * blockDim.x + threadIdx.x;
       i < N;
       i += blockDim.x * gridDim.x)
    {
       _complex_minus(s1[i],s1[i]);
    }
}

/* s1+=c*g5*s2 c complex  */
template< typename COMPLEX>
__global__ void spinor_field_g5_mulc_add_assign_gpu(COMPLEX *s1, COMPLEX c, COMPLEX *s2,int N)
{
  for (int i = blockIdx.x * blockDim.x + threadIdx.x;
       i < N;
       i += blockDim.x * gridDim.x)
    {
       _complex_mul_assign(s1[i],-c,s2[i]);
    }
}

/* tools per eva.c  */
template< typename COMPLEX , typename REAL >
__global__ void spinor_field_lc1_gpu(REAL r, COMPLEX *s1, COMPLEX *s2, int N)
{
  for (int i = blockIdx.x * blockDim.x + threadIdx.x;
       i < N;
       i += blockDim.x * gridDim.x)
    {
       _complex_mulr_assign(s1[i],r,s2[i]);
    }
}

template< typename COMPLEX, typename REAL >
__global__ void spinor_field_lc2_gpu(REAL r1, REAL r2, COMPLEX *s1, COMPLEX *s2, int N)
{
  for (int i = blockIdx.x * blockDim.x + threadIdx.x;
       i < N;
       i += blockDim.x * gridDim.x)
    {
       _complex_rlc(s1[i],r1,s1[i],r2,s2[i]);
    }
}

template< typename COMPLEX , typename REAL >
__global__ void spinor_field_lc3_gpu(REAL r1,REAL r2, COMPLEX *s1, COMPLEX *s2, COMPLEX *s3, int N)
{
  for (int i = blockIdx.x * blockDim.x + threadIdx.x;
       i < N;
       i += blockDim.x * gridDim.x)
    {
       _complex_rlc_assign(s3[i],r1,s1[i],r2,s2[i]);
       _complex_minus(s3[i],s3[i]);
    }
}

#endif
