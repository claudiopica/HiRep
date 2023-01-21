#ifndef LINEAR_ALGEBRA_GPU_CU
#define LINEAR_ALGEBRA_GPU_CU

#define _CONCAT(_name,_suffix) _name ## _suffix
#define _F_NAME(_name,_suffix) _CONCAT(_name,_suffix)

#define spinor_idx(_idx) (_idx/(4*NF*THREADSIZE))*THREADSIZE + _idx % 4*NF*THREADSIZE % THREADSIZE

// Loop over current kernel thread block
#define CUDA_BLK_FOR(_i,_N) \
   for (int _i = blockIdx.x * blockDim.x + threadIdx.x; \
        (_i < ((_N-1)/THREADSIZE + 1)*THREADSIZE) && (spinor_idx(_i) < _N/(4*NF));\
        _i += blockDim.x * gridDim.x)

#ifdef FIXED_STRIDE
   #define _is_in_second_hspinor(_idx, _piece_size) (_idx/(2*NF*THREADSIZE)) % 2 == 1
#else
   #define _is_in_second_hspinor(_idx, _piece_size) _idx > (_piece_size/2)
#endif

// These functions do not consider the case of real reps (SAM)
template<typename SITE_TYPE, typename COMPLEX>
__device__ void read_gpu(int stride, SITE_TYPE *s, SITE_TYPE *in, int ix, int dim) {
    #ifdef FIXED_STRIDE
        int iz = ((ix / THREADSIZE) * THREADSIZE) * dim * 4*NF  + ix % THREADSIZE;
    #else
        int iz = ix;
    #endif
    COMPLEX* in_cpx = (COMPLEX*)in;
    for (int s_comp = 0; s_comp < 4; ++s_comp) {
        for (int vec_comp = 0; vec_comp < NF; ++vec_comp) {
            (*s).c[s_comp].c[vec_comp] = in_cpx[iz];
            iz+=THREADSIZE; 
        }
    }
}

template <typename SITE_TYPE, typename COMPLEX>
__device__ void write_gpu(int stride, SITE_TYPE *s, SITE_TYPE *out, int ix, int dim) {
   #ifdef FIXED_STRIDE
        int iz = ((ix / THREADSIZE) * THREADSIZE) * dim * 4*NF  + ix % THREADSIZE;
    #else
        int iz = ix;
    #endif
    COMPLEX* out_cpx = (COMPLEX*)out;
    for (int s_comp = 0; s_comp < 4; ++s_comp) {
      for (int vec_comp = 0; vec_comp < NF; ++vec_comp) {
         out_cpx[iz] = (*s).c[s_comp].c[vec_comp];
         iz+=THREADSIZE;
      }
    }
}

/* Re <s1,s2> */
template<typename SITE_TYPE, typename REAL, typename COMPLEX>
__global__ void spinor_field_prod_re_gpu(SITE_TYPE* s1, SITE_TYPE* s2, double* resField, int N)
{
   SITE_TYPE site1;
   SITE_TYPE site2;
   int ix = blockIdx.x * BLOCK_SIZE + threadIdx.x;
   if (ix < N) {
      read_gpu<SITE_TYPE, COMPLEX>(N, &site1, s1, ix, 1);
      read_gpu<SITE_TYPE, COMPLEX>(N, &site2, s2, ix, 1);
      _spinor_prod_re_f(resField[ix], site1, site2);
   }
}

/* Im <s1,s2> */
template<typename SITE_TYPE, typename REAL, typename COMPLEX>
__global__ void spinor_field_prod_im_gpu(SITE_TYPE *s1, SITE_TYPE *s2, double* resField, int N)
{
   SITE_TYPE site1;
   SITE_TYPE site2;
   int ix = blockIdx.x * BLOCK_SIZE + threadIdx.x;
   if (ix < N) {
      read_gpu<SITE_TYPE, COMPLEX>(N, &site1, s1, ix, 1);
      read_gpu<SITE_TYPE, COMPLEX>(N, &site2, s2, ix, 1);
      _spinor_prod_im_f(resField[ix], site1, site2);
   }
}

/* <s1,s2> */
template<typename SITE_TYPE, typename COMPLEX>
__global__ void spinor_field_prod_gpu(SITE_TYPE* s1, SITE_TYPE* s2, hr_complex* resField, int N)
{
   SITE_TYPE site1;
   SITE_TYPE site2;
   int ix = blockIdx.x * BLOCK_SIZE + threadIdx.x;
   if (ix < N) {
      read_gpu<SITE_TYPE, COMPLEX>(N, &site1, s1, ix, 1);
      read_gpu<SITE_TYPE, COMPLEX>(N, &site2, s2, ix, 1);
      _spinor_prod_f(resField[ix], site1, site2);
   }
}

/* Re <g5*s1,s2> */
template<typename SITE_TYPE, typename COMPLEX, typename REAL>
__global__ void spinor_field_g5_prod_re_gpu(SITE_TYPE* s1, SITE_TYPE* s2, double* resField, int N)
{
   SITE_TYPE site1;
   SITE_TYPE site2;
   int ix = blockIdx.x * BLOCK_SIZE + threadIdx.x;
   if (ix < N) {
      read_gpu<SITE_TYPE, COMPLEX>(N, &site1, s1, ix, 1);
      read_gpu<SITE_TYPE, COMPLEX>(N, &site2, s2, ix, 1);
      _spinor_g5_prod_re_f(resField[ix], site1, site2);
   }
}

/* Im <g5*s1,s2> */
template<typename SITE_TYPE, typename COMPLEX, typename REAL>
 __global__ void spinor_field_g5_prod_im_gpu(SITE_TYPE* s1, SITE_TYPE* s2, double* resField, int N)
{
   SITE_TYPE site1;
   SITE_TYPE site2;
   int ix = blockIdx.x * BLOCK_SIZE + threadIdx.x;
   if (ix < N) {
      read_gpu<SITE_TYPE, COMPLEX>(N, &site1, s1, ix, 1);
      read_gpu<SITE_TYPE, COMPLEX>(N, &site2, s2, ix, 1);
      _spinor_g5_prod_im_f(resField[ix], site1, site2);
   }
}

/* Re <s1,s1> */
template<typename SITE_TYPE, typename COMPLEX, typename REAL>
__global__ void spinor_field_sqnorm_gpu(SITE_TYPE* s1, double* resField, int N)
{
   SITE_TYPE site1;
   int ix = blockIdx.x * BLOCK_SIZE + threadIdx.x;
   if (ix < N) {
      read_gpu<SITE_TYPE, COMPLEX>(N, &site1, s1, ix, 1);
      _spinor_prod_re_f(resField[ix], site1, site1);
   }
}

/* s1+=r*s2 r real */
template<typename SITE_TYPE, typename COMPLEX, typename REAL>
__global__ void spinor_field_mul_add_assign_gpu(SITE_TYPE *s1, REAL r, SITE_TYPE *s2, int N)
{
   SITE_TYPE site1;
   SITE_TYPE site2;
   int ix = blockIdx.x * BLOCK_SIZE + threadIdx.x;
   if (ix < N) {
      read_gpu<SITE_TYPE, COMPLEX>(N, &site1, s1, ix, 1);
      read_gpu<SITE_TYPE, COMPLEX>(N, &site2, s2, ix, 1);
      _spinor_mul_add_assign_f(site1, r, site2);
      write_gpu<SITE_TYPE, COMPLEX>(N, &site1, s1, ix, 1);
   }
}

/* s1+=c*s2 c complex */
template<typename SITE_TYPE, typename COMPLEX>
__global__ void spinor_field_mulc_add_assign_gpu(SITE_TYPE *s1, COMPLEX c, SITE_TYPE *s2, int N)
{
   SITE_TYPE site1;
   SITE_TYPE site2;
   int ix = blockIdx.x * BLOCK_SIZE + threadIdx.x;
   if (ix < N) {
      read_gpu<SITE_TYPE, COMPLEX>(N, &site1, s1, ix, 1);
      read_gpu<SITE_TYPE, COMPLEX>(N, &site2, s2, ix, 1);
      _spinor_mulc_add_assign_f(site1, c, site2);
      write_gpu<SITE_TYPE, COMPLEX>(N, &site1, s1, ix, 1);
   }
}

/* s1=r*s2 r real */
template<typename SITE_TYPE, typename COMPLEX , typename REAL>
__global__ void spinor_field_mul_gpu(SITE_TYPE *s1, REAL r, SITE_TYPE *s2, int N)
{
   SITE_TYPE site1;
   SITE_TYPE site2;
   int ix = blockIdx.x * BLOCK_SIZE + threadIdx.x;
   if (ix < N) {
      read_gpu<SITE_TYPE, COMPLEX>(N, &site2, s2, ix, 1);
      _spinor_mul_f(site1, r, site2);
      write_gpu<SITE_TYPE, COMPLEX>(N, &site1, s1, ix, 1);
   }
}

/* s1=c*s2 c complex */
template<typename SITE_TYPE, typename COMPLEX>
__global__ void spinor_field_mulc_gpu(SITE_TYPE *s1, COMPLEX c, SITE_TYPE *s2, int N)
{
   SITE_TYPE site1;
   SITE_TYPE site2;
   int ix = blockIdx.x * BLOCK_SIZE + threadIdx.x;
   if (ix < N) {
      read_gpu<SITE_TYPE, COMPLEX>(N, &site2, s2, ix, 1);
      _spinor_mulc_f(site1, c, site2);
      write_gpu<SITE_TYPE, COMPLEX>(N, &site1, s1, ix, 1);
   }
}

/* r=s1+s2 */
template<typename SITE_TYPE, typename COMPLEX>
__global__ void spinor_field_add_gpu(SITE_TYPE *r, SITE_TYPE *s1, SITE_TYPE *s2, int N)
{
   SITE_TYPE site1;
   SITE_TYPE site2;
   SITE_TYPE res;
   int ix = blockIdx.x * BLOCK_SIZE + threadIdx.x;
   if (ix < N) {
      read_gpu<SITE_TYPE, COMPLEX>(N, &site1, s1, ix, 1);
      read_gpu<SITE_TYPE, COMPLEX>(N, &site2, s2, ix, 1);
      _spinor_add_f(res, site1, site2);
      write_gpu<SITE_TYPE, COMPLEX>(N, &res, r, ix, 1);
   }
}

/* r=s1-s2 */
template<typename SITE_TYPE, typename COMPLEX>
__global__ void spinor_field_sub_gpu(SITE_TYPE *r, SITE_TYPE * s1, SITE_TYPE *s2, int N)
{
   SITE_TYPE site1;
   SITE_TYPE site2;
   SITE_TYPE res;
   int ix = blockIdx.x * BLOCK_SIZE + threadIdx.x;
   if (ix < N) {
      read_gpu<SITE_TYPE, COMPLEX>(N, &site1, s1, ix, 1);
      read_gpu<SITE_TYPE, COMPLEX>(N, &site2, s2, ix, 1);
      _spinor_sub_f(res, site1, site2);
      write_gpu<SITE_TYPE, COMPLEX>(N, &res, r, ix, 1);
   }
}

/* s1+=s2 */
template<typename SITE_TYPE, typename COMPLEX>
__global__ void spinor_field_add_assign_gpu(SITE_TYPE *s1, SITE_TYPE *s2, int N)
{
   SITE_TYPE site1;
   SITE_TYPE site2;
   int ix = blockIdx.x * BLOCK_SIZE + threadIdx.x;
   if (ix < N) {
      read_gpu<SITE_TYPE, COMPLEX>(N, &site1, s1, ix, 1);
      read_gpu<SITE_TYPE, COMPLEX>(N, &site2, s2, ix, 1);
      _spinor_add_assign_f(site1, site2);
      write_gpu<SITE_TYPE, COMPLEX>(N, &site1, s1, ix, 1);
   }
}

/* s1-=s2 */
template<typename SITE_TYPE, typename COMPLEX>
__global__ void spinor_field_sub_assign_gpu(SITE_TYPE* s1, SITE_TYPE *s2, int N)
{
   SITE_TYPE site1;
   SITE_TYPE site2;
   int ix = blockIdx.x * BLOCK_SIZE + threadIdx.x;
   if (ix < N) {
      read_gpu<SITE_TYPE, COMPLEX>(N, &site1, s1, ix, 1);
      read_gpu<SITE_TYPE, COMPLEX>(N, &site2, s2, ix, 1);
      _spinor_sub_assign_f(site1, site2);
      write_gpu<SITE_TYPE, COMPLEX>(N, &site1, s1, ix, 1);
   }
}

/* s1=0 */
template<typename SITE_TYPE, typename COMPLEX>
__global__ void spinor_field_zero_gpu(SITE_TYPE *s1, int N)
{
   SITE_TYPE site;  
   int ix = blockIdx.x * BLOCK_SIZE + threadIdx.x;
   if (ix < N) {
      _spinor_zero_f(site);
      write_gpu<SITE_TYPE, COMPLEX>(N, &site, s1, ix, 1);
   }
}

/* s1=-s2 */
template<typename SITE_TYPE, typename COMPLEX>
__global__ void spinor_field_minus_gpu(SITE_TYPE * s1, SITE_TYPE * s2, int N)
{
   SITE_TYPE site1;
   SITE_TYPE site2;
   int ix = blockIdx.x * BLOCK_SIZE + threadIdx.x;
   if (ix < N) {
      read_gpu<SITE_TYPE, COMPLEX>(N, &site2, s2, ix, 1);
      _spinor_minus_f(site1, site2);
      write_gpu<SITE_TYPE, COMPLEX>(N, &site1, s1, ix, 1);
   }
}

/* s1=-s1 */
template<typename SITE_TYPE, typename COMPLEX>
__global__ void spinor_field_minus_assign_gpu(SITE_TYPE* s1, int N)
{
   SITE_TYPE site1;
   int ix = blockIdx.x * BLOCK_SIZE + threadIdx.x;
   if (ix < N) {
      read_gpu<SITE_TYPE, COMPLEX>(N, &site1, s1, ix, 1);
      _spinor_minus_f(site1, site1);
      write_gpu<SITE_TYPE, COMPLEX>(N, &site1, s1, ix, 1);
   }
}

/* s1=r1*s2+r2*s3 r1,r2 real */
template<typename SITE_TYPE, typename COMPLEX , typename REAL>
__global__ void spinor_field_lc_gpu(SITE_TYPE *s1, REAL r1, SITE_TYPE *s2, REAL r2, SITE_TYPE *s3, int N)
{
   SITE_TYPE site1;
   SITE_TYPE site2;
   SITE_TYPE site3;
   int ix = blockIdx.x * BLOCK_SIZE + threadIdx.x;
   if (ix < N) {
      read_gpu<SITE_TYPE, COMPLEX>(N, &site2, s2, ix, 1);
      read_gpu<SITE_TYPE, COMPLEX>(N, &site3, s3, ix, 1);
      _spinor_lc_f(site1, r1, site2, r2, site3);
      write_gpu<SITE_TYPE, COMPLEX>(N, &site1, s1, ix, 1);
   }
}

/* s1+=r1*s2+r2*s3 r1,r2 real */
template<typename SITE_TYPE, typename COMPLEX, typename REAL>
__global__ void spinor_field_lc_add_assign_gpu(SITE_TYPE *s1, REAL r1, SITE_TYPE *s2, REAL r2, SITE_TYPE *s3, int N)
{
   SITE_TYPE site1;
   SITE_TYPE site2;
   SITE_TYPE site3;
   int ix = blockIdx.x * BLOCK_SIZE + threadIdx.x;
   if (ix < N) {
      read_gpu<SITE_TYPE, COMPLEX>(N, &site1, s1, ix, 1);
      read_gpu<SITE_TYPE, COMPLEX>(N, &site2, s2, ix, 1);
      read_gpu<SITE_TYPE, COMPLEX>(N, &site3, s3, ix, 1);
      _spinor_lc_add_assign_f(site1, r1, site2, r2, site3);
      write_gpu<SITE_TYPE, COMPLEX>(N, &site1, s1, ix, 1);
   }
}

/* s1=cd1*s2+cd2*s3 cd1, cd2 complex */
template<typename SITE_TYPE, typename COMPLEX>
__global__ void spinor_field_clc_gpu(SITE_TYPE *s1, COMPLEX c1, SITE_TYPE *s2, COMPLEX c2, SITE_TYPE *s3, int N)
{
   SITE_TYPE site1;
   SITE_TYPE site2;
   SITE_TYPE site3;
   int ix = blockIdx.x * BLOCK_SIZE + threadIdx.x;
   if (ix < N) {
      read_gpu<SITE_TYPE, COMPLEX>(N, &site2, s2, ix, 1);
      read_gpu<SITE_TYPE, COMPLEX>(N, &site3, s3, ix, 1);
      _spinor_clc_f(site1, c1, site2, c2, site3);
      write_gpu<SITE_TYPE, COMPLEX>(N, &site1, s1, ix, 1);
   }
}

/* s1+=cd1*s2+cd2*s3 cd1,cd2 complex */
template<typename SITE_TYPE, typename COMPLEX>
__global__ void spinor_field_clc_add_assign_gpu(SITE_TYPE *s1, COMPLEX c1, SITE_TYPE *s2, COMPLEX c2, SITE_TYPE *s3, int N)
{
   SITE_TYPE site1;
   SITE_TYPE site2;
   SITE_TYPE site3;
   int ix = blockIdx.x * BLOCK_SIZE + threadIdx.x;
   if (ix < N) {
      read_gpu<SITE_TYPE, COMPLEX>(N, &site1, s1, ix, 1);
      read_gpu<SITE_TYPE, COMPLEX>(N, &site2, s2, ix, 1);
      read_gpu<SITE_TYPE, COMPLEX>(N, &site3, s3, ix, 1);
      _spinor_clc_add_assign_f(site1, c1, site2, c2, site3);
      write_gpu<SITE_TYPE, COMPLEX>(N, &site1, s1, ix, 1);
   }
}

/* s1=g5*s2  */
template<typename SITE_TYPE, typename COMPLEX>
__global__ void spinor_field_g5_gpu(SITE_TYPE *s1, SITE_TYPE *s2, int N)
{
   SITE_TYPE site1;
   SITE_TYPE site2;
   int ix = blockIdx.x * BLOCK_SIZE + threadIdx.x;
   if (ix < N) {
      read_gpu<SITE_TYPE, COMPLEX>(N, &site2, s2, ix, 1);
      _spinor_g5_f(site1, site2);
      write_gpu<SITE_TYPE, COMPLEX>(N, &site1, s1, ix, 1);
   }
}

/* s1=g5*s1 */
template<typename SITE_TYPE, typename COMPLEX>
__global__ void spinor_field_g5_assign_gpu(SITE_TYPE* s1, int N)
{
   SITE_TYPE site1;
   int ix = blockIdx.x * BLOCK_SIZE + threadIdx.x;
   if (ix < N) {
      read_gpu<SITE_TYPE, COMPLEX>(N, &site1, s1, ix, 1);
      _spinor_g5_assign_f(site1);
      write_gpu<SITE_TYPE, COMPLEX>(N, &site1, s1, ix, 1);
   }
}

/* s1+=c*g5*s2 c complex  */
template<typename SITE_TYPE, typename COMPLEX>
__global__ void spinor_field_g5_mulc_add_assign_gpu(SITE_TYPE *s1, COMPLEX c, SITE_TYPE *s2, int N)
{
   SITE_TYPE site1;
   SITE_TYPE site2;
   SITE_TYPE tmp;
   int ix = blockIdx.x * BLOCK_SIZE + threadIdx.x;
   if (ix < N) {
      read_gpu<SITE_TYPE, COMPLEX>(N, &site1, s1, ix, 1);
      read_gpu<SITE_TYPE, COMPLEX>(N, &site2, s2, ix, 1);
      _spinor_mulc_f(tmp, c, site2);
      _spinor_g5_f(site2, tmp);
      _spinor_add_assign_f(site1, site2);
      write_gpu<SITE_TYPE, COMPLEX>(N, &site1, s1, ix, 1);
   }
}

/* tools per eva.c  */
template<typename SITE_TYPE, typename COMPLEX , typename REAL>
__global__ void spinor_field_lc1_gpu(REAL r, SITE_TYPE *s1, SITE_TYPE *s2, int N)
{
   SITE_TYPE site1;
   SITE_TYPE site2;
   int ix = blockIdx.x * BLOCK_SIZE + threadIdx.x;
   if (ix < N) {
      read_gpu<SITE_TYPE, COMPLEX>(N, &site1, s1, ix, 1);
      read_gpu<SITE_TYPE, COMPLEX>(N, &site2, s2, ix, 1);
      _spinor_mul_f(site1, r, site2);
      write_gpu<SITE_TYPE, COMPLEX>(N, &site1, s1, ix, 1);
   }
}

template<typename SITE_TYPE, typename COMPLEX, typename REAL>
__global__ void spinor_field_lc2_gpu(REAL r1, REAL r2, SITE_TYPE *s1, SITE_TYPE *s2, int N)
{
   SITE_TYPE site1;
   SITE_TYPE site2;
   int ix = blockIdx.x * BLOCK_SIZE + threadIdx.x;
   if (ix < N) {
      read_gpu<SITE_TYPE, COMPLEX>(N, &site1, s1, ix, 1);
      read_gpu<SITE_TYPE, COMPLEX>(N, &site2, s2, ix, 1);
      _spinor_lc_f(site1, r1, site1, r2, site2);
   }
}

template<typename SITE_TYPE, typename COMPLEX , typename REAL>
__global__ void spinor_field_lc3_gpu(REAL r1, REAL r2, SITE_TYPE *s1, SITE_TYPE *s2, SITE_TYPE *s3, int N)
{
   SITE_TYPE site1;
   SITE_TYPE site2;
   SITE_TYPE site3;
   int ix = blockIdx.x * BLOCK_SIZE + threadIdx.x;
   if (ix < N) {
      read_gpu<SITE_TYPE, COMPLEX>(N, &site1, s1, ix, 1);
      read_gpu<SITE_TYPE, COMPLEX>(N, &site2, s2, ix, 1);
      _spinor_lc_f(site3, r1, site1, r2, site2);
      _spinor_minus_f(site3, site3);
      write_gpu<SITE_TYPE, COMPLEX>(N, &site3, s3, ix, 1);
   }
}

#endif
