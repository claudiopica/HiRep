#ifndef LINEAR_ALGEBRA_GPU_CU
#define LINEAR_ALGEBRA_GPU_CU

#include "geometry.h"
#include "Utils/generics.h"
#include "inverters.h"

/* Re <s1,s2> */
template <unsigned int FIELD_DIM, typename REAL, typename SITE_TYPE>
__global__ void prod_re_gpu(SITE_TYPE *s1, SITE_TYPE *s2, quad_double *resField, int N) {
    SITE_TYPE site1;
    SITE_TYPE site2;
    for (int ix = blockIdx.x * blockDim.x + threadIdx.x; ix < N; ix += blockDim.x * gridDim.x) {
        resField[ix] = quad_double(0.0);
        for (int mu = 0; mu < FIELD_DIM; ++mu) {
            read_gpu<REAL>(N, &site1, s1, ix, mu, FIELD_DIM);
            read_gpu<REAL>(N, &site2, s2, ix, mu, FIELD_DIM);
            resField[ix].val += prod_re(&site1, &site2);
        }
    }
}

/* Im <s1,s2> */
template <unsigned int FIELD_DIM, typename REAL, typename SITE_TYPE>
__global__ void prod_im_gpu(SITE_TYPE *s1, SITE_TYPE *s2, quad_double *resField, int N) {
    SITE_TYPE site1;
    SITE_TYPE site2;
    for (int ix = blockIdx.x * blockDim.x + threadIdx.x; ix < N; ix += blockDim.x * gridDim.x) {
        resField[ix] = quad_double(0.0);
        for (int mu = 0; mu < FIELD_DIM; ++mu) {
            read_gpu<REAL>(N, &site1, s1, ix, mu, FIELD_DIM);
            read_gpu<REAL>(N, &site2, s2, ix, mu, FIELD_DIM);
            resField[ix].val += prod_im(&site1, &site2);
        }
    }
}

/* <s1,s2> */
template <unsigned int FIELD_DIM, typename REAL, typename SITE_TYPE>
__global__ void prod_gpu(SITE_TYPE *s1, SITE_TYPE *s2, hr_complex *resField, int N) {
    SITE_TYPE site1;
    SITE_TYPE site2;
    for (int ix = blockIdx.x * blockDim.x + threadIdx.x; ix < N; ix += blockDim.x * gridDim.x) {
        _complex_0(resField[ix]);
        for (int mu = 0; mu < FIELD_DIM; mu++) {
            read_gpu<REAL>(N, &site1, s1, ix, mu, FIELD_DIM);
            read_gpu<REAL>(N, &site2, s2, ix, mu, FIELD_DIM);
            _complex_add(resField[ix], resField[ix], prod(&site1, &site2));
        }
    }
}

/* Re <g5*s1,s2> */
template <unsigned int FIELD_DIM, typename REAL, typename SITE_TYPE>
__global__ void g5_prod_re_gpu(SITE_TYPE *s1, SITE_TYPE *s2, double *resField, int N) {
    SITE_TYPE site1;
    SITE_TYPE site2;
    for (int ix = blockIdx.x * blockDim.x + threadIdx.x; ix < N; ix += blockDim.x * gridDim.x) {
        resField[ix] = 0.0;
        for (int mu = 0; mu < FIELD_DIM; mu++) {
            read_gpu<REAL>(N, &site1, s1, ix, mu, FIELD_DIM);
            read_gpu<REAL>(N, &site2, s2, ix, mu, FIELD_DIM);
            resField[ix] += g5_prod_re(&site1, &site2);
        }
    }
}

/* Im <g5*s1,s2> */
template <unsigned int FIELD_DIM, typename REAL, typename SITE_TYPE>
__global__ void g5_prod_im_gpu(SITE_TYPE *s1, SITE_TYPE *s2, double *resField, int N) {
    SITE_TYPE site1;
    SITE_TYPE site2;
    for (int ix = blockIdx.x * blockDim.x + threadIdx.x; ix < N; ix += blockDim.x * gridDim.x) {
        resField[ix] = 0.0;
        for (int mu = 0; mu < FIELD_DIM; mu++) {
            read_gpu<REAL>(N, &site1, s1, ix, mu, FIELD_DIM);
            read_gpu<REAL>(N, &site2, s2, ix, mu, FIELD_DIM);
            resField[ix] += g5_prod_im(&site1, &site2);
        }
    }
}

/* Re <s1,s1> */
template <unsigned int FIELD_DIM, typename REAL, typename SITE_TYPE>
__global__ void sqnorm_gpu(SITE_TYPE *s1, quad_double *resField, int N) {
    SITE_TYPE site1;
    for (int ix = blockIdx.x * blockDim.x + threadIdx.x; ix < N; ix += blockDim.x * gridDim.x) {
        resField[ix] = quad_double(0.0);
        for (int mu = 0; mu < FIELD_DIM; mu++) {
            read_gpu<REAL>(N, &site1, s1, ix, mu, FIELD_DIM);
            resField[ix].val += prod_re(&site1, &site1);
        }
    }
}

/* Re <s1,s1> */
template <unsigned int FIELD_DIM, typename REAL, typename SITE_TYPE>
__global__ void max_gpu(SITE_TYPE *s1, double *resField, int N) {
    SITE_TYPE site1;
    for (int ix = blockIdx.x * blockDim.x + threadIdx.x; ix < N; ix += blockDim.x * gridDim.x) {
        resField[ix] = 0.0;
        for (int mu = 0; mu < FIELD_DIM; mu++) {
            read_gpu<REAL>(N, &site1, s1, ix, mu, FIELD_DIM);
            resField[ix] = resField[ix] > max(&site1) ? resField[ix] : max(&site1);
        }
    }
}

/* s1=s2 */
template <unsigned int FIELD_DIM, typename REAL, typename SITE_TYPE>
__global__ void id_gpu(SITE_TYPE *s1, SITE_TYPE *s2, int N) {
    SITE_TYPE site;
    for (int ix = blockIdx.x * blockDim.x + threadIdx.x; ix < N; ix += blockDim.x * gridDim.x) {
        for (int mu = 0; mu < FIELD_DIM; ++mu) {
            read_gpu<REAL>(N, &site, s2, ix, mu, FIELD_DIM);
            write_gpu<REAL>(N, &site, s1, ix, mu, FIELD_DIM);
        }
    }
}

/* s1+=r*s2 r real */
template <unsigned int FIELD_DIM, typename REAL, typename SITE_TYPE>
__global__ void mul_add_assign_gpu(SITE_TYPE *s1, REAL r, SITE_TYPE *s2, int N) {
    SITE_TYPE site1;
    SITE_TYPE site2;
    for (int ix = blockIdx.x * blockDim.x + threadIdx.x; ix < N; ix += blockDim.x * gridDim.x) {
        for (int mu = 0; mu < FIELD_DIM; ++mu) {
            read_gpu<REAL>(N, &site1, s1, ix, mu, FIELD_DIM);
            read_gpu<REAL>(N, &site2, s2, ix, mu, FIELD_DIM);
            mul_add_assign(&site1, r, &site2);
            write_gpu<REAL>(N, &site1, s1, ix, mu, FIELD_DIM);
        }
    }
}

/* s1+=c*s2 c complex */
template <unsigned int FIELD_DIM, typename REAL, typename SITE_TYPE, typename COMPLEX>
__global__ void mulc_add_assign_gpu(SITE_TYPE *s1, COMPLEX c, SITE_TYPE *s2, int N) {
    SITE_TYPE site1;
    SITE_TYPE site2;
    for (int ix = blockIdx.x * blockDim.x + threadIdx.x; ix < N; ix += blockDim.x * gridDim.x) {
        for (int mu = 0; mu < FIELD_DIM; ++mu) {
            read_gpu<REAL>(N, &site1, s1, ix, mu, FIELD_DIM);
            read_gpu<REAL>(N, &site2, s2, ix, mu, FIELD_DIM);
            mulc_add_assign(&site1, c, &site2);
            write_gpu<REAL>(N, &site1, s1, ix, mu, FIELD_DIM);
        }
    }
}

/* s1+=c*g5*s2 c complex  */
template <unsigned int FIELD_DIM, typename REAL, typename SITE_TYPE, typename COMPLEX>
__global__ void g5_mulc_add_assign_gpu(SITE_TYPE *s1, COMPLEX c, SITE_TYPE *s2, int N) {
    SITE_TYPE site1;
    SITE_TYPE site2;
    SITE_TYPE tmp;
    for (int ix = blockIdx.x * blockDim.x + threadIdx.x; ix < N; ix += blockDim.x * gridDim.x) {
        for (int mu = 0; mu < FIELD_DIM; mu++) {
            read_gpu<REAL>(N, &site1, s1, ix, mu, FIELD_DIM);
            read_gpu<REAL>(N, &site2, s2, ix, mu, FIELD_DIM);
            mulc(&tmp, c, &site2);
            g5(&site2, &tmp);
            add_assign(&site1, &site2);
            write_gpu<REAL>(N, &site1, s1, ix, mu, FIELD_DIM);
        }
    }
}

/* s1=r*s2 r real */
template <unsigned int FIELD_DIM, typename SITE_TYPE, typename REAL>
__global__ void mul_gpu(SITE_TYPE *s1, REAL r, SITE_TYPE *s2, int N) {
    SITE_TYPE site1;
    SITE_TYPE site2;
    for (int ix = blockIdx.x * blockDim.x + threadIdx.x; ix < N; ix += blockDim.x * gridDim.x) {
        for (int mu = 0; mu < FIELD_DIM; mu++) {
            read_gpu<REAL>(N, &site2, s2, ix, mu, FIELD_DIM);
            mul(&site1, r, &site2);
            write_gpu<REAL>(N, &site1, s1, ix, mu, FIELD_DIM);
        }
    }
}

/* s1=c*s2 c complex */
template <unsigned int FIELD_DIM, typename REAL, typename SITE_TYPE, typename COMPLEX>
__global__ void mulc_gpu(SITE_TYPE *s1, COMPLEX c, SITE_TYPE *s2, int N) {
    SITE_TYPE site1;
    SITE_TYPE site2;
    for (int ix = blockIdx.x * blockDim.x + threadIdx.x; ix < N; ix += blockDim.x * gridDim.x) {
        for (int mu = 0; mu < FIELD_DIM; mu++) {
            read_gpu<REAL>(N, &site2, s2, ix, mu, FIELD_DIM);
            mulc(&site1, c, &site2);
            write_gpu<REAL>(N, &site1, s1, ix, mu, FIELD_DIM);
        }
    }
}

/* r=s1+s2 */
template <unsigned int FIELD_DIM, typename REAL, typename SITE_TYPE>
__global__ void add_gpu(SITE_TYPE *r, SITE_TYPE *s1, SITE_TYPE *s2, int N) {
    SITE_TYPE site1;
    SITE_TYPE site2;
    SITE_TYPE res;
    for (int ix = blockIdx.x * blockDim.x + threadIdx.x; ix < N; ix += blockDim.x * gridDim.x) {
        for (int mu = 0; mu < FIELD_DIM; mu++) {
            read_gpu<REAL>(N, &site1, s1, ix, mu, FIELD_DIM);
            read_gpu<REAL>(N, &site2, s2, ix, mu, FIELD_DIM);
            add(&res, &site1, &site2);
            write_gpu<REAL>(N, &res, r, ix, mu, FIELD_DIM);
        }
    }
}

/* r=s1-s2 */
template <unsigned int FIELD_DIM, typename REAL, typename SITE_TYPE>
__global__ void sub_gpu(SITE_TYPE *r, SITE_TYPE *s1, SITE_TYPE *s2, int N) {
    SITE_TYPE site1;
    SITE_TYPE site2;
    SITE_TYPE res;
    for (int ix = blockIdx.x * blockDim.x + threadIdx.x; ix < N; ix += blockDim.x * gridDim.x) {
        for (int mu = 0; mu < FIELD_DIM; mu++) {
            read_gpu<REAL>(N, &site1, s1, ix, mu, FIELD_DIM);
            read_gpu<REAL>(N, &site2, s2, ix, mu, FIELD_DIM);
            sub(&res, &site1, &site2);
            write_gpu<REAL>(N, &res, r, ix, mu, FIELD_DIM);
        }
    }
}

/* s1+=s2 */
template <unsigned int FIELD_DIM, typename REAL, typename SITE_TYPE>
__global__ void add_assign_gpu(SITE_TYPE *s1, SITE_TYPE *s2, int N) {
    SITE_TYPE site1;
    SITE_TYPE site2;
    for (int ix = blockIdx.x * blockDim.x + threadIdx.x; ix < N; ix += blockDim.x * gridDim.x) {
        for (int mu = 0; mu < FIELD_DIM; mu++) {
            read_gpu<REAL>(N, &site1, s1, ix, mu, FIELD_DIM);
            read_gpu<REAL>(N, &site2, s2, ix, mu, FIELD_DIM);
            add_assign(&site1, &site2);
            write_gpu<REAL>(N, &site1, s1, ix, mu, FIELD_DIM);
        }
    }
}

/* s1-=s2 */
template <unsigned int FIELD_DIM, typename REAL, typename SITE_TYPE>
__global__ void sub_assign_gpu(SITE_TYPE *s1, SITE_TYPE *s2, int N) {
    SITE_TYPE site1;
    SITE_TYPE site2;
    for (int ix = blockIdx.x * blockDim.x + threadIdx.x; ix < N; ix += blockDim.x * gridDim.x) {
        for (int mu = 0; mu < FIELD_DIM; mu++) {
            read_gpu<REAL>(N, &site1, s1, ix, mu, FIELD_DIM);
            read_gpu<REAL>(N, &site2, s2, ix, mu, FIELD_DIM);
            sub_assign(&site1, &site2);
            write_gpu<REAL>(N, &site1, s1, ix, mu, FIELD_DIM);
        }
    }
}

/* s1=0 */
template <unsigned int FIELD_DIM, typename REAL, typename SITE_TYPE> __global__ void zero_gpu(SITE_TYPE *s1, int N) {
    SITE_TYPE site;
    for (int ix = blockIdx.x * blockDim.x + threadIdx.x; ix < N; ix += blockDim.x * gridDim.x) {
        for (int mu = 0; mu < FIELD_DIM; mu++) {
            zero(&site);
            write_gpu<REAL>(N, &site, s1, ix, mu, FIELD_DIM);
        }
    }
}

/* s1=-s2 */
template <unsigned int FIELD_DIM, typename REAL, typename SITE_TYPE>
__global__ void minus_gpu(SITE_TYPE *s1, SITE_TYPE *s2, int N) {
    SITE_TYPE site1;
    SITE_TYPE site2;
    for (int ix = blockIdx.x * blockDim.x + threadIdx.x; ix < N; ix += blockDim.x * gridDim.x) {
        for (int mu = 0; mu < FIELD_DIM; mu++) {
            read_gpu<REAL>(N, &site2, s2, ix, mu, FIELD_DIM);
            minus(&site1, &site2);
            write_gpu<REAL>(N, &site1, s1, ix, mu, FIELD_DIM);
        }
    }
}

/* s1=r1*s2+r2*s3 r1,r2 real */
template <unsigned int FIELD_DIM, typename REAL, typename SITE_TYPE>
__global__ void lc_gpu(SITE_TYPE *s1, REAL r1, SITE_TYPE *s2, REAL r2, SITE_TYPE *s3, int N) {
    SITE_TYPE site1;
    SITE_TYPE site2;
    SITE_TYPE site3;
    for (int ix = blockIdx.x * blockDim.x + threadIdx.x; ix < N; ix += blockDim.x * gridDim.x) {
        for (int mu = 0; mu < FIELD_DIM; mu++) {
            read_gpu<REAL>(N, &site2, s2, ix, mu, FIELD_DIM);
            read_gpu<REAL>(N, &site3, s3, ix, mu, FIELD_DIM);
            lc(&site1, r1, &site2, r2, &site3);
            write_gpu<REAL>(N, &site1, s1, ix, mu, FIELD_DIM);
        }
    }
}

/* s1+=r1*s2+r2*s3 r1,r2 real */
template <unsigned int FIELD_DIM, typename SITE_TYPE, typename REAL>
__global__ void lc_add_assign_gpu(SITE_TYPE *s1, REAL r1, SITE_TYPE *s2, REAL r2, SITE_TYPE *s3, int N) {
    SITE_TYPE site1;
    SITE_TYPE site2;
    SITE_TYPE site3;
    for (int ix = blockIdx.x * blockDim.x + threadIdx.x; ix < N; ix += blockDim.x * gridDim.x) {
        for (int mu = 0; mu < FIELD_DIM; mu++) {
            read_gpu<REAL>(N, &site1, s1, ix, mu, FIELD_DIM);
            read_gpu<REAL>(N, &site2, s2, ix, mu, FIELD_DIM);
            read_gpu<REAL>(N, &site3, s3, ix, mu, FIELD_DIM);
            lc_add_assign(&site1, r1, &site2, r2, &site3);
            write_gpu<REAL>(N, &site1, s1, ix, mu, FIELD_DIM);
        }
    }
}

/* s1=cd1*s2+cd2*s3 cd1, cd2 complex */
template <unsigned int FIELD_DIM, typename REAL, typename SITE_TYPE, typename COMPLEX>
__global__ void clc_gpu(SITE_TYPE *s1, COMPLEX c1, SITE_TYPE *s2, COMPLEX c2, SITE_TYPE *s3, int N) {
    SITE_TYPE site1;
    SITE_TYPE site2;
    SITE_TYPE site3;
    for (int ix = blockIdx.x * blockDim.x + threadIdx.x; ix < N; ix += blockDim.x * gridDim.x) {
        for (int mu = 0; mu < FIELD_DIM; mu++) {
            read_gpu<REAL>(N, &site2, s2, ix, mu, FIELD_DIM);
            read_gpu<REAL>(N, &site3, s3, ix, mu, FIELD_DIM);
            clc(&site1, c1, &site2, c2, &site3);
            write_gpu<REAL>(N, &site1, s1, ix, mu, FIELD_DIM);
        }
    }
}

/* s1+=cd1*s2+cd2*s3 cd1,cd2 complex */
template <unsigned int FIELD_DIM, typename REAL, typename SITE_TYPE, typename COMPLEX>
__global__ void clc_add_assign_gpu(SITE_TYPE *s1, COMPLEX c1, SITE_TYPE *s2, COMPLEX c2, SITE_TYPE *s3, int N) {
    SITE_TYPE site1;
    SITE_TYPE site2;
    SITE_TYPE site3;
    for (int ix = blockIdx.x * blockDim.x + threadIdx.x; ix < N; ix += blockDim.x * gridDim.x) {
        for (int mu = 0; mu < FIELD_DIM; mu++) {
            read_gpu<REAL>(N, &site1, s1, ix, mu, FIELD_DIM);
            read_gpu<REAL>(N, &site2, s2, ix, mu, FIELD_DIM);
            read_gpu<REAL>(N, &site3, s3, ix, mu, FIELD_DIM);
            clc_add_assign(&site1, c1, &site2, c2, &site3);
            write_gpu<REAL>(N, &site1, s1, ix, mu, FIELD_DIM);
        }
    }
}

/* s1=g5*s2  */
template <unsigned int FIELD_DIM, typename REAL, typename SITE_TYPE>
__global__ void g5_gpu(SITE_TYPE *s1, SITE_TYPE *s2, int N) {
    SITE_TYPE site1;
    SITE_TYPE site2;
    for (int ix = blockIdx.x * blockDim.x + threadIdx.x; ix < N; ix += blockDim.x * gridDim.x) {
        for (int mu = 0; mu < FIELD_DIM; mu++) {
            read_gpu<REAL>(N, &site2, s2, ix, mu, FIELD_DIM);
            g5(&site1, &site2);
            write_gpu<REAL>(N, &site1, s1, ix, mu, FIELD_DIM);
        }
    }
}

/* s1=g5*s1 */
template <unsigned int FIELD_DIM, typename REAL, typename SITE_TYPE> __global__ void g5_assign_gpu(SITE_TYPE *s1, int N) {
    SITE_TYPE site1;
    for (int ix = blockIdx.x * blockDim.x + threadIdx.x; ix < N; ix += blockDim.x * gridDim.x) {
        for (int mu = 0; mu < FIELD_DIM; mu++) {
            read_gpu<REAL>(N, &site1, s1, ix, mu, FIELD_DIM);
            g5_assign(&site1);
            write_gpu<REAL>(N, &site1, s1, ix, mu, FIELD_DIM);
        }
    }
}

template <unsigned int FIELD_DIM, typename REAL, typename SITE_TYPE>
__global__ void g0_gpu(SITE_TYPE *s1, SITE_TYPE *s2, int N) {
    SITE_TYPE site1;
    SITE_TYPE site2;
    for (int ix = blockIdx.x * blockDim.x + threadIdx.x; ix < N; ix += blockDim.x * gridDim.x) {
        for (int mu = 0; mu < FIELD_DIM; mu++) {
            read_gpu<REAL>(N, &site2, s2, ix, mu, FIELD_DIM);
            g0(&site1, &site2);
            write_gpu<REAL>(N, &site1, s1, ix, mu, FIELD_DIM);
        }
    }
}

template <unsigned int FIELD_DIM, typename REAL, typename SITE_TYPE>
__global__ void g1_gpu(SITE_TYPE *s1, SITE_TYPE *s2, int N) {
    SITE_TYPE site1;
    SITE_TYPE site2;
    for (int ix = blockIdx.x * blockDim.x + threadIdx.x; ix < N; ix += blockDim.x * gridDim.x) {
        for (int mu = 0; mu < FIELD_DIM; mu++) {
            read_gpu<REAL>(N, &site2, s2, ix, mu, FIELD_DIM);
            g1(&site1, &site2);
            write_gpu<REAL>(N, &site1, s1, ix, mu, FIELD_DIM);
        }
    }
}

template <unsigned int FIELD_DIM, typename REAL, typename SITE_TYPE>
__global__ void g2_gpu(SITE_TYPE *s1, SITE_TYPE *s2, int N) {
    SITE_TYPE site1;
    SITE_TYPE site2;
    for (int ix = blockIdx.x * blockDim.x + threadIdx.x; ix < N; ix += blockDim.x * gridDim.x) {
        for (int mu = 0; mu < FIELD_DIM; mu++) {
            read_gpu<REAL>(N, &site2, s2, ix, mu, FIELD_DIM);
            g2(&site1, &site2);
            write_gpu<REAL>(N, &site1, s1, ix, mu, FIELD_DIM);
        }
    }
}

template <unsigned int FIELD_DIM, typename REAL, typename SITE_TYPE>
__global__ void g3_gpu(SITE_TYPE *s1, SITE_TYPE *s2, int N) {
    SITE_TYPE site1;
    SITE_TYPE site2;
    for (int ix = blockIdx.x * blockDim.x + threadIdx.x; ix < N; ix += blockDim.x * gridDim.x) {
        for (int mu = 0; mu < FIELD_DIM; mu++) {
            read_gpu<REAL>(N, &site2, s2, ix, mu, FIELD_DIM);
            g3(&site1, &site2);
            write_gpu<REAL>(N, &site1, s1, ix, mu, FIELD_DIM);
        }
    }
}

/* tools per eva.c  */
template <unsigned int FIELD_DIM, typename SITE_TYPE, typename REAL>
__global__ void lc1_gpu(REAL r, SITE_TYPE *s1, SITE_TYPE *s2, int N) {
    SITE_TYPE site1;
    SITE_TYPE site2;
    for (int ix = blockIdx.x * blockDim.x + threadIdx.x; ix < N; ix += blockDim.x * gridDim.x) {
        for (int mu = 0; mu < FIELD_DIM; mu++) {
            read_gpu<REAL>(N, &site1, s1, ix, mu, FIELD_DIM);
            read_gpu<REAL>(N, &site2, s2, ix, mu, FIELD_DIM);
            mul_add_assign(&site1, r, &site2);
            write_gpu<REAL>(N, &site1, s1, ix, mu, FIELD_DIM);
        }
    }
}

template <unsigned int FIELD_DIM, typename SITE_TYPE, typename REAL>
__global__ void lc2_gpu(REAL r1, REAL r2, SITE_TYPE *s1, SITE_TYPE *s2, int N) {
    SITE_TYPE site1;
    SITE_TYPE site2;
    for (int ix = blockIdx.x * blockDim.x + threadIdx.x; ix < N; ix += blockDim.x * gridDim.x) {
        for (int mu = 0; mu < FIELD_DIM; mu++) {
            read_gpu<REAL>(N, &site1, s1, ix, mu, FIELD_DIM);
            read_gpu<REAL>(N, &site2, s2, ix, mu, FIELD_DIM);
            lc(&site1, r1, &site1, r2, &site2);
            write_gpu<REAL>(N, &site1, s1, ix, mu, FIELD_DIM);
        }
    }
}

template <unsigned int FIELD_DIM, typename SITE_TYPE, typename REAL>
__global__ void lc3_gpu(REAL r1, REAL r2, SITE_TYPE *s1, SITE_TYPE *s2, SITE_TYPE *s3, int N) {
    SITE_TYPE site1;
    SITE_TYPE site2;
    SITE_TYPE site3;
    for (int ix = blockIdx.x * blockDim.x + threadIdx.x; ix < N; ix += blockDim.x * gridDim.x) {
        for (int mu = 0; mu < FIELD_DIM; mu++) {
            read_gpu<REAL>(N, &site1, s1, ix, mu, FIELD_DIM);
            read_gpu<REAL>(N, &site2, s2, ix, mu, FIELD_DIM);
            read_gpu<REAL>(N, &site3, s3, ix, mu, FIELD_DIM);
            lc_add_assign(&site3, r1, &site1, r2, &site2);
            minus(&site3, &site3);
            write_gpu<REAL>(N, &site3, s3, ix, mu, FIELD_DIM);
        }
    }
}

#endif
