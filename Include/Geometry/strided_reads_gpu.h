// These reading functions do not work with real representations yet (SAM)

#ifndef STRIDED_READS_GPU_HPP
#define STRIDED_READS_GPU_HPP

#ifdef FIXED_STRIDE
    #define THREADSIZE 32
#else
    #define THREADSIZE 1
#endif

//#include "libhr_core.h"
#include "geometry.h"

#ifdef __cplusplus

template<typename COMPLEX, typename FIELD_TYPE, typename SITE_TYPE>
__host__ __device__ void read_gpu(int stride, SITE_TYPE *s, const FIELD_TYPE *in, int ix, int comp, int dim) {
    const int field_dim = sizeof(FIELD_TYPE)/sizeof(COMPLEX);
    const int n_components = sizeof(SITE_TYPE)/sizeof(COMPLEX);
    #ifdef FIXED_STRIDE
        int iz = ((ix / THREADSIZE) * THREADSIZE) * dim * field_dim  + (ix % THREADSIZE) + ((comp)*n_components)*(THREADSIZE);
    #else
        int iz = ix + ((comp)*n_components)*(THREADSIZE);
    #endif
    COMPLEX* in_cpx = (COMPLEX*)in;
    COMPLEX* in_comp_cpx = (COMPLEX*)s;
    for (int i = 0; i < n_components; ++i) {
         in_comp_cpx[i] = in_cpx[iz];
         iz+=THREADSIZE;
    }
}

template <typename COMPLEX, typename FIELD_TYPE, typename SITE_TYPE>
__host__ __device__ void write_gpu(int stride, SITE_TYPE *s, FIELD_TYPE *out, int ix, int comp, int dim) {
    const int field_dim = sizeof(FIELD_TYPE)/sizeof(COMPLEX);
    const int n_components = sizeof(SITE_TYPE)/sizeof(COMPLEX);
    #ifdef FIXED_STRIDE
        int iz = ((ix / THREADSIZE) * THREADSIZE) * dim * field_dim  + (ix % THREADSIZE) + ((comp)*n_components)*(THREADSIZE);
    #else
        int iz = ix + ((comp)*n_components)*(THREADSIZE);
    #endif
    COMPLEX* out_cpx = (COMPLEX*)out;
    COMPLEX* out_comp_cpx = (COMPLEX*)s;
    for (int i = 0; i < n_components; ++i) {
        out_cpx[iz] = out_comp_cpx[i];
        iz += THREADSIZE;
    }
}

template<typename COMPLEX, typename VECTOR_TYPE, typename SITE_TYPE>
__device__ void in_spinor_field(VECTOR_TYPE *v, SITE_TYPE *in, int iy, int comp) {
    read_gpu<COMPLEX>(0, v, in, iy, comp, 1);
}

template<typename COMPLEX, typename GAUGE_TYPE>
__device__ void in_gauge_field(GAUGE_TYPE *u, const GAUGE_TYPE *in, int ix, int iy, int comp, int dir) {
    if (dir == UP) {
        read_gpu<COMPLEX>(0, u, in, ix, comp, 4);
    } else if (dir == DOWN) {
        read_gpu<COMPLEX>(0, u, in, iy, comp, 4);
    }
}

template<typename COMPLEX, typename SITE_TYPE>
__device__ void write_out_spinor_field(SITE_TYPE *r, SITE_TYPE *in, int ix) {
    write_gpu<COMPLEX>(0, r, in, ix, 0, 1);
}

#endif

#ifdef __cplusplus
    extern "C" {
#endif

#define _FIELD_NAME spinor_field_f
#define _FIELD_TYPE spinor_field
#define _SITE_TYPE suNf_spinor
#define _FIELD_DIM 1
#define _REAL double
#include "TMPL/strided_reads_gpu.hpp.tmpl"

#define _FIELD_NAME spinor_field_f_flt
#define _FIELD_TYPE spinor_field_flt
#define _SITE_TYPE suNf_spinor_flt
#define _FIELD_DIM 1
#define _REAL float
#include "TMPL/strided_reads_gpu.hpp.tmpl"

#define _FIELD_NAME sfield
#define _FIELD_TYPE scalar_field
#define _SITE_TYPE double
#define _FIELD_DIM 1
#define _REAL double
#include "TMPL/strided_reads_gpu.hpp.tmpl"

#define _FIELD_NAME gfield
#define _FIELD_TYPE suNg_field
#define _SITE_TYPE suNg
#define _FIELD_DIM 4
#define _REAL double
#include "TMPL/strided_reads_gpu.hpp.tmpl"

#define _FIELD_NAME gfield_flt
#define _FIELD_TYPE suNg_field_flt
#define _SITE_TYPE suNg_flt
#define _FIELD_DIM 4
#define _REAL float
#include "TMPL/strided_reads_gpu.hpp.tmpl"

#define _FIELD_NAME gfield_f
#define _FIELD_TYPE suNf_field
#define _SITE_TYPE suNf
#define _FIELD_DIM 4
#define _REAL double
#include "TMPL/strided_reads_gpu.hpp.tmpl"

#define _FIELD_NAME gfield_f_flt
#define _FIELD_TYPE suNf_field_flt
#define _SITE_TYPE suNf_flt
#define _FIELD_DIM 4
#define _REAL float
#include "TMPL/strided_reads_gpu.hpp.tmpl"

#define _FIELD_NAME suNg_scalar_field
#define _FIELD_TYPE suNg_scalar_field
#define _SITE_TYPE suNg_vector
#define _FIELD_DIM 1
#define _REAL double
#include "TMPL/strided_reads_gpu.hpp.tmpl"

#define _FIELD_NAME avfield
#define _FIELD_TYPE suNg_av_field
#define _SITE_TYPE suNg_algebra_vector
#define _FIELD_DIM 4
#define _REAL double
#include "TMPL/strided_reads_gpu.hpp.tmpl"

#define _FIELD_NAME gtransf
#define _FIELD_TYPE suNg_field
#define _SITE_TYPE suNg
#define _FIELD_DIM 1
#define _REAL double
#include "TMPL/strided_reads_gpu.hpp.tmpl"

#define _FIELD_NAME clover_ldl
#define _FIELD_TYPE ldl_field
#define _SITE_TYPE ldl_t
#define _FIELD_DIM 1
#define _REAL double
#include "TMPL/strided_reads_gpu.hpp.tmpl"

#define _FIELD_NAME clover_term
#define _FIELD_TYPE suNfc_field
#define _SITE_TYPE suNfc
#define _FIELD_DIM 4
#define _REAL double
#include "TMPL/strided_reads_gpu.hpp.tmpl"

#define _FIELD_NAME clover_force
#define _FIELD_TYPE suNf_field
#define _SITE_TYPE suNf
#define _FIELD_DIM 6
#define _REAL double
#include "TMPL/strided_reads_gpu.hpp.tmpl"

#define _FIELD_NAME staple_field
#define _FIELD_TYPE suNg_field
#define _SITE_TYPE suNg
#define _FIELD_DIM 3
#define _REAL double
#include "TMPL/strided_reads_gpu.hpp.tmpl"

#ifdef __cplusplus
    }
#endif

#endif