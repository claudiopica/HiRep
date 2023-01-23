// These reading functions do not work with real representations yet (SAM)

#ifndef STRIDED_READS_GPU_H
#define STRIDED_READS_GPU_H

#define THREADSIZE 32

#include "libhr_core.h"
#include "geometry.h"

template<typename COMPLEX, typename FIELD_TYPE, typename SITE_TYPE>
__device__ void read_gpu(int stride, SITE_TYPE *s, const FIELD_TYPE *in, int ix, int comp, int dim) {
    const int field_dim = sizeof(FIELD_TYPE)/sizeof(COMPLEX);
    const int n_components = sizeof(SITE_TYPE)/sizeof(COMPLEX);
    #ifdef FIXED_STRIDE
        int iz = ((ix / THREADSIZE) * THREADSIZE) * dim * field_dim  + (ix % THREADSIZE) + ((comp)*n_components)*(THREADSIZE);
    #else
        int iz = ix ((comp)*n_components)*(THREADSIZE);
    #endif
    COMPLEX* in_cpx = (COMPLEX*)in;
    COMPLEX* in_comp_cpx = (COMPLEX*)s;
    for (int i = 0; i < n_components; ++i) {
         in_comp_cpx[i] = in_cpx[iz];
         iz+=THREADSIZE;
    }
}

template <typename COMPLEX, typename FIELD_TYPE, typename SITE_TYPE>
__device__ void write_gpu(int stride, SITE_TYPE *s, FIELD_TYPE *out, int ix, int comp, int dim) {
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