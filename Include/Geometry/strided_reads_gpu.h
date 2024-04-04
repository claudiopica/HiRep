/***************************************************************************\
* Copyright (c) 2023, Sofie Martins                                         *
* All rights reserved.                                                      *
\***************************************************************************/

/// Headerfile for:
/// - strided_reads.cu

#ifndef STRIDED_READS_GPU_HPP
#define STRIDED_READS_GPU_HPP

#ifdef FIXED_STRIDE
#define THREADSIZE 32
#else
#define THREADSIZE 1
#endif

//#include "libhr_core.h"
//#include "geometry.h"

enum DIRECTION { UP = 0, DOWN = 1 };

#ifdef __cplusplus

template <typename REAL, typename FIELD_TYPE, typename SITE_TYPE>
__host__ __device__ __forceinline__ void read_gpu(int stride, SITE_TYPE *s, const FIELD_TYPE *in, size_t ix, int comp,
                                                  int dim) {
    const int field_dim = sizeof(FIELD_TYPE) / sizeof(REAL);
    const int n_components = sizeof(SITE_TYPE) / sizeof(REAL);
#ifdef FIXED_STRIDE
    size_t iz = ((ix / THREADSIZE) * THREADSIZE) * dim * field_dim + (ix % THREADSIZE) + ((comp)*n_components) * (THREADSIZE);
    const int _stride = THREADSIZE;
#else
    size_t iz = ix + ((comp)*n_components) * (THREADSIZE);
    const int _stride = stride;
#endif
    REAL *in_cpx = (REAL *)in;
    REAL *in_comp_cpx = (REAL *)s;
    for (int i = 0; i < n_components; ++i) {
        in_comp_cpx[i] = in_cpx[iz];
        iz += _stride;
    }
}

template <typename REAL, typename FIELD_TYPE, typename SITE_TYPE>
__host__ __device__ __forceinline__ void write_gpu(int stride, SITE_TYPE *s, FIELD_TYPE *out, size_t ix, int comp, int dim) {
    const int field_dim = sizeof(FIELD_TYPE) / sizeof(REAL);
    const int n_components = sizeof(SITE_TYPE) / sizeof(REAL);
#ifdef FIXED_STRIDE
    size_t iz = ((ix / THREADSIZE) * THREADSIZE) * dim * field_dim + (ix % THREADSIZE) + (comp)*n_components * (THREADSIZE);
    const int _stride = THREADSIZE;
#else
    size_t iz = ix + ((comp)*n_components) * (THREADSIZE);
    const int _stride = stride;
#endif
    REAL *out_cpx = (REAL *)out;
    REAL *out_comp_cpx = (REAL *)s;
    for (int i = 0; i < n_components; ++i) {
        out_cpx[iz] = out_comp_cpx[i];
        iz += _stride;
    }
}

template <typename REAL, typename FIELD_TYPE, typename SITE_TYPE>
__host__ __device__ __forceinline__ void read_assign_gpu(int stride, SITE_TYPE *s, const FIELD_TYPE *in, int ix, int comp,
                                                         int dim) {
    const int field_dim = sizeof(FIELD_TYPE) / sizeof(REAL);
    const int n_components = sizeof(SITE_TYPE) / sizeof(REAL);
#ifdef FIXED_STRIDE
    int iz = ((ix / THREADSIZE) * THREADSIZE) * dim * field_dim + (ix % THREADSIZE) + ((comp)*n_components) * (THREADSIZE);
    const int _stride = THREADSIZE;
#else
    int iz = ix + ((comp)*n_components) * (THREADSIZE);
    const int _stride = stride;
#endif
    REAL *in_cpx = (REAL *)in;
    REAL *in_comp_cpx = (REAL *)s;
    for (int i = 0; i < n_components; ++i) {
        in_comp_cpx[i] += in_cpx[iz];
        iz += _stride;
    }
}

template <typename REAL, typename FIELD_TYPE, typename SITE_TYPE>
__host__ __device__ __forceinline__ void write_assign_gpu(int stride, SITE_TYPE *s, FIELD_TYPE *out, int ix, int comp,
                                                          int dim) {
    const int field_dim = sizeof(FIELD_TYPE) / sizeof(REAL);
    const int n_components = sizeof(SITE_TYPE) / sizeof(REAL);
#ifdef FIXED_STRIDE
    int iz = ((ix / THREADSIZE) * THREADSIZE) * dim * field_dim + (ix % THREADSIZE) + (comp)*n_components * (THREADSIZE);
    const int _stride = THREADSIZE;
#else
    int iz = ix + ((comp)*n_components) * (THREADSIZE);
    const int _stride = stride;
#endif
    REAL *out_cpx = (REAL *)out;
    REAL *out_comp_cpx = (REAL *)s;
    for (int i = 0; i < n_components; ++i) {
        out_cpx[iz] += out_comp_cpx[i];
        iz += _stride;
    }
}

__device__ __forceinline__ void read_clover(hr_complex *c, const ldl_t *in, int ix, int dn, int j) {
    // Number of doubles per site
    const int n_components = sizeof(ldl_t) / sizeof(double);

    // offset: Round down to next lowest multiple of 32,
    // find the offset at double level
    int iz = ((ix / THREADSIZE) * THREADSIZE) * n_components;

    // find the first double in the array of the given index
    // components can now be found at strides of 32 away
    iz += ix % THREADSIZE;

    // Move by j hr complexes (containing 2 doubles)
    // also offset by half the components if you want to read down
    iz += 2 * j * THREADSIZE + dn * n_components * THREADSIZE / 2;

    double *in_cpx = (double *)in;
    double *in_comp_cpx = (double *)c;
    for (int i = 0; i < 2; ++i) {
        in_comp_cpx[i] = in_cpx[iz];
        iz += THREADSIZE;
    }
}

__device__ __forceinline__ void write_clover(hr_complex *c, ldl_t *out, int ix, int dn, int j) {
    const int n_components = sizeof(ldl_t) / sizeof(double);
    int iz = ((ix / THREADSIZE) * THREADSIZE) * n_components;
    iz += ix % THREADSIZE;
    iz += 2 * j * THREADSIZE + dn * n_components * THREADSIZE / 2;
    double *out_cpx = (double *)out;
    double *out_comp_cpx = (double *)c;
    for (int i = 0; i < 2; ++i) {
        out_cpx[iz] = out_comp_cpx[i];
        iz += THREADSIZE;
    }
}

__device__ __forceinline__ void read_force(hr_complex *c, const suNf *in, int ix, int comp, int i, int j) {
    const int n_components = sizeof(suNf) / sizeof(double);
    int iz = ((ix / THREADSIZE) * THREADSIZE) * n_components * 6;
    iz += ix % THREADSIZE;
    // Move by j hr complexes (containing 2 doubles)
    // also offset by half the components if you want to read down
    iz += 2 * (i * NF + j) * THREADSIZE + comp * n_components * THREADSIZE;

    double *in_cpx = (double *)in;
    double *in_comp_cpx = (double *)c;
    for (int i = 0; i < 2; ++i) {
        in_comp_cpx[i] = in_cpx[iz];
        iz += THREADSIZE;
    }
}

__device__ __forceinline__ void write_force(hr_complex *c, suNf *out, int ix, int comp, int i, int j) {
    const int n_components = sizeof(suNf) / sizeof(double);
    int iz = ((ix / THREADSIZE) * THREADSIZE) * n_components * 6;
    iz += ix % THREADSIZE;
    iz += 2 * (i * NF + j) * THREADSIZE + comp * n_components * THREADSIZE;
    double *out_cpx = (double *)out;
    double *out_comp_cpx = (double *)c;
    for (int i = 0; i < 2; ++i) {
        out_cpx[iz] = out_comp_cpx[i];
        iz += THREADSIZE;
    }
}

__device__ __forceinline__ void read_clover_term_comp(hr_complex *c, const suNfc *in, int ix, int comp, int i, int j) {
    const int n_components = sizeof(suNfc) / sizeof(double);
    int iz = ((ix / THREADSIZE) * THREADSIZE) * 4 * n_components;
    iz += ix % (THREADSIZE);
    iz += comp * n_components * (THREADSIZE);
    iz += (i * NF + j) * 2 * (THREADSIZE);
    double *in_cpx = (double *)in;
    double *in_comp_cpx = (double *)c;
    for (int i = 0; i < 2; ++i) {
        in_comp_cpx[i] = in_cpx[iz];
        iz += THREADSIZE;
    }
}

__device__ __forceinline__ void write_clover_term_comp(hr_complex *c, suNfc *out, int ix, int comp, int i, int j) {
    const int n_components = sizeof(suNfc) / sizeof(double);
    int iz = ((ix / THREADSIZE) * THREADSIZE) * 4 * n_components;
    iz += ix % (THREADSIZE);
    iz += comp * n_components * (THREADSIZE);
    iz += (i * NF + j) * 2 * (THREADSIZE);
    double *out_cpx = (double *)out;
    double *out_comp_cpx = (double *)c;
    for (int i = 0; i < 2; ++i) {
        out_cpx[iz] = out_comp_cpx[i];
        iz += THREADSIZE;
    }
    template <typename REAL, typename FIELD_TYPE, typename SITE_TYPE>
    __host__ __device__ void read_sub_assign_gpu(int stride, SITE_TYPE *s, const FIELD_TYPE *in, int ix, int comp, int dim) {
        const int field_dim = sizeof(FIELD_TYPE) / sizeof(REAL);
        const int n_components = sizeof(SITE_TYPE) / sizeof(REAL);
        template <typename REAL, typename FIELD_TYPE, typename SITE_TYPE>
        __host__ __device__ __forceinline__ void read_sub_assign_gpu(int stride, SITE_TYPE *s, const FIELD_TYPE *in, int ix,
                                                                     int comp, int dim) {
            const int field_dim = sizeof(FIELD_TYPE) / sizeof(REAL);
            const int n_components = sizeof(SITE_TYPE) / sizeof(REAL);
#ifdef FIXED_STRIDE
            int iz =
                ((ix / THREADSIZE) * THREADSIZE) * dim * field_dim + (ix % THREADSIZE) + ((comp)*n_components) * (THREADSIZE);
            const int _stride = THREADSIZE;
#else
            int iz = ix + ((comp)*n_components) * (THREADSIZE);
            const int _stride = stride;
#endif
            REAL *in_cpx = (REAL *)in;
            REAL *in_comp_cpx = (REAL *)s;
            for (int i = 0; i < n_components; ++i) {
                in_comp_cpx[i] -= in_cpx[iz];
                iz += _stride;
            }
        }

        template <typename REAL, typename VECTOR_TYPE, typename SITE_TYPE>
        __device__ __forceinline__ void in_spinor_field(VECTOR_TYPE * v, SITE_TYPE * in, int iy, int comp) {
            read_gpu<REAL>(0, v, in, iy, comp, 1);
        }

        template <typename REAL, typename GAUGE_TYPE>
        __device__ __forceinline__ void in_gauge_field(GAUGE_TYPE * u, const GAUGE_TYPE *in, int ix, int iy, int comp,
                                                       int dir) {
            if (dir == UP) {
                read_gpu<REAL>(0, u, in, ix, comp, 4);
            } else if (dir == DOWN) {
                read_gpu<REAL>(0, u, in, iy, comp, 4);
            }

            template <typename REAL, typename SITE_TYPE>
            __device__ __forceinline__ void write_out_spinor_field(SITE_TYPE * r, SITE_TYPE * in, int ix) {
                write_gpu<REAL>(0, r, in, ix, 0, 1);
            }

            template <typename REAL, typename FIELD_TYPE, typename SITE_TYPE>
            __host__ __device__ __forceinline__ void write_assign_gpu(int stride, SITE_TYPE *s, FIELD_TYPE *out, int ix,
                                                                      int comp, int dim) {
                const int field_dim = sizeof(FIELD_TYPE) / sizeof(REAL);
                const int n_components = sizeof(SITE_TYPE) / sizeof(REAL);
#ifdef FIXED_STRIDE
                int iz =
                    ((ix / THREADSIZE) * THREADSIZE) * dim * field_dim + (ix % THREADSIZE) + (comp)*n_components * (THREADSIZE);
                const int _stride = THREADSIZE;
#else
                int iz = ix + ((comp)*n_components) * (THREADSIZE);
                const int _stride = stride;
#endif
                REAL *out_cpx = (REAL *)out;
                REAL *out_comp_cpx = (REAL *)s;
                for (int i = 0; i < n_components; ++i) {
                    out_cpx[iz] += out_comp_cpx[i];
                    iz += _stride;
                }
            }

            __device__ __forceinline__ void read_clover(hr_complex * c, const ldl_t *in, int ix, int dn, int j) {
                // Number of doubles per site
                const int n_components = sizeof(ldl_t) / sizeof(double);

                // offset: Round down to next lowest multiple of 32,
                // find the offset at double level
                int iz = ((ix / THREADSIZE) * THREADSIZE) * n_components;

                // find the first double in the array of the given index
                // components can now be found at strides of 32 away
                iz += ix % THREADSIZE;

                // Move by j hr complexes (containing 2 doubles)
                // also offset by half the components if you want to read down
                iz += 2 * j * THREADSIZE + dn * n_components * THREADSIZE / 2;

                double *in_cpx = (double *)in;
                double *in_comp_cpx = (double *)c;
                for (int i = 0; i < 2; ++i) {
                    in_comp_cpx[i] = in_cpx[iz];
                    iz += THREADSIZE;
                }
            }

            __device__ __forceinline__ void write_clover(hr_complex * c, ldl_t * out, int ix, int dn, int j) {
                const int n_components = sizeof(ldl_t) / sizeof(double);
                int iz = ((ix / THREADSIZE) * THREADSIZE) * n_components;
                iz += ix % THREADSIZE;
                iz += 2 * j * THREADSIZE + dn * n_components * THREADSIZE / 2;
                double *out_cpx = (double *)out;
                double *out_comp_cpx = (double *)c;
                for (int i = 0; i < 2; ++i) {
                    out_cpx[iz] = out_comp_cpx[i];
                    iz += THREADSIZE;
                }
            }

            __device__ __forceinline__ void read_force(hr_complex * c, const suNf *in, int ix, int comp, int i, int j) {
                const int n_components = sizeof(suNf) / sizeof(double);
                int iz = ((ix / THREADSIZE) * THREADSIZE) * n_components * 6;
                iz += ix % THREADSIZE;
                // Move by j hr complexes (containing 2 doubles)
                // also offset by half the components if you want to read down
                iz += 2 * (i * NF + j) * THREADSIZE + comp * n_components * THREADSIZE;

                double *in_cpx = (double *)in;
                double *in_comp_cpx = (double *)c;
                for (int i = 0; i < 2; ++i) {
                    in_comp_cpx[i] = in_cpx[iz];
                    iz += THREADSIZE;
                }
            }

            __device__ __forceinline__ void write_force(hr_complex * c, suNf * out, int ix, int comp, int i, int j) {
                const int n_components = sizeof(suNf) / sizeof(double);
                int iz = ((ix / THREADSIZE) * THREADSIZE) * n_components * 6;
                iz += ix % THREADSIZE;
                iz += 2 * (i * NF + j) * THREADSIZE + comp * n_components * THREADSIZE;
                double *out_cpx = (double *)out;
                double *out_comp_cpx = (double *)c;
                for (int i = 0; i < 2; ++i) {
                    out_cpx[iz] = out_comp_cpx[i];
                    iz += THREADSIZE;
                }
            }

            __device__ __forceinline__ void read_clover_term_comp(hr_complex * c, const suNfc *in, int ix, int comp, int i,
                                                                  int j) {
                const int n_components = sizeof(suNfc) / sizeof(double);
                int iz = ((ix / THREADSIZE) * THREADSIZE) * 4 * n_components;
                iz += ix % (THREADSIZE);
                iz += comp * n_components * (THREADSIZE);
                iz += (i * NF + j) * 2 * (THREADSIZE);
                double *in_cpx = (double *)in;
                double *in_comp_cpx = (double *)c;
                for (int i = 0; i < 2; ++i) {
                    in_comp_cpx[i] = in_cpx[iz];
                    iz += THREADSIZE;
                }
            }

            __device__ __forceinline__ void write_clover_term_comp(hr_complex * c, suNfc * out, int ix, int comp, int i,
                                                                   int j) {
                const int n_components = sizeof(suNfc) / sizeof(double);
                int iz = ((ix / THREADSIZE) * THREADSIZE) * 4 * n_components;
                iz += ix % (THREADSIZE);
                iz += comp * n_components * (THREADSIZE);
                iz += (i * NF + j) * 2 * (THREADSIZE);
                double *out_cpx = (double *)out;
                double *out_comp_cpx = (double *)c;
                for (int i = 0; i < 2; ++i) {
                    out_cpx[iz] = out_comp_cpx[i];
                    iz += THREADSIZE;
                }
            }

#endif

#ifdef __cplusplus
            extern "C" {
#endif

#define _FIELD_TYPE spinor_field
#define _SITE_TYPE suNf_spinor
#define _FIELD_DIM 1
#define _REAL double
#include "TMPL/strided_reads_gpu.h.tmpl"

#define _FIELD_TYPE spinor_field_flt
#define _SITE_TYPE suNf_spinor_flt
#define _FIELD_DIM 1
#define _REAL float
#include "TMPL/strided_reads_gpu.h.tmpl"

#define _FIELD_TYPE scalar_field
#define _SITE_TYPE double
#define _FIELD_DIM 1
#define _REAL double
#include "TMPL/strided_reads_gpu.h.tmpl"

#define _FIELD_TYPE suNg_field
#define _SITE_TYPE suNg
#define _FIELD_DIM 4
#define _REAL double
#include "TMPL/strided_reads_gpu.h.tmpl"

#define _FIELD_TYPE suNg_field_flt
#define _SITE_TYPE suNg_flt
#define _FIELD_DIM 4
#define _REAL float
#include "TMPL/strided_reads_gpu.h.tmpl"

#define _FIELD_TYPE suNf_field
#define _SITE_TYPE suNf
#define _FIELD_DIM 4
#define _REAL double
#include "TMPL/strided_reads_gpu.h.tmpl"

#define _FIELD_TYPE suNf_field_flt
#define _SITE_TYPE suNf_flt
#define _FIELD_DIM 4
#define _REAL float
#include "TMPL/strided_reads_gpu.h.tmpl"

#define _FIELD_TYPE suNg_scalar_field
#define _SITE_TYPE suNg_vector
#define _FIELD_DIM 1
#define _REAL double
#include "TMPL/strided_reads_gpu.h.tmpl"

#define _FIELD_TYPE suNg_av_field
#define _SITE_TYPE suNg_algebra_vector
#define _FIELD_DIM 4
#define _REAL double
#include "TMPL/strided_reads_gpu.h.tmpl"

#define _FIELD_TYPE gtransf
#define _SITE_TYPE suNg
#define _FIELD_DIM 1
#define _REAL double
#include "TMPL/strided_reads_gpu.h.tmpl"

#define _FIELD_TYPE ldl_field
#define _SITE_TYPE ldl_t
#define _FIELD_DIM 1
#define _REAL double
#include "TMPL/strided_reads_gpu.h.tmpl"

#define _FIELD_TYPE clover_term
#define _SITE_TYPE suNfc
#define _FIELD_DIM 4
#define _REAL double
#include "TMPL/strided_reads_gpu.h.tmpl"

#define _FIELD_TYPE clover_force
#define _SITE_TYPE suNf
#define _FIELD_DIM 6
#define _REAL double
#include "TMPL/strided_reads_gpu.h.tmpl"

#define _FIELD_TYPE staple_field
#define _SITE_TYPE suNg
#define _FIELD_DIM 3
#define _REAL double
#include "TMPL/strided_reads_gpu.h.tmpl"

#ifdef __cplusplus
            }
#endif

#endif