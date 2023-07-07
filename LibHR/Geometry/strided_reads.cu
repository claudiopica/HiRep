/***************************************************************************\
* Copyright (c) 2022, Sofie Martins                                         *   
* All rights reserved.                                                      * 
\***************************************************************************/

#ifdef WITH_GPU

#include "geometry.h"
#include "libhr_core.h"

__device__ void read_clover(hr_complex *c, const ldl_t *in, int ix, int dn, int j) {
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

__device__ void write_clover(hr_complex *c, ldl_t *out, int ix, int dn, int j) {
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

__device__ void read_force(hr_complex *c, const suNf *in, int ix, int comp, int i, int j) {
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

__device__ void write_force(hr_complex *c, suNf *out, int ix, int comp, int i, int j) {
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

__device__ void read_clover_term_comp(hr_complex *c, const suNfc *in, int ix, int comp, int i, int j) {
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

__device__ void write_clover_term_comp(hr_complex *c, suNfc *out, int ix, int comp, int i, int j) {
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

// Note: C functions are less flexible than C++ functions
// they are basically only needed for the geometry convert

#define _FIELD_NAME spinor_field_f
#define _FIELD_TYPE spinor_field
#define _SITE_TYPE suNf_spinor
#define _FIELD_DIM 1
#define _REAL double
#include "TMPL/strided_reads.cu.tmpl"

#define _FIELD_NAME spinor_field_f_flt
#define _FIELD_TYPE spinor_field_flt
#define _SITE_TYPE suNf_spinor_flt
#define _FIELD_DIM 1
#define _REAL float
#include "TMPL/strided_reads.cu.tmpl"

#define _FIELD_NAME sfield
#define _FIELD_TYPE scalar_field
#define _SITE_TYPE double
#define _FIELD_DIM 1
#define _REAL double
#include "TMPL/strided_reads.cu.tmpl"

#define _FIELD_NAME gfield
#define _FIELD_TYPE suNg_field
#define _SITE_TYPE suNg
#define _FIELD_DIM 4
#define _REAL double
#include "TMPL/strided_reads.cu.tmpl"

#define _FIELD_NAME gfield_flt
#define _FIELD_TYPE suNg_field_flt
#define _SITE_TYPE suNg_flt
#define _FIELD_DIM 4
#define _REAL float
#include "TMPL/strided_reads.cu.tmpl"

#define _FIELD_NAME gfield_f
#define _FIELD_TYPE suNf_field
#define _SITE_TYPE suNf
#define _FIELD_DIM 4
#define _REAL double
#include "TMPL/strided_reads.cu.tmpl"

#define _FIELD_NAME gfield_f_flt
#define _FIELD_TYPE suNf_field_flt
#define _SITE_TYPE suNf_flt
#define _FIELD_DIM 4
#define _REAL float
#include "TMPL/strided_reads.cu.tmpl"

#define _FIELD_NAME suNg_scalar_field
#define _FIELD_TYPE suNg_scalar_field
#define _SITE_TYPE suNg_vector
#define _FIELD_DIM 1
#define _REAL double
#include "TMPL/strided_reads.cu.tmpl"

#define _FIELD_NAME avfield
#define _FIELD_TYPE suNg_av_field
#define _SITE_TYPE suNg_algebra_vector
#define _FIELD_DIM 4
#define _REAL double
#include "TMPL/strided_reads.cu.tmpl"

#define _FIELD_NAME gtransf
#define _FIELD_TYPE suNg_field
#define _SITE_TYPE suNg
#define _FIELD_DIM 1
#define _REAL double
#include "TMPL/strided_reads.cu.tmpl"

#define _FIELD_NAME clover_ldl
#define _FIELD_TYPE ldl_field
#define _SITE_TYPE ldl_t
#define _FIELD_DIM 1
#define _REAL double
#include "TMPL/strided_reads.cu.tmpl"

#define _FIELD_NAME clover_term
#define _FIELD_TYPE suNfc_field
#define _SITE_TYPE suNfc
#define _FIELD_DIM 4
#define _REAL double
#include "TMPL/strided_reads.cu.tmpl"

#define _FIELD_NAME clover_force
#define _FIELD_TYPE suNf_field
#define _SITE_TYPE suNf
#define _FIELD_DIM 6
#define _REAL double
#include "TMPL/strided_reads.cu.tmpl"

#define _FIELD_NAME staple_field
#define _FIELD_TYPE suNg_field
#define _SITE_TYPE suNg
#define _FIELD_DIM 3
#define _REAL double
#include "TMPL/strided_reads.cu.tmpl"

#endif