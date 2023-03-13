/***************************************************************************\
* Copyright (c) 2008-2023, Claudio Pica, Sofie Martins                      *   
* All rights reserved.                                                      * 
\***************************************************************************/

/**
 * @file single_double_utils.c
 * @brief Functions for conversion from single to double precision and
 *        vice versa. 
 */

#include "utils.h"
#include "libhr_core.h"
#include "memory.h"

void assign_ud2u_cpu(void) {
    if (u_gauge_flt != NULL) {
        double *d;
        float *f;

        d = (double *)(u_gauge->ptr);
        f = (float *)(u_gauge_flt->ptr);

        const int ndoubles_per_site = sizeof(suNg) / sizeof(double);
        const size_t ndoubles = 4 * glattice.gsize_gauge * ndoubles_per_site;

        _OMP_PRAGMA(_omp_parallel)
        _OMP_PRAGMA(_omp_for)
        for (int i = 0; i < ndoubles; i++) {
            *(f + i) = (float)(*(d + i));
        }
    } else {
        error(0, 1, __func__,
              "Trying to assign to single precision field "
              "that was not previously allocated. Try to "
              "recompile with DPHI_FLT.\n");
    }
}

void assign_u2ud_cpu(void) {
    if (u_gauge_flt != NULL) {
        double *d;
        float *f;

        d = (double *)(u_gauge->ptr);
        f = (float *)(u_gauge_flt->ptr);

        const int nfloats_per_site = sizeof(suNg_flt) / sizeof(float);
        const size_t nfloats = 4 * glattice.gsize_gauge * nfloats_per_site;

        _OMP_PRAGMA(_omp_parallel)
        _OMP_PRAGMA(_omp_for)
        for (int i = 0; i < nfloats; i++) {
            *(f + i) = (float)(*(d + i));
        }
    } else {
        error(0, 1, __func__,
              "Trying to assign values from a single precision "
              "field that was not previously allocated. "
              "Try to recompile with DPHI_FLT.\n");
    }
}

void assign_ud2u_f_cpu(void) {
    if (u_gauge_f_flt != NULL) {
        double *d;
        float *f;

        d = (double *)(u_gauge_f->ptr);
        f = (float *)(u_gauge_f_flt->ptr);

        const int ndoubles_per_site = sizeof(suNf) / sizeof(double);
        const size_t ndoubles = 4 * glattice.gsize_gauge * ndoubles_per_site;

        _OMP_PRAGMA(_omp_parallel)
        _OMP_PRAGMA(_omp_for)
        for (int i = 0; i < ndoubles; i++) {
            *(f + i) = (float)(*(d + i));
        }
    } else {
        error(0, 1, __func__,
              "Trying to assign to single precision field "
              "that was not previously allocated. Try to "
              "recompile with DPHI_FLT.\n");
    }
}

void assign_u2ud_f_cpu(void) {
    if (u_gauge_f_flt != NULL) {
        double *d;
        float *f;

        d = (double *)(u_gauge_f->ptr);
        f = (float *)(u_gauge_f_flt->ptr);

        const int nfloats_per_site = sizeof(suNf_flt) / sizeof(float);
        const size_t nfloats = 4 * glattice.gsize_gauge * nfloats_per_site;

        _OMP_PRAGMA(_omp_parallel)
        _OMP_PRAGMA(_omp_for)
        for (int i = 0; i < nfloats; i++) {
            *(d + i) = (double)(*(f + i));
        }
    } else {
        error(0, 1, __func__,
              "Trying to assign values from a single precision "
              "field that was not previously allocated. "
              "Try to recompile with DPHI_FLT.\n");
    }
}

void assign_s2sd_cpu(spinor_field *out, spinor_field_flt *in) {
    _TWO_SPINORS_FOR(out, in) {
        double *o = (double *)_SPINOR_PTR(out);
        float *i = (float *)_SPINOR_PTR(in);
        for (int n = 0; n < (8 * NF); n++) {
            *(o++) = (double)*(i++);
        }
    }
}

void assign_sd2s_cpu(spinor_field_flt *out, spinor_field *in) {
    _TWO_SPINORS_FOR(out, in) {
        float *o = (float *)_SPINOR_PTR(out);
        double *i = (double *)_SPINOR_PTR(in);
        for (int n = 0; n < (8 * NF); n++) {
            *(o++) = (float)*(i++);
        }
    }
}

void add_assign_s2sd_cpu(spinor_field *out, spinor_field_flt *in) {
    _TWO_SPINORS_FOR(out, in) {
        double *o = (double *)_SPINOR_PTR(out);
        float *i = (float *)_SPINOR_PTR(in);
        for (int n = 0; n < (8 * NF); n++) {
            *(o++) += (double)*(i++);
        }
    }
}

void add_assign_sd2s_cpu(spinor_field_flt *out, spinor_field *in) {
    _TWO_SPINORS_FOR(out, in) {
        float *o = (float *)_SPINOR_PTR(out);
        double *i = (double *)_SPINOR_PTR(in);
        for (int n = 0; n < (8 * NF); n++) {
            *(o++) += (float)*(i++);
        }
    }
}

#ifndef WITH_GPU
void (*assign_ud2u)(void) = assign_ud2u_cpu;
void (*assign_u2ud)(void) = assign_u2ud_cpu;
void (*assign_ud2u_f)(void) = assign_ud2u_f_cpu;
void (*assign_u2ud_f)(void) = assign_u2ud_f_cpu;
void (*assign_s2sd)(spinor_field *out, spinor_field_flt *in) = assign_s2sd_cpu;
void (*assign_sd2s)(spinor_field_flt *out, spinor_field *in) = assign_sd2s_cpu;
void (*add_assign_s2sd)(spinor_field *out, spinor_field_flt *in) = add_assign_s2sd_cpu;
void (*add_assign_sd2s)(spinor_field_flt *out, spinor_field *in) = add_assign_sd2s_cpu;
#endif
