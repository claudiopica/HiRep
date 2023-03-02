/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *   
* All rights reserved.                                                      * 
\***************************************************************************/

#include "random.h"
#include "libhr_core.h"
#include "Inverters/linear_algebra.h"
#include "Utils/boundary_conditions.h"
#include "memory.h"

void gaussian_spinor_field(spinor_field *s) {
#ifdef WITH_FUSE_MASTER_FOR
    _FUSE_MASTER_FOR(s->type, ix) {
        _FUSE_IDX(s->type, ix);
#else
    _MASTER_FOR(s->type, ix) {
#endif
        suNf_spinor *r = _FIELD_AT(s, ix);
        gauss((double *)(r), sizeof(suNf_spinor) / sizeof(double));
    }
    const double c1 = 1. / sqrt(2.);
    spinor_field_mul_f_cpu(s, c1, s);

#ifdef WITH_GPU
    copy_to_gpu_spinor_field_f(s);
#endif

    apply_BCs_on_spinor_field(s);
}

void gaussian_spinor_field_flt(spinor_field_flt *s) {
#ifdef WITH_FUSE_MASTER_FOR
    _FUSE_MASTER_FOR(s->type, ix) {
        _FUSE_IDX(s->type, ix);
#else
    _MASTER_FOR(s->type, ix) {
#endif
        suNf_spinor_flt *r = _FIELD_AT(s, ix);
        gauss_flt((float *)(r), sizeof(suNf_spinor_flt) / sizeof(float));
    }
    const float c1 = (float)(1. / sqrt(2.));
    spinor_field_mul_f_flt_cpu(s, c1, s);

#ifdef WITH_GPU
    copy_to_gpu_spinor_field_f_flt(s);
#endif

    apply_BCs_on_spinor_field_flt(s);
}

void z2_spinor_field(spinor_field *s) {
#ifdef WITH_FUSE_MASTER_FOR
    _FUSE_MASTER_FOR(s->type, ix) {
        _FUSE_IDX(s->type, ix);
#else
    _MASTER_FOR(s->type, ix) {
#endif
        suNf_spinor *r = _FIELD_AT(s, ix);
        ranz2((double *)(r), sizeof(suNf_spinor) / sizeof(double));
    }
    apply_BCs_on_spinor_field(s);
}
