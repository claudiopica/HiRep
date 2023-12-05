/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *
* All rights reserved.                                                      *
\***************************************************************************/

#include "inverters.h"
#include "libhr_core.h"
#include "Utils/generics.h"

// Linear Algebra functions are generic
// They are parametrized over the input types for double/single precision

// double precision
#define _FIELD_TYPE spinor_field
#define _SITE_TYPE suNf_spinor
#define _REAL double
#define _COMPLEX hr_complex
#define _FIELD_DIM 1
#include "TMPL/linear_algebra_reduction.c.tmpl"
#include "TMPL/linear_algebra_base_operations.c.tmpl"
#include "TMPL/linear_algebra_lc.c.tmpl"
#include "TMPL/linear_algebra_gamma.c.tmpl"
#include "TMPL/linear_algebra_base.c.tmpl" // This one has to go last

// single precision
#define _FIELD_TYPE spinor_field_flt
#define _SITE_TYPE suNf_spinor_flt
#define _REAL float
#define _COMPLEX hr_complex_flt
#define _FIELD_DIM 1
#include "TMPL/linear_algebra_reduction.c.tmpl"
#include "TMPL/linear_algebra_base_operations.c.tmpl"
#include "TMPL/linear_algebra_lc.c.tmpl"
#include "TMPL/linear_algebra_gamma.c.tmpl"
#include "TMPL/linear_algebra_base.c.tmpl"

// double precision
#define _FIELD_TYPE scalar_field
#define _SITE_TYPE double
#define _REAL double
#define _COMPLEX hr_complex
#define _FIELD_DIM 1
#include "TMPL/linear_algebra_reduction.c.tmpl"
#include "TMPL/linear_algebra_base_operations.c.tmpl"
#include "TMPL/linear_algebra_base.c.tmpl"

#define _FIELD_TYPE suNg_field
#define _SITE_TYPE suNg
#define _REAL double
#define _COMPLEX hr_complex
#define _FIELD_DIM 4
#include "TMPL/linear_algebra_reduction.c.tmpl"
#include "TMPL/linear_algebra_base_operations.c.tmpl"
#include "TMPL/linear_algebra_base.c.tmpl"

#define _FIELD_TYPE suNf_field
#define _SITE_TYPE suNf
#define _REAL double
#define _COMPLEX hr_complex
#define _FIELD_DIM 4
#include "TMPL/linear_algebra_reduction.c.tmpl"
#include "TMPL/linear_algebra_base_operations.c.tmpl"
#include "TMPL/linear_algebra_base.c.tmpl"

#define _FIELD_TYPE suNfc_field
#define _SITE_TYPE suNfc
#define _REAL double
#define _COMPLEX hr_complex
#define _FIELD_DIM 4
#include "TMPL/linear_algebra_reduction.c.tmpl"
#include "TMPL/linear_algebra_base_operations.c.tmpl"
#include "TMPL/linear_algebra_base.c.tmpl"

#define _FIELD_TYPE suNg_field_flt
#define _SITE_TYPE suNg_flt
#define _REAL float
#define _COMPLEX hr_complex_flt
#define _FIELD_DIM 4
#include "TMPL/linear_algebra_reduction.c.tmpl"
#include "TMPL/linear_algebra_base_operations.c.tmpl"
#include "TMPL/linear_algebra_base.c.tmpl"

#define _FIELD_TYPE suNf_field_flt
#define _SITE_TYPE suNf_flt
#define _REAL float
#define _COMPLEX hr_complex_flt
#define _FIELD_DIM 4
#include "TMPL/linear_algebra_reduction.c.tmpl"
#include "TMPL/linear_algebra_base_operations.c.tmpl"
#include "TMPL/linear_algebra_base.c.tmpl"

#define _FIELD_TYPE suNg_scalar_field
#define _SITE_TYPE suNg_vector
#define _REAL double
#define _COMPLEX hr_complex
#define _FIELD_DIM 1
#include "TMPL/linear_algebra_reduction.c.tmpl"
#include "TMPL/linear_algebra_base_operations.c.tmpl"
#include "TMPL/linear_algebra_base.c.tmpl"

#define _FIELD_TYPE suNg_av_field
#define _SITE_TYPE suNg_algebra_vector
#define _REAL double
#define _COMPLEX hr_complex
#define _FIELD_DIM 4
#include "TMPL/linear_algebra_reduction.c.tmpl"
#include "TMPL/linear_algebra_base_operations.c.tmpl"
#include "TMPL/linear_algebra_base.c.tmpl"

#define _FIELD_TYPE gtransf
#define _SITE_TYPE suNg
#define _REAL double
#define _COMPLEX hr_complex
#define _FIELD_DIM 1
#include "TMPL/linear_algebra_reduction.c.tmpl"
#include "TMPL/linear_algebra_base_operations.c.tmpl"
#include "TMPL/linear_algebra_base.c.tmpl"

#define _FIELD_TYPE ldl_field
#define _SITE_TYPE ldl_t
#define _REAL double
#define _COMPLEX hr_complex
#define _FIELD_DIM 1
#include "TMPL/linear_algebra_base.c.tmpl"

#define _FIELD_TYPE clover_term
#define _SITE_TYPE suNfc
#define _REAL double
#define _COMPLEX hr_complex
#define _FIELD_DIM 4
#include "TMPL/linear_algebra_reduction.c.tmpl"
#include "TMPL/linear_algebra_base_operations.c.tmpl"
#include "TMPL/linear_algebra_base.c.tmpl"

#define _FIELD_TYPE clover_force
#define _SITE_TYPE suNf
#define _REAL double
#define _COMPLEX hr_complex
#define _FIELD_DIM 6
#include "TMPL/linear_algebra_reduction.c.tmpl"
#include "TMPL/linear_algebra_base_operations.c.tmpl"
#include "TMPL/linear_algebra_base.c.tmpl"

#define _FIELD_TYPE staple_field
#define _SITE_TYPE suNg
#define _REAL double
#define _COMPLEX hr_complex
#define _FIELD_DIM 3
#include "TMPL/linear_algebra_reduction.c.tmpl"
#include "TMPL/linear_algebra_base_operations.c.tmpl"
#include "TMPL/linear_algebra_base.c.tmpl"

double spinor_prod_re_f(suNf_spinor *r, suNf_spinor *s) {
    double k = 0;
    _spinor_prod_re_f(k, *r, *s);
    return k;
}

double spinor_prod_re_f_flt(suNf_spinor_flt *r, suNf_spinor_flt *s) {
    double k = 0;
    _spinor_prod_re_f(k, *r, *s);
    return k;
}

double suNg_prod_re(suNg *u, suNg *v) {
    suNg tmp;
    double res = 0.0;
    _suNg_dagger_times_suNg(tmp, *u, *v);
    _suNg_trace_re(res, tmp);
    return res;
}

double suNf_prod_re(suNf *u, suNf *v) {
    suNf tmp;
    double res = 0.0;
    _suNf_dagger_times_suNf(tmp, *u, *v);
    _suNf_trace_re(res, tmp);
    return res;
}

double suNg_flt_prod_re(suNg_flt *u, suNg_flt *v) {
    suNg_flt tmp;
    double res = 0.0;
    _suNg_dagger_times_suNg(tmp, *u, *v);
    _suNg_trace_re(res, tmp);
    return res;
}

double suNf_flt_prod_re(suNf_flt *u, suNf_flt *v) {
    suNf_flt tmp;
    double res = 0.0;
    _suNf_dagger_times_suNf(tmp, *u, *v);
    _suNf_trace_re(res, tmp);
    return res;
}

double suNgc_prod_re(suNgc *u, suNgc *v) {
    suNgc tmp1, tmp2;
    hr_complex res = 0.0;
    _suNgc_dagger(tmp1, *u);
    _suNgc_times_suNgc(tmp2, tmp1, *v);
    _suNgc_trace(res, tmp2);
    return creal(res);
}

double suNfc_prod_re(suNfc *u, suNfc *v) {
    suNfc tmp1, tmp2;
    hr_complex res = 0.0;
    _suNfc_dagger(tmp1, *u);
    _suNfc_times_suNfc(tmp2, tmp1, *v);
    _suNfc_trace(res, tmp2);
    return creal(res);
}

double suNf_vector_prod_re(suNf_vector *r, suNf_vector *s) {
    double prod;
    _vector_prod_re_f(prod, *r, *s);
    return prod;
}

double suNg_vector_prod_re(suNg_vector *r, suNg_vector *s) {
    double prod;
    _vector_prod_re_g(prod, *r, *s);
    return prod;
}

double suNg_algebra_vector_prod_re(suNg_algebra_vector *r, suNg_algebra_vector *s) {
    double prod;
    _algebra_vector_prod_g(prod, *s, *r);
    return prod;
}

double spinor_prod_im_f(suNf_spinor *r, suNf_spinor *s) {
    double k = 0;
    _spinor_prod_im_f(k, *r, *s);
    return k;
}

double spinor_prod_im_f_flt(suNf_spinor_flt *r, suNf_spinor_flt *s) {
    double k = 0;
    _spinor_prod_im_f(k, *r, *s);
    return k;
}

double suNg_prod_im(suNg *u, suNg *v) {
    suNg tmp;
    double res = 0.0;
    _suNg_dagger_times_suNg(tmp, *u, *v);
    _suNg_trace_im(res, tmp);
    return res;
}

double suNf_prod_im(suNf *u, suNf *v) {
#ifdef REPR_IS_REAL
    return 0.0;
#else
    suNf tmp;
    double res = 0.0;
    _suNf_dagger_times_suNf(tmp, *u, *v);
    _suNf_trace_im(res, tmp);
    return res;
#endif
}

double suNg_flt_prod_im(suNg_flt *u, suNg_flt *v) {
    suNg_flt tmp;
    double res = 0.0;
    _suNg_dagger_times_suNg(tmp, *u, *v);
    _suNg_trace_im(res, tmp);
    return res;
}

double suNf_flt_prod_im(suNf_flt *u, suNf_flt *v) {
#ifdef REPR_IS_REAL
    return 0.0;
#else
    suNf_flt tmp;
    double res = 0.0;
    _suNf_dagger_times_suNf(tmp, *u, *v);
    _suNf_trace_im(res, tmp);
    return res;
#endif
}

double suNgc_prod_im(suNgc *u, suNgc *v) {
    suNgc tmp1;
    suNgc tmp2;
    _suNgc_unit(tmp2);
    hr_complex res = 0.0;
    _suNgc_dagger(tmp1, *u);
    _suNgc_multiply(tmp2, tmp1, *v);
    _suNgc_trace(res, tmp2);
    return cimag(res);
}

double suNfc_prod_im(suNfc *u, suNfc *v) {
    suNfc tmp1;
    suNfc tmp2;
    _suNfc_unit(tmp2);
    hr_complex res = 0.0;
    _suNfc_dagger(tmp1, *u);
    _suNfc_multiply(tmp2, tmp1, *v);
    _suNfc_trace(res, tmp2);
    return cimag(res);
}

double suNf_vector_prod_im(suNf_vector *r, suNf_vector *s) {
    double prod;
    _vector_prod_im_g(prod, *r, *s);
    return prod;
}

double suNg_vector_prod_im(suNg_vector *r, suNg_vector *s) {
    double prod;
    _vector_prod_im_g(prod, *r, *s);
    return prod;
}

double suNg_algebra_vector_prod_im(suNg_algebra_vector *r, suNg_algebra_vector *s) {
    return 0; // All real, please don't call this.
}

hr_complex spinor_prod_f(suNf_spinor *r, suNf_spinor *s) {
    hr_complex z = 0;
    _spinor_prod_f(z, *r, *s);
    return z;
}

hr_complex spinor_prod_f_flt(suNf_spinor_flt *r, suNf_spinor_flt *s) {
    hr_complex_flt z = 0;
    _spinor_prod_f(z, *r, *s);
    //hr_complex res = (double)creal(z) + I * (double)cimag(z);
    hr_complex res = spinor_prod_re_f_flt(r, s) + I * spinor_prod_im_f_flt(r, s);
    return res;
}

hr_complex suNg_prod(suNg *u, suNg *v) {
    suNg tmp;
    hr_complex res;
    _suNg_dagger_times_suNg(tmp, *u, *v);
    _suNg_trace(res, tmp);
    return res;
}

hr_complex suNf_prod(suNf *u, suNf *v) {
    suNf tmp;
    hr_complex res;
    _suNf_dagger_times_suNf(tmp, *u, *v);
    _suNf_trace(res, tmp);
    return res;
}

hr_complex suNg_flt_prod(suNg_flt *u, suNg_flt *v) {
    suNg_flt tmp;
    hr_complex res;
    _suNg_dagger_times_suNg(tmp, *u, *v);
    _suNg_trace_re(res, tmp);
    return res;
}

hr_complex suNf_flt_prod(suNf_flt *u, suNf_flt *v) {
    suNf_flt tmp;
    hr_complex res;
    _suNf_dagger_times_suNf(tmp, *u, *v);
    _suNf_trace_re(res, tmp);
    return res;
}

hr_complex suNgc_prod(suNgc *u, suNgc *v) {
    suNgc tmp1, tmp2;
    hr_complex res = 0.0;
    _suNg_dagger(tmp1, *u);
    _suNg_multiply(tmp2, tmp1, *v);
    _suNgc_trace(res, tmp2);
    return res;
}

hr_complex suNfc_prod(suNfc *u, suNfc *v) {
    suNfc tmp1, tmp2;
    hr_complex res = 0.0;
    _suNfc_dagger(tmp1, *u);
    _suNfc_multiply(tmp2, tmp1, *v);
    _suNfc_trace(res, tmp2);
    return res;
}

hr_complex suNf_vector_prod(suNf_vector *r, suNf_vector *s) {
    hr_complex prod;
    _vector_prod_f(prod, *r, *s);
    return prod;
}

hr_complex suNg_vector_prod(suNg_vector *r, suNg_vector *s) {
    hr_complex prod;
    _vector_prod_g(prod, *r, *s);
    return prod;
}

hr_complex suNg_algebra_vector_prod(suNg_algebra_vector *r, suNg_algebra_vector *s) {
    hr_complex prod;
    _algebra_vector_prod_g(prod, *s, *r);
    return prod;
}

double spinor_sqnorm_f(suNf_spinor *r) {
    double z = 0;
    _spinor_prod_re_f(z, *r, *r);
    return z;
}

double spinor_sqnorm_f_flt(suNf_spinor_flt *r) {
    double z = 0;
    _spinor_prod_re_f(z, *r, *r);
    return z;
}

double suNg_sqnorm(suNg *u) {
    double norm;
    _suNg_sqnorm(norm, *u);
    return norm;
}

double suNf_sqnorm(suNf *u) {
    double norm;
    _suNf_sqnorm(norm, *u);
    return norm;
}

double suNg_flt_sqnorm(suNg_flt *u) {
    double norm;
    _suNg_sqnorm(norm, *u);
    return norm;
}

double suNf_flt_sqnorm(suNf_flt *u) {
    double norm;
    _suNf_sqnorm(norm, *u);
    return norm;
}

double suNgc_sqnorm(suNgc *u) {
    double norm;
    _suNgc_sqnorm(norm, *u);
    return norm;
}

double suNfc_sqnorm(suNfc *u) {
    double norm;
    _suNfc_sqnorm(norm, *u);
    return norm;
}

double suNf_vector_sqnorm(suNf_vector *r) {
    double prod;
    _vector_prod_re_f(prod, *r, *r);
    return prod;
}

double suNg_vector_sqnorm(suNg_vector *r) {
    double prod;
    _vector_prod_re_g(prod, *r, *r);
    return prod;
}

double suNg_algebra_vector_sqnorm(suNg_algebra_vector *r) {
    double sqnorm;
    _algebra_vector_sqnorm_g(sqnorm, *r);
    return sqnorm;
}

double spinor_max_f(suNf_spinor *r) {
    double *a = (double *)r;
    double max = 0.;
    for (int i = 0; i < sizeof(suNf_spinor) / sizeof(*a); i++) {
        double v = fabs(a[i]);
        if (max < v) { max = v; }
    }
    return max;
}

double spinor_max_f_flt(suNf_spinor_flt *r) {
    double *a = (double *)r;
    double max = 0.;
    for (int i = 0; i < sizeof(suNf_spinor_flt) / sizeof(*a); i++) {
        double v = fabs(a[i]);
        if (max < v) { max = v; }
    }
    return max;
}

double suNg_max(suNg *r) {
    double *a = (double *)r;
    double max = 0.;
    for (int i = 0; i < sizeof(suNg) / sizeof(*a); i++) {
        double v = fabs(a[i]);
        if (max < v) { max = v; }
    }
    return max;
}

double suNf_max(suNf *r) {
    double *a = (double *)r;
    double max = 0.;
    for (int i = 0; i < sizeof(suNf) / sizeof(*a); i++) {
        double v = fabs(a[i]);
        if (max < v) { max = v; }
    }
    return max;
}

double suNg_flt_max(suNg_flt *r) {
    double *a = (double *)r;
    double max = 0.;
    for (int i = 0; i < sizeof(suNg_flt) / sizeof(*a); i++) {
        double v = fabs(a[i]);
        if (max < v) { max = v; }
    }
    return max;
}

double suNf_flt_max(suNf_flt *r) {
    double *a = (double *)r;
    double max = 0.;
    for (int i = 0; i < sizeof(suNf_flt) / sizeof(*a); i++) {
        double v = fabs(a[i]);
        if (max < v) { max = v; }
    }
    return max;
}

double suNgc_max(suNgc *r) {
    double *a = (double *)r;
    double max = 0.;
    for (int i = 0; i < sizeof(suNgc) / sizeof(*a); i++) {
        double v = fabs(a[i]);
        if (max < v) { max = v; }
    }
    return max;
}

double suNfc_max(suNfc *r) {
    double *a = (double *)r;
    double max = 0.;
    for (int i = 0; i < sizeof(suNfc) / sizeof(*a); i++) {
        double v = fabs(a[i]);
        if (max < v) { max = v; }
    }
    return max;
}

double suNf_vector_max(suNf_vector *r) {
    double *a = (double *)r;
    double max = 0.;
    for (int i = 0; i < sizeof(suNf_vector) / sizeof(*a); i++) {
        double v = fabs(a[i]);
        if (max < v) { max = v; }
    }
    return max;
}

double suNg_vector_max(suNg_vector *r) {
    double *a = (double *)r;
    double max = 0.;
    for (int i = 0; i < sizeof(suNg_vector) / sizeof(*a); i++) {
        double v = fabs(a[i]);
        if (max < v) { max = v; }
    }
    return max;
}

double suNg_algebra_vector_max(suNg_algebra_vector *r) {
    double *a = (double *)r;
    double max = 0.;
    for (int i = 0; i < sizeof(suNg_algebra_vector) / sizeof(*a); i++) {
        double v = fabs(a[i]);
        if (max < v) { max = v; }
    }
    return max;
}

double spinor_g5_prod_re_f(suNf_spinor *r, suNf_spinor *s) {
    double k = 0;
    _spinor_g5_prod_re_f(k, *r, *s);
    return k;
}

double spinor_g5_prod_re_f_flt(suNf_spinor_flt *r, suNf_spinor_flt *s) {
    double z = 0;
    _spinor_g5_prod_re_f(z, *r, *s);
    return z;
}

double spinor_g5_prod_im_f(suNf_spinor *r, suNf_spinor *s) {
    double k = 0;
    _spinor_g5_prod_im_f(k, *r, *s);
    return k;
}

double spinor_g5_prod_im_f_flt(suNf_spinor_flt *r, suNf_spinor_flt *s) {
    double z = 0;
    _spinor_g5_prod_im_f(z, *r, *s);
    return z;
}
