#ifndef LINEAR_ALGEBRA_GENERIC_GPU_HPP
#define LINEAR_ALGEBRA_GENERIC_GPU_HPP

#ifdef __cplusplus

#include "inverters.h"
#include "libhr_core.h"

#define _FUNC_GENERIC(_type, _name, _args) _type _name _args
#define _DECLARE_LINA_HEADER(a, b, c) a b c;

// double precision
#define _FIELD_TYPE spinor_field
#define _REAL double
#define _COMPLEX hr_complex
#include "TMPL/linear_algebra_reduction.h.tmpl"
#include "TMPL/linear_algebra_base_operations.h.tmpl"
#include "TMPL/linear_algebra_lc.h.tmpl"
#include "TMPL/linear_algebra_gamma.h.tmpl"
#include "TMPL/linear_algebra_base.h.tmpl" // This has to come last

// single precision
#define _FIELD_TYPE spinor_field_flt
#define _REAL float
#define _COMPLEX hr_complex_flt
#include "TMPL/linear_algebra_reduction.h.tmpl"
#include "TMPL/linear_algebra_base_operations.h.tmpl"
#include "TMPL/linear_algebra_lc.h.tmpl"
#include "TMPL/linear_algebra_gamma.h.tmpl"
#include "TMPL/linear_algebra_base.h.tmpl"

// double precision
#define _FIELD_TYPE scalar_field
#define _REAL double
#define _COMPLEX hr_complex
#include "TMPL/linear_algebra_reduction.h.tmpl"
#include "TMPL/linear_algebra_base_operations.h.tmpl"
#include "TMPL/linear_algebra_base.h.tmpl"

#define _FIELD_TYPE suNg_field
#define _REAL double
#define _COMPLEX hr_complex
#include "TMPL/linear_algebra_reduction.h.tmpl"
#include "TMPL/linear_algebra_base_operations.h.tmpl"
#include "TMPL/linear_algebra_base.h.tmpl"

#define _FIELD_TYPE suNf_field
#define _REAL double
#define _COMPLEX hr_complex
#include "TMPL/linear_algebra_reduction.h.tmpl"
#include "TMPL/linear_algebra_base_operations.h.tmpl"
#include "TMPL/linear_algebra_base.h.tmpl"

#define _FIELD_TYPE suNfc_field
#define _REAL double
#define _COMPLEX hr_complex
#include "TMPL/linear_algebra_reduction.h.tmpl"
#include "TMPL/linear_algebra_base_operations.h.tmpl"
#include "TMPL/linear_algebra_base.h.tmpl"

#define _FIELD_TYPE suNg_field_flt
#define _REAL float
#define _COMPLEX hr_complex_flt
#include "TMPL/linear_algebra_reduction.h.tmpl"
#include "TMPL/linear_algebra_base_operations.h.tmpl"
#include "TMPL/linear_algebra_base.h.tmpl"

#define _FIELD_TYPE suNf_field_flt
#define _REAL float
#define _COMPLEX hr_complex_flt
#include "TMPL/linear_algebra_reduction.h.tmpl"
#include "TMPL/linear_algebra_base_operations.h.tmpl"
#include "TMPL/linear_algebra_base.h.tmpl"

#define _FIELD_TYPE suNg_scalar_field
#define _REAL double
#define _COMPLEX hr_complex
#include "TMPL/linear_algebra_reduction.h.tmpl"
#include "TMPL/linear_algebra_base_operations.h.tmpl"
#include "TMPL/linear_algebra_base.h.tmpl"

#define _FIELD_TYPE suNg_av_field
#define _REAL double
#define _COMPLEX hr_complex
#include "TMPL/linear_algebra_reduction.h.tmpl"
#include "TMPL/linear_algebra_base_operations.h.tmpl"
#include "TMPL/linear_algebra_base.h.tmpl"

#define _FIELD_TYPE gtransf
#define _REAL double
#define _COMPLEX hr_complex
#include "TMPL/linear_algebra_reduction.h.tmpl"
#include "TMPL/linear_algebra_base_operations.h.tmpl"
#include "TMPL/linear_algebra_base.h.tmpl"

#define _FIELD_TYPE ldl_field
#define _REAL double
#define _COMPLEX hr_complex
#include "TMPL/linear_algebra_base.h.tmpl"

#define _FIELD_TYPE clover_term
#define _REAL double
#define _COMPLEX hr_complex
#include "TMPL/linear_algebra_reduction.h.tmpl"
#include "TMPL/linear_algebra_base_operations.h.tmpl"
#include "TMPL/linear_algebra_base.h.tmpl"

#define _FIELD_TYPE clover_force
#define _REAL double
#define _COMPLEX hr_complex
#include "TMPL/linear_algebra_reduction.h.tmpl"
#include "TMPL/linear_algebra_base_operations.h.tmpl"
#include "TMPL/linear_algebra_base.h.tmpl"

#define _FIELD_TYPE staple_field
#define _REAL double
#define _COMPLEX hr_complex
#include "TMPL/linear_algebra_reduction.h.tmpl"
#include "TMPL/linear_algebra_base_operations.h.tmpl"
#include "TMPL/linear_algebra_base.h.tmpl"

// Cannot use error here because its undefined in device code
//#define template_error error(1, 1, __func__, "Complex multiplication of real-valued field \n");
#define template_error

visible __forceinline__ void mul_add_assign(suNf_spinor *s1, double rho, suNf_spinor *s2) {
    _spinor_mul_add_assign_f(*s1, rho, *s2);
}

visible __forceinline__ void mul_add_assign(suNf_spinor_flt *s1, double rho, suNf_spinor_flt *s2) {
    _spinor_mul_add_assign_f(*s1, rho, *s2);
}

visible __forceinline__ void mul_add_assign(suNf *s1, double rho, suNf *s2) {
    _suNf_mul_add(*s1, 1.0, *s1, rho, *s2);
}

#ifdef REPR_IS_REAL
visible __forceinline__ void mul_add_assign(suNfc *s1, double rho, suNfc *s2) {
    _suNfc_mul_add(*s1, 1.0, *s1, rho, *s2);
}
#endif

visible __forceinline__ void mul_add_assign(suNg *s1, double rho, suNg *s2) {
    _suNg_mul_add(*s1, 1.0, *s1, rho, *s2);
}

visible __forceinline__ void mul_add_assign(suNf_flt *s1, double rho, suNf_flt *s2) {
    _suNf_mul_add(*s1, 1.0, *s1, rho, *s2);
}

visible __forceinline__ void mul_add_assign(suNg_flt *s1, double rho, suNg_flt *s2) {
    _suNg_mul_add(*s1, 1.0, *s1, rho, *s2);
}

visible __forceinline__ void mul_add_assign(suNf_vector *s1, double rho, suNf_vector *s2) {
    _vector_mul_add_assign_f(*s1, rho, *s2);
}

visible __forceinline__ void mul_add_assign(suNg_vector *s1, double rho, suNg_vector *s2) {
    _vector_mul_add_assign_g(*s1, rho, *s2);
}

visible __forceinline__ void mul_add_assign(suNg_algebra_vector *s1, double rho, suNg_algebra_vector *s2) {
    _algebra_vector_mul_add_assign_g(*s1, rho, *s2);
}

visible __forceinline__ void mul_add_assign(double *s1, double rho, double *s2) {
    (*s1) += rho * (*s2);
}

visible __forceinline__ void mul_add_assign(float *s1, float rho, float *s2) {
    (*s1) += rho * (*s2);
}

visible __forceinline__ void mulc_add_assign(suNf_spinor *s1, hr_complex rho, suNf_spinor *s2) {
    _spinor_mulc_add_assign_f(*s1, rho, *s2);
}

visible __forceinline__ void mulc_add_assign(suNf_spinor_flt *s1, hr_complex_flt rho, suNf_spinor_flt *s2) {
    _spinor_mulc_add_assign_f(*s1, rho, *s2);
}

visible __forceinline__ void mulc_add_assign(suNf *s1, hr_complex rho, suNf *s2) {
#ifdef REPR_IS_REAL
    template_error;
#else
    suNf tmp;
    _suNf_mulc(tmp, rho, *s2);
    _suNf_add_assign(*s1, tmp);
#endif
}

#ifdef REPR_IS_REAL
visible __forceinline__ void mulc_add_assign(suNfc *s1, hr_complex rho, suNfc *s2) {
    suNfc tmp;
    _suNfc_mulc(tmp, rho, *s2);
    _suNfc_add_assign(*s1, tmp);
}
#endif

visible __forceinline__ void mulc_add_assign(suNg *s1, hr_complex rho, suNg *s2) {
    suNg tmp;
    _suNg_mulc(tmp, rho, *s2);
    _suNg_add_assign(*s1, tmp);
}

visible __forceinline__ void mulc_add_assign(suNf_flt *s1, hr_complex_flt rho, suNf_flt *s2) {
#ifdef REPR_IS_REAL
    template_error;
#else
    suNf_flt tmp;
    _suNf_mulc(tmp, rho, *s2);
    _suNf_add_assign(*s1, tmp);
#endif
}

visible __forceinline__ void mulc_add_assign(suNg_flt *s1, hr_complex_flt rho, suNg_flt *s2) {
    suNg_flt tmp;
    _suNg_mulc(tmp, rho, *s2);
    _suNg_add_assign(*s1, tmp);
}

visible __forceinline__ void mulc_add_assign(suNf_vector *s1, hr_complex rho, suNf_vector *s2) {
    _vector_mulc_add_assign_f(*s1, rho, *s2);
}

visible __forceinline__ void mulc_add_assign(suNg_vector *s1, hr_complex rho, suNg_vector *s2) {
    _vector_mulc_add_assign_f(*s1, rho, *s2);
}

visible __forceinline__ void mulc_add_assign(suNg_algebra_vector *s1, hr_complex rho, suNg_algebra_vector *s2) {
    _algebra_vector_mul_g(*s1, creal(rho), *s2);
}

visible __forceinline__ void mulc_add_assign(double *s1, hr_complex rho, double *s2) {
    // TODO: this needs to throw an error instead.
    (*s1) += creal(rho) * (*s2);
}

visible __forceinline__ void mulc_add_assign(float *s1, hr_complex_flt rho, float *s2) {
    // ERROR
    (*s1) += creal(rho) * (*s2);
}

visible __forceinline__ void mul(suNf_spinor *s1, double rho, suNf_spinor *s2) {
    _spinor_mul_f(*s1, rho, *s2);
}

visible __forceinline__ void mul(suNf_spinor_flt *s1, double rho, suNf_spinor_flt *s2) {
    _spinor_mul_f(*s1, rho, *s2);
}

visible __forceinline__ void mul(suNf *s1, double rho, suNf *s2) {
    _suNf_mul(*s1, rho, *s2);
}

#ifdef REPR_IS_REAL
visible __forceinline__ void mul(suNfc *s1, double rho, suNfc *s2) {
    _suNfc_mul(*s1, rho, *s2);
}
#endif

visible __forceinline__ void mul(suNg *s1, double rho, suNg *s2) {
    _suNg_mul(*s1, rho, *s2);
}

visible __forceinline__ void mul(suNf_flt *s1, double rho, suNf_flt *s2) {
    _suNf_mul(*s1, rho, *s2);
}

visible __forceinline__ void mul(suNg_flt *s1, double rho, suNg_flt *s2) {
    _suNg_mul(*s1, rho, *s2);
}

visible __forceinline__ void mul(suNf_vector *s1, double rho, suNf_vector *s2) {
    _vector_mul_f(*s1, rho, *s2);
}

visible __forceinline__ void mul(suNg_vector *s1, double rho, suNg_vector *s2) {
    _vector_mul_g(*s1, rho, *s2);
}

visible __forceinline__ void mul(suNg_algebra_vector *s1, double rho, suNg_algebra_vector *s2) {
    _algebra_vector_mul_g(*s1, rho, *s2);
}

visible __forceinline__ void mul(double *s1, double rho, double *s2) {
    (*s1) = rho * (*s2);
}

visible __forceinline__ void mul(float *s1, float rho, float *s2) {
    (*s1) = rho * (*s2);
}

visible __forceinline__ void mulc(suNf_spinor *s1, hr_complex rho, suNf_spinor *s2) {
    _spinor_mulc_f(*s1, rho, *s2);
}

visible __forceinline__ void mulc(suNf_spinor_flt *s1, hr_complex_flt rho, suNf_spinor_flt *s2) {
    _spinor_mulc_f(*s1, rho, *s2);
}

visible __forceinline__ void mulc(suNf *s1, hr_complex rho, suNf *s2) {
#ifdef REPR_IS_REAL
    template_error;
#else
    _suNf_mulc(*s1, rho, *s2);
#endif
}

#ifdef REPR_IS_REAL

visible __forceinline__ void mulc(suNfc *s1, hr_complex rho, suNfc *s2) {
    _suNfc_mulc(*s1, rho, *s2);
}
#endif

visible __forceinline__ void mulc(suNg *s1, hr_complex rho, suNg *s2) {
    _suNg_mulc(*s1, rho, *s2);
}

visible __forceinline__ void mulc(suNf_flt *s1, hr_complex_flt rho, suNf_flt *s2) {
#ifdef REPR_IS_REAL
    template_error;
#else
    _suNf_mulc(*s1, rho, *s2);
#endif
}

visible __forceinline__ void mulc(suNg_flt *s1, hr_complex_flt rho, suNg_flt *s2) {
    _suNg_mulc(*s1, rho, *s2);
}

visible __forceinline__ void mulc(suNf_vector *s1, hr_complex rho, suNf_vector *s2) {
    _vector_mulc_f(*s1, rho, *s2);
}

visible __forceinline__ void mulc(suNg_vector *s1, hr_complex rho, suNg_vector *s2) {
    _vector_mulc_g(*s1, rho, *s2);
}

visible __forceinline__ void mulc(suNg_algebra_vector *s1, hr_complex rho, suNg_algebra_vector *s2) {
    _algebra_vector_mul_g(*s1, creal(rho), *s2);
}

visible __forceinline__ void mulc(double *s1, hr_complex rho, double *s2) {
    // Error!
}

visible __forceinline__ void mulc(float *s1, hr_complex rho, float *s2) {
    // Error
}

visible __forceinline__ void add(suNf_spinor *r, suNf_spinor *s1, suNf_spinor *s2) {
    _spinor_add_f(*r, *s1, *s2);
}

visible __forceinline__ void add(suNf_spinor_flt *r, suNf_spinor_flt *s1, suNf_spinor_flt *s2) {
    _spinor_add_f(*r, *s1, *s2);
}

visible __forceinline__ void add(suNf *r, suNf *s1, suNf *s2) {
    _suNf_mul_add(*r, 1.0, *s1, 1.0, *s2);
}

#ifdef REPR_IS_REAL
visible __forceinline__ void add(suNfc *r, suNfc *s1, suNfc *s2) {
    _suNf_mul_add(*r, 1.0, *s1, 1.0, *s2);
}
#endif

visible __forceinline__ void add(suNg *r, suNg *s1, suNg *s2) {
    _suNg_mul_add(*r, 1.0, *s1, 1.0, *s2);
}

visible __forceinline__ void add(suNf_flt *r, suNf_flt *s1, suNf_flt *s2) {
    _suNf_mul_add(*r, 1.0, *s1, 1.0, *s2);
}

visible __forceinline__ void add(suNg_flt *r, suNg_flt *s1, suNg_flt *s2) {
    _suNg_mul_add(*r, 1.0, *s1, 1.0, *s2);
}

visible __forceinline__ void add(suNf_vector *r, suNf_vector *s1, suNf_vector *s2) {
    _vector_add_f(*r, *s1, *s2);
}

visible __forceinline__ void add(suNg_vector *r, suNg_vector *s1, suNg_vector *s2) {
    _vector_add_g(*r, *s1, *s2);
}

visible __forceinline__ void add(suNg_algebra_vector *r, suNg_algebra_vector *s1, suNg_algebra_vector *s2) {
    _algebra_vector_zero_g(*r);
    _algebra_vector_add_assign_g(*r, *s1);
    _algebra_vector_add_assign_g(*r, *s2);
}

visible __forceinline__ void add(double *r, double *s1, double *s2) {
    (*r) = (*s1) + (*s2);
}

visible __forceinline__ void add(float *r, float *s1, float *s2) {
    (*r) = (*s1) + (*s2);
}

visible __forceinline__ void sub(suNf_spinor *r, suNf_spinor *s1, suNf_spinor *s2) {
    _spinor_sub_f(*r, *s1, *s2);
}

visible __forceinline__ void sub(suNf_spinor_flt *r, suNf_spinor_flt *s1, suNf_spinor_flt *s2) {
    _spinor_sub_f(*r, *s1, *s2);
}

visible __forceinline__ void sub(suNf *r, suNf *s1, suNf *s2) {
    _suNf_mul_add(*r, 1.0, *s1, -1.0, *s2);
}

#ifdef REPR_IS_REAL
visible __forceinline__ void sub(suNfc *r, suNfc *s1, suNfc *s2) {
    _suNfc_mul_add(*r, 1.0, *s1, -1.0, *s2);
}
#endif

visible __forceinline__ void sub(suNg *r, suNg *s1, suNg *s2) {
    _suNg_mul_add(*r, 1.0, *s1, -1.0, *s2);
}

visible __forceinline__ void sub(suNf_flt *r, suNf_flt *s1, suNf_flt *s2) {
    _suNf_mul_add(*r, 1.0, *s1, -1.0, *s2);
}

visible __forceinline__ void sub(suNg_flt *r, suNg_flt *s1, suNg_flt *s2) {
    _suNg_mul_add(*r, 1.0, *s1, -1.0, *s2);
}

visible __forceinline__ void sub(suNf_vector *r, suNf_vector *s1, suNf_vector *s2) {
    _vector_sub_f(*r, *s1, *s2);
}

visible __forceinline__ void sub(suNg_vector *r, suNg_vector *s1, suNg_vector *s2) {
    _vector_sub_g(*r, *s1, *s2);
}

visible __forceinline__ void sub(suNg_algebra_vector *r, suNg_algebra_vector *s1, suNg_algebra_vector *s2) {
    _algebra_vector_zero_g(*r);
    _algebra_vector_add_assign_g(*r, *s1);
    _algebra_vector_sub_assign_g(*r, *s2);
}

visible __forceinline__ void sub(double *r, double *s1, double *s2) {
    (*r) = (*s1) - (*s2);
}

visible __forceinline__ void sub(float *r, float *s1, float *s2) {
    (*r) = (*s1) - (*s2);
}

visible __forceinline__ void sub_assign(suNf_spinor *s1, suNf_spinor *s2) {
    _spinor_sub_assign_f(*s1, *s2);
}

visible __forceinline__ void sub_assign(suNf_spinor_flt *s1, suNf_spinor_flt *s2) {
    _spinor_sub_assign_f(*s1, *s2);
}

visible __forceinline__ void sub_assign(suNf *s1, suNf *s2) {
    _suNf_sub_assign(*s1, *s2);
}

#ifdef REPR_IS_REAL
visible __forceinline__ void sub_assign(suNfc *s1, suNfc *s2) {
    _suNfc_sub_assign(*s1, *s2);
}
#endif

visible __forceinline__ void sub_assign(suNg *s1, suNg *s2) {
    _suNg_sub_assign(*s1, *s2);
}

visible __forceinline__ void sub_assign(suNf_flt *s1, suNf_flt *s2) {
    _suNf_sub_assign(*s1, *s2);
}

visible __forceinline__ void sub_assign(suNg_flt *s1, suNg_flt *s2) {
    _suNg_sub_assign(*s1, *s2);
}

visible __forceinline__ void sub_assign(suNf_vector *s1, suNf_vector *s2) {
    _vector_sub_assign_f(*s1, *s2);
}

visible __forceinline__ void sub_assign(suNg_vector *s1, suNg_vector *s2) {
    _vector_sub_assign_g(*s1, *s2);
}

visible __forceinline__ void sub_assign(suNg_algebra_vector *s1, suNg_algebra_vector *s2) {
    _algebra_vector_sub_assign_g(*s1, *s2);
}

visible __forceinline__ void sub_assign(double *s1, double *s2) {
    (*s1) -= (*s2);
}

visible __forceinline__ void sub_assign(float *s1, float *s2) {
    (*s1) -= (*s2);
}

visible __forceinline__ void minus(suNf_spinor *s1, suNf_spinor *s2) {
    _spinor_minus_f(*s1, *s2);
}

visible __forceinline__ void minus(suNf_spinor_flt *s1, suNf_spinor_flt *s2) {
    _spinor_minus_f(*s1, *s2);
}

visible __forceinline__ void minus(suNf *s1, suNf *s2) {
    _suNf_minus(*s1, *s2);
}

#ifdef REPR_IS_REAL
visible __forceinline__ void minus(suNfc *s1, suNfc *s2) {
    _suNfc_minus(*s1, *s2);
}
#endif

visible __forceinline__ void minus(suNg *s1, suNg *s2) {
    _suNg_minus(*s1, *s2);
}

visible __forceinline__ void minus(suNf_flt *s1, suNf_flt *s2) {
    _suNf_minus(*s1, *s2);
}

visible __forceinline__ void minus(suNg_flt *s1, suNg_flt *s2) {
    _suNg_minus(*s1, *s2);
}

visible __forceinline__ void minus(suNf_vector *s1, suNf_vector *s2) {
    _vector_minus_f(*s1, *s2);
}

visible __forceinline__ void minus(suNg_vector *s1, suNg_vector *s2) {
    _vector_minus_g(*s1, *s2);
}

visible __forceinline__ void minus(suNg_algebra_vector *s1, suNg_algebra_vector *s2) {
    _algebra_vector_mul_g(*s1, -1.0, *s2);
}

visible __forceinline__ void minus(double *s1, double *s2) {
    (*s1) = -(*s2);
}

visible __forceinline__ void minus(float *s1, float *s2) {
    (*s1) = -(*s2);
}

visible __forceinline__ void add_assign(suNf_spinor *s1, suNf_spinor *s2) {
    _spinor_add_assign_f(*s1, *s2);
}

visible __forceinline__ void add_assign(suNf_spinor_flt *s1, suNf_spinor_flt *s2) {
    _spinor_add_assign_f(*s1, *s2);
}

visible __forceinline__ void add_assign(suNf *s1, suNf *s2) {
    _suNf_add_assign(*s1, *s2);
}

#ifdef REPR_IS_REAL
visible __forceinline__ void add_assign(suNfc *s1, suNfc *s2) {
    _suNfc_add_assign(*s1, *s2);
}
#endif

visible __forceinline__ void add_assign(suNg *s1, suNg *s2) {
    _suNg_add_assign(*s1, *s2);
}

visible __forceinline__ void add_assign(suNf_flt *s1, suNf_flt *s2) {
    _suNf_add_assign(*s1, *s2);
}

visible __forceinline__ void add_assign(suNg_flt *s1, suNg_flt *s2) {
    _suNg_add_assign(*s1, *s2);
}

visible __forceinline__ void add_assign(suNf_vector *s1, suNf_vector *s2) {
    _vector_add_assign_f(*s1, *s2);
}

visible __forceinline__ void add_assign(suNg_vector *s1, suNg_vector *s2) {
    _vector_add_assign_g(*s1, *s2);
}

visible __forceinline__ void add_assign(suNg_algebra_vector *s1, suNg_algebra_vector *s2) {
    _algebra_vector_add_assign_g(*s1, *s2);
}

visible __forceinline__ void add_assign(double *s1, double *s2) {
    (*s1) += (*s2);
}

visible __forceinline__ void add_assign(float *s1, float *s2) {
    (*s1) += (*s2);
}

visible __forceinline__ double prod_re(suNf_spinor *s1, suNf_spinor *s2) {
    double k = 0;
    _spinor_prod_re_f(k, *s1, *s2);
    return k;
}

visible __forceinline__ double prod_re(suNf_spinor_flt *s1, suNf_spinor_flt *s2) {
    double k = 0;
    _spinor_prod_re_f(k, *s1, *s2);
    return k;
}

visible __forceinline__ double prod_re(suNf *s1, suNf *s2) {
    suNf tmp;
    double res = 0.0;
    _suNf_dagger_times_suNf(tmp, *s1, *s2);
    _suNf_trace_re(res, tmp);
    return res;
}

#ifdef REPR_IS_REAL
visible __forceinline__ double prod_re(suNfc *s1, suNfc *s2) {
    suNfc tmp;
    hr_complex res = 0.0;
    _suNfc_dagger_times_suNfc(tmp, *s1, *s2);
    _suNfc_trace(res, tmp);
    return creal(res);
}
#endif

visible __forceinline__ double prod_re(suNg *s1, suNg *s2) {
    suNg tmp;
    double res = 0.0;
    _suNg_dagger_times_suNg(tmp, *s1, *s2);
    _suNg_trace_re(res, tmp);
    return res;
}

visible __forceinline__ double prod_re(suNf_flt *s1, suNf_flt *s2) {
    suNf_flt tmp;
    double res = 0.0;
    _suNf_dagger_times_suNf(tmp, *s1, *s2);
    _suNf_trace_re(res, tmp);
    return res;
}

visible __forceinline__ double prod_re(suNg_flt *s1, suNg_flt *s2) {
    suNg_flt tmp;
    double res = 0.0;
    _suNg_dagger_times_suNg(tmp, *s1, *s2);
    _suNg_trace_re(res, tmp);
    return res;
}

visible __forceinline__ double prod_re(suNf_vector *s1, suNf_vector *s2) {
    double prod;
    _vector_prod_re_f(prod, *s1, *s2);
    return prod;
}

visible __forceinline__ double prod_re(suNg_vector *s1, suNg_vector *s2) {
    double prod;
    _vector_prod_re_g(prod, *s1, *s2);
    return prod;
}

visible __forceinline__ double prod_re(suNg_algebra_vector *s1, suNg_algebra_vector *s2) {
    double prod;
    _algebra_vector_prod_g(prod, *s1, *s2);
    return prod;
}

visible __forceinline__ double prod_re(double *s1, double *s2) {
    return (*s1) * (*s2);
}

visible __forceinline__ double prod_re(float *s1, float *s2) {
    return (*s1) * (*s2);
}

visible __forceinline__ double prod_im(suNf_spinor *s1, suNf_spinor *s2) {
    double k = 0;
    _spinor_prod_im_f(k, *s1, *s2);
    return k;
}

visible __forceinline__ double prod_im(suNf_spinor_flt *s1, suNf_spinor_flt *s2) {
    double k = 0;
    _spinor_prod_im_f(k, *s1, *s2);
    return k;
}

visible __forceinline__ double prod_im(suNf *s1, suNf *s2) {
#ifdef REPR_IS_REAL
    return 0.0;
#else
    suNf tmp;
    double res = 0.0;
    _suNf_dagger_times_suNf(tmp, *s1, *s2);
    _suNf_trace_im(res, tmp);
    return res;
#endif
}

#ifdef REPR_IS_REAL
visible __forceinline__ double prod_im(suNfc *s1, suNfc *s2) {
    return 0.0;
    /*suNfc tmp;
    double res = 0.0;
    _suNfc_dagger_times_suNfc(tmp, *s1, *s2);
    _suNfc_trace_im(res, tmp);
    return res;*/
}
#endif

visible __forceinline__ double prod_im(suNg *s1, suNg *s2) {
    suNg tmp;
    double res = 0.0;
    _suNg_dagger_times_suNg(tmp, *s1, *s2);
    _suNg_trace_im(res, tmp);
    return res;
}

visible __forceinline__ double prod_im(suNf_flt *s1, suNf_flt *s2) {
#ifdef REPR_IS_REAL
    return 0.0;
#else
    suNf_flt tmp;
    double res = 0.0;
    _suNf_dagger_times_suNf(tmp, *s1, *s2);
    _suNf_trace_im(res, tmp);
    return res;
#endif
}

visible __forceinline__ double prod_im(suNg_flt *s1, suNg_flt *s2) {
    suNg_flt tmp;
    double res = 0.0;
    _suNg_dagger_times_suNg(tmp, *s1, *s2);
    _suNg_trace_im(res, tmp);
    return res;
}

visible __forceinline__ double prod_im(suNf_vector *s1, suNf_vector *s2) {
    double prod;
    _vector_prod_im_g(prod, *s1, *s2);
    return prod;
}

visible __forceinline__ double prod_im(suNg_vector *s1, suNg_vector *s2) {
    double prod;
    _vector_prod_im_g(prod, *s1, *s2);
    return prod;
}

visible __forceinline__ double prod_im(suNg_algebra_vector *s1, suNg_algebra_vector *s2) {
    return 0;
}

visible __forceinline__ double prod_im(double *s1, double *s2) {
    return 0.0;
}

visible __forceinline__ double prod_im(float *s1, float *s2) {
    return 0.0;
}

visible __forceinline__ hr_complex prod(suNf_spinor *s1, suNf_spinor *s2) {
    hr_complex z = 0;
    _spinor_prod_f(z, *s1, *s2);
    return z;
}

visible __forceinline__ hr_complex prod(suNf_spinor_flt *s1, suNf_spinor_flt *s2) {
    hr_complex_flt k = 0;
    _spinor_prod_re_f(k, *s1, *s2);
    hr_complex k1 = (double)creal(k) + I * (double)cimag(k);
    return k1;
}

visible __forceinline__ hr_complex prod(suNf *s1, suNf *s2) {
    suNf tmp;
    hr_complex res;
    _suNf_dagger_times_suNf(tmp, *s1, *s2);
    _suNf_trace(res, tmp);
    return res;
}

#ifdef REPR_IS_REAL
visible __forceinline__ hr_complex prod(suNfc *s1, suNfc *s2) {
    suNfc tmp;
    hr_complex res;
    _suNfc_dagger_times_suNfc(tmp, *s1, *s2);
    _suNfc_trace(res, tmp);
    return res;
}
#endif

visible __forceinline__ hr_complex prod(suNg *s1, suNg *s2) {
    suNg tmp;
    hr_complex res;
    _suNg_dagger_times_suNg(tmp, *s1, *s2);
    _suNg_trace(res, tmp);
    return res;
}

visible __forceinline__ hr_complex prod(suNf_flt *s1, suNf_flt *s2) {
    suNf_flt tmp;
    hr_complex res;
    _suNf_dagger_times_suNf(tmp, *s1, *s2);
    _suNf_trace_re(res, tmp);
    return res;
}

visible __forceinline__ hr_complex prod(suNg_flt *s1, suNg_flt *s2) {
    suNg_flt tmp;
    hr_complex res;
    _suNg_dagger_times_suNg(tmp, *s1, *s2);
    _suNg_trace_re(res, tmp);
    return res;
}

visible __forceinline__ hr_complex prod(suNf_vector *s1, suNf_vector *s2) {
    hr_complex prod;
    _vector_prod_f(prod, *s1, *s2);
    return prod;
}

visible __forceinline__ hr_complex prod(suNg_vector *s1, suNg_vector *s2) {
    hr_complex prod;
    _vector_prod_g(prod, *s1, *s2);
    return prod;
}

visible __forceinline__ hr_complex prod(suNg_algebra_vector *s1, suNg_algebra_vector *s2) {
    hr_complex prod;
    _algebra_vector_prod_g(prod, *s1, *s2);
    return prod;
}

visible __forceinline__ hr_complex prod(double *s1, double *s2) {
    hr_complex prod = (*s1) * (*s2);
    return prod;
}

visible __forceinline__ hr_complex prod(float *s1, float *s2) {
    hr_complex prod = ((double)((*s1) * (*s2)));
    return prod;
}

visible __forceinline__ double sqnorm(suNf_spinor *r) {
    double z = 0;
    _spinor_prod_re_f(z, *r, *r);
    return z;
}

visible __forceinline__ double sqnorm(suNf_spinor_flt *r) {
    double z = 0;
    _spinor_prod_re_f(z, *r, *r);
    return z;
}

visible __forceinline__ double sqnorm(suNf *r) {
    double norm;
    _suNf_sqnorm(norm, *r);
    return norm;
}

#ifdef REPR_IS_REAL
visible __forceinline__ double sqnorm(suNfc *r) {
    double norm;
    _suNfc_sqnorm(norm, *r);
    return norm;
}
#endif

visible __forceinline__ double sqnorm(suNg *r) {
    double norm;
    _suNg_sqnorm(norm, *r);
    return norm;
}

visible __forceinline__ double sqnorm(suNf_flt *r) {
    double norm;
    _suNf_sqnorm(norm, *r);
    return norm;
}

visible __forceinline__ double sqnorm(suNg_flt *r) {
    double norm;
    _suNg_sqnorm(norm, *r);
    return norm;
}

visible __forceinline__ double sqnorm(suNf_vector *r) {
    double prod;
    _vector_prod_re_f(prod, *r, *r);
    return prod;
}

visible __forceinline__ double sqnorm(suNg_vector *r) {
    double prod;
    _vector_prod_re_g(prod, *r, *r);
    return prod;
}

visible __forceinline__ double sqnorm(suNg_algebra_vector *r) {
    double sqnorm;
    _algebra_vector_sqnorm_g(sqnorm, *r);
    return sqnorm;
}

visible __forceinline__ double sqnorm(double *r) {
    return (*r) * (*r);
}

visible __forceinline__ double sqnorm(float *r) {
    return ((double)((*r) * (*r)));
}

visible __forceinline__ double max(suNf_spinor *r) {
    double *a = (double *)r;
    double max = 0.;
    for (int i = 0; i < sizeof(suNf_spinor) / sizeof(*a); i++) {
        double v = fabs(a[i]);
        if (max < v) { max = v; }
    }
    return max;
}

visible __forceinline__ double max(suNf_spinor_flt *r) {
    double *a = (double *)r;
    double max = 0.;
    for (int i = 0; i < sizeof(suNf_spinor_flt) / sizeof(*a); i++) {
        double v = fabs(a[i]);
        if (max < v) { max = v; }
    }
    return max;
}
visible __forceinline__ double max(suNf *r) {
    double *a = (double *)r;
    double max = 0.;
    for (int i = 0; i < sizeof(suNf) / sizeof(*a); i++) {
        double v = fabs(a[i]);
        if (max < v) { max = v; }
    }
    return max;
}
#ifdef REPR_IS_REAL
visible __forceinline__ double max(suNfc *r) {
    double *a = (double *)r;
    double max = 0.;
    for (int i = 0; i < sizeof(suNfc) / sizeof(*a); i++) {
        double v = fabs(a[i]);
        if (max < v) { max = v; }
    }
    return max;
}
#endif
visible __forceinline__ double max(suNg *r) {
    double *a = (double *)r;
    double max = 0.;
    for (int i = 0; i < sizeof(suNg) / sizeof(*a); i++) {
        double v = fabs(a[i]);
        if (max < v) { max = v; }
    }
    return max;
}
visible __forceinline__ double max(suNf_flt *r) {
    double *a = (double *)r;
    double max = 0.;
    for (int i = 0; i < sizeof(suNf_flt) / sizeof(*a); i++) {
        double v = fabs(a[i]);
        if (max < v) { max = v; }
    }
    return max;
}
visible __forceinline__ double max(suNg_flt *r) {
    double *a = (double *)r;
    double max = 0.;
    for (int i = 0; i < sizeof(suNg_flt) / sizeof(*a); i++) {
        double v = fabs(a[i]);
        if (max < v) { max = v; }
    }
    return max;
}
visible __forceinline__ double max(suNf_vector *r) {
    double *a = (double *)r;
    double max = 0.;
    for (int i = 0; i < sizeof(suNf_vector) / sizeof(*a); i++) {
        double v = fabs(a[i]);
        if (max < v) { max = v; }
    }
    return max;
}
visible __forceinline__ double max(suNg_vector *r) {
    double *a = (double *)r;
    double max = 0.;
    for (int i = 0; i < sizeof(suNg_vector) / sizeof(*a); i++) {
        double v = fabs(a[i]);
        if (max < v) { max = v; }
    }
    return max;
}
visible __forceinline__ double max(suNg_algebra_vector *r) {
    double *a = (double *)r;
    double max = 0.;
    for (int i = 0; i < sizeof(suNg_algebra_vector) / sizeof(*a); i++) {
        double v = fabs(a[i]);
        if (max < v) { max = v; }
    }
    return max;
}
visible __forceinline__ double max(double *r) {
    return *r;
}
visible __forceinline__ double max(float *r) {
    return (double)(*r);
}

visible __forceinline__ double g5_prod_re(suNf_spinor *s1, suNf_spinor *s2) {
    double k = 0;
    _spinor_g5_prod_re_f(k, *s1, *s2);
    return k;
}

visible __forceinline__ double g5_prod_re(suNf_spinor_flt *s1, suNf_spinor_flt *s2) {
    double z = 0;
    _spinor_g5_prod_re_f(z, *s1, *s2);
    return z;
}

visible __forceinline__ double g5_prod_im(suNf_spinor *s1, suNf_spinor *s2) {
    double z = 0;
    _spinor_g5_prod_im_f(z, *s1, *s2);
    return z;
}

visible __forceinline__ double g5_prod_im(suNf_spinor_flt *s1, suNf_spinor_flt *s2) {
    double z = 0;
    _spinor_g5_prod_im_f(z, *s1, *s2);
    return z;
}

visible __forceinline__ void g5(suNf_spinor *s1, suNf_spinor *s2) {
    _spinor_g5_f(*s1, *s2);
}

visible __forceinline__ void g5(suNf_spinor_flt *s1, suNf_spinor_flt *s2) {
    _spinor_g5_f(*s1, *s2);
}

visible __forceinline__ void g5_assign(suNf_spinor *s) {
    _spinor_g5_assign_f(*s);
}

visible __forceinline__ void g5_assign(suNf_spinor_flt *s) {
    _spinor_g5_assign_f(*s);
}

visible __forceinline__ void g0(suNf_spinor *s1, suNf_spinor *s2) {
    _spinor_g0_f(*s1, *s2);
}

visible __forceinline__ void g0(suNf_spinor_flt *s1, suNf_spinor_flt *s2) {
    _spinor_g0_f(*s1, *s2);
}

visible __forceinline__ void g1(suNf_spinor *s1, suNf_spinor *s2) {
    _spinor_g1_f(*s1, *s2);
}

visible __forceinline__ void g1(suNf_spinor_flt *s1, suNf_spinor_flt *s2) {
    _spinor_g1_f(*s1, *s2);
}

visible __forceinline__ void g2(suNf_spinor *s1, suNf_spinor *s2) {
    _spinor_g2_f(*s1, *s2);
}

visible __forceinline__ void g2(suNf_spinor_flt *s1, suNf_spinor_flt *s2) {
    _spinor_g2_f(*s1, *s2);
}

visible __forceinline__ void g3(suNf_spinor *s1, suNf_spinor *s2) {
    _spinor_g3_f(*s1, *s2);
}

visible __forceinline__ void g3(suNf_spinor_flt *s1, suNf_spinor_flt *s2) {
    _spinor_g3_f(*s1, *s2);
}

visible __forceinline__ void lc(suNf_spinor *r, double k1, suNf_spinor *s1, double k2, suNf_spinor *s2) {
    _spinor_lc_f(*r, k1, *s1, k2, *s2);
}

visible __forceinline__ void lc(suNf_spinor_flt *r, double k1, suNf_spinor_flt *s1, double k2, suNf_spinor_flt *s2) {
    _spinor_lc_f(*r, k1, *s1, k2, *s2);
}

visible __forceinline__ void lc_add_assign(suNf_spinor *r, double k1, suNf_spinor *s1, double k2, suNf_spinor *s2) {
    _spinor_lc_add_assign_f(*r, k1, *s1, k2, *s2);
}

visible __forceinline__ void lc_add_assign(suNf_spinor_flt *r, double k1, suNf_spinor_flt *s1, double k2, suNf_spinor_flt *s2) {
    _spinor_lc_add_assign_f(*r, k1, *s1, k2, *s2);
}

visible __forceinline__ void clc(suNf_spinor *r, hr_complex k1, suNf_spinor *s1, hr_complex k2, suNf_spinor *s2) {
    _spinor_clc_f(*r, k1, *s1, k2, *s2);
}

visible __forceinline__ void clc(suNf_spinor_flt *r, hr_complex k1, suNf_spinor_flt *s1, hr_complex k2, suNf_spinor_flt *s2) {
    _spinor_clc_f(*r, k1, *s1, k2, *s2);
}

visible __forceinline__ void clc_add_assign(suNf_spinor *r, hr_complex k1, suNf_spinor *s1, hr_complex k2, suNf_spinor *s2) {
    _spinor_clc_add_assign_f(*r, k1, *s1, k2, *s2);
}

visible __forceinline__ void clc_add_assign(suNf_spinor_flt *r, hr_complex k1, suNf_spinor_flt *s1, hr_complex k2,
                                            suNf_spinor_flt *s2) {
    _spinor_clc_add_assign_f(*r, k1, *s1, k2, *s2);
}

#undef _DECLARE_LINA_HEADER

#endif

#endif
