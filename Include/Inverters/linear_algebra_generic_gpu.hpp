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

// in CUDA kernels, functions are inlined by the compiler (visible but not inline)
visible void mul_add_assign(suNf_spinor *s1, double rho, suNf_spinor *s2);
visible void mul_add_assign(suNf_spinor_flt *s1, double rho, suNf_spinor_flt *s2);
visible void mul_add_assign(suNf *s1, double rho, suNf *s2);
#ifdef REPR_IS_REAL
visible void mul_add_assign(suNfc *s1, double rho, suNfc *s2);
#endif
visible void mul_add_assign(suNg *s1, double rho, suNg *s2);
visible void mul_add_assign(suNf_flt *s1, double rho, suNf_flt *s2);
visible void mul_add_assign(suNg_flt *s1, double rho, suNg_flt *s2);
visible void mul_add_assign(suNf_vector *s1, double rho, suNf_vector *s2);
visible void mul_add_assign(suNg_vector *s1, double rho, suNg_vector *s2);
visible void mul_add_assign(suNg_algebra_vector *s1, double rho, suNg_algebra_vector *s2);
visible void mul_add_assign(double *s1, double rho, double *s2);
visible void mul_add_assign(float *s1, float rho, float *s2);

visible void mulc_add_assign(suNf_spinor *s1, hr_complex rho, suNf_spinor *s2);
visible void mulc_add_assign(suNf_spinor_flt *s1, hr_complex_flt rho, suNf_spinor_flt *s2);
visible void mulc_add_assign(suNf *s1, hr_complex rho, suNf *s2);
#ifdef REPR_IS_REAL
visible void mulc_add_assign(suNfc *s1, hr_complex rho, suNfc *s2);
#endif
visible void mulc_add_assign(suNg *s1, hr_complex rho, suNg *s2);
visible void mulc_add_assign(suNf_flt *s1, hr_complex_flt rho, suNf_flt *s2);
visible void mulc_add_assign(suNg_flt *s1, hr_complex_flt rho, suNg_flt *s2);
visible void mulc_add_assign(suNf_vector *s1, hr_complex rho, suNf_vector *s2);
visible void mulc_add_assign(suNg_vector *s1, hr_complex rho, suNg_vector *s2);
visible void mulc_add_assign(suNg_algebra_vector *s1, hr_complex rho, suNg_algebra_vector *s2);
visible void mulc_add_assign(double *s1, hr_complex rho, double *s2);
visible void mulc_add_assign(float *s1, hr_complex_flt rho, float *s2);

visible void mul(suNf_spinor *s1, double rho, suNf_spinor *s2);
visible void mul(suNf_spinor_flt *s1, double rho, suNf_spinor_flt *s2);
visible void mul(suNf *s1, double rho, suNf *s2);
#ifdef REPR_IS_REAL
visible void mul(suNfc *s1, double rho, suNfc *s2);
#endif
visible void mul(suNg *s1, double rho, suNg *s2);
visible void mul(suNf_flt *s1, double rho, suNf_flt *s2);
visible void mul(suNg_flt *s1, double rho, suNg_flt *s2);
visible void mul(suNf_vector *s1, double rho, suNf_vector *s2);
visible void mul(suNg_vector *s1, double rho, suNg_vector *s2);
visible void mul(suNg_algebra_vector *s1, double rho, suNg_algebra_vector *s2);
visible void mul(double *s1, double rho, double *s2);
visible void mul(float *s1, float rho, float *s2);

visible void mulc(suNf_spinor *s1, hr_complex rho, suNf_spinor *s2);
visible void mulc(suNf_spinor_flt *s1, hr_complex_flt rho, suNf_spinor_flt *s2);
visible void mulc(suNf *s1, hr_complex rho, suNf *s2);
#ifdef REPR_IS_REAL
visible void mulc(suNfc *s1, hr_complex rho, suNfc *s2);
#endif
visible void mulc(suNg *s1, hr_complex rho, suNg *s2);
visible void mulc(suNf_flt *s1, hr_complex_flt rho, suNf_flt *s2);
visible void mulc(suNg_flt *s1, hr_complex_flt rho, suNg_flt *s2);
visible void mulc(suNf_vector *s1, hr_complex rho, suNf_vector *s2);
visible void mulc(suNg_vector *s1, hr_complex rho, suNg_vector *s2);
visible void mulc(suNg_algebra_vector *s1, hr_complex rho, suNg_algebra_vector *s2);
visible void mulc(double *s1, hr_complex rho, double *s2);
visible void mulc(float *s1, hr_complex_flt rho, float *s2);

visible void add(suNf_spinor *r, suNf_spinor *s1, suNf_spinor *s2);
visible void add(suNf_spinor_flt *r, suNf_spinor_flt *s1, suNf_spinor_flt *s2);
visible void add(suNf *r, suNf *s1, suNf *s2);
#ifdef REPR_IS_REAL
visible void add(suNfc *r, suNfc *s1, suNfc *s2);
#endif
visible void add(suNg *r, suNg *s1, suNg *s2);
visible void add(suNf_flt *r, suNf_flt *s1, suNf_flt *s2);
visible void add(suNg_flt *r, suNg_flt *s1, suNg_flt *s2);
visible void add(suNf_vector *r, suNf_vector *s1, suNf_vector *s2);
visible void add(suNg_vector *r, suNg_vector *s1, suNg_vector *s2);
visible void add(suNg_algebra_vector *r, suNg_algebra_vector *s1, suNg_algebra_vector *s2);
visible void add(double *r, double *s1, double *s2);
visible void add(float *r, float *s1, float *s2);

visible void sub(suNf_spinor *r, suNf_spinor *s1, suNf_spinor *s2);
visible void sub(suNf_spinor_flt *r, suNf_spinor_flt *s1, suNf_spinor_flt *s2);
visible void sub(suNf *r, suNf *s1, suNf *s2);
#ifdef REPR_IS_REAL
visible void sub(suNfc *r, suNfc *s1, suNfc *s2);
#endif
visible void sub(suNg *r, suNg *s1, suNg *s2);
visible void sub(suNf_flt *r, suNf_flt *s1, suNf_flt *s2);
visible void sub(suNg_flt *r, suNg_flt *s1, suNg_flt *s2);
visible void sub(suNf_vector *r, suNf_vector *s1, suNf_vector *s2);
visible void sub(suNg_vector *r, suNg_vector *s1, suNg_vector *s2);
visible void sub(suNg_algebra_vector *r, suNg_algebra_vector *s1, suNg_algebra_vector *s2);
visible void sub(double *r, double *s1, double *s2);
visible void sub(float *r, float *s1, float *s2);

visible void sub_assign(suNf_spinor *s1, suNf_spinor *s2);
visible void sub_assign(suNf_spinor_flt *s1, suNf_spinor_flt *s2);
visible void sub_assign(suNf *s1, suNf *s2);
#ifdef REPR_IS_REAL
visible void sub_assign(suNfc *s1, suNfc *s2);
#endif
visible void sub_assign(suNg *s1, suNg *s2);
visible void sub_assign(suNf_flt *s1, suNf_flt *s2);
visible void sub_assign(suNg_flt *s1, suNg_flt *s2);
visible void sub_assign(suNf_vector *s1, suNf_vector *s2);
visible void sub_assign(suNg_vector *s1, suNg_vector *s2);
visible void sub_assign(suNg_algebra_vector *s1, suNg_algebra_vector *s2);
visible void sub_assign(double *s1, double *s2);
visible void sub_assign(float *s1, float *s2);
visible void minus(suNf_spinor *s1, suNf_spinor *s2);
visible void minus(suNf_spinor_flt *s1, suNf_spinor_flt *s2);
visible void minus(suNf *s1, suNf *s2);
#ifdef REPR_IS_REAL
visible void minus(suNfc *s1, suNfc *s2);
#endif
visible void minus(suNg *s1, suNg *s2);
visible void minus(suNf_flt *s1, suNf_flt *s2);
visible void minus(suNg_flt *s1, suNg_flt *s2);
visible void minus(suNf_vector *s1, suNf_vector *s2);
visible void minus(suNg_vector *s1, suNg_vector *s2);
visible void minus(suNg_algebra_vector *s1, suNg_algebra_vector *s2);
visible void minus(double *s1, double *s2);
visible void minus(float *s1, float *s2);
visible void add_assign(suNf_spinor *s1, suNf_spinor *s2);
visible void add_assign(suNf_spinor_flt *s1, suNf_spinor_flt *s2);
visible void add_assign(suNf *s1, suNf *s2);
#ifdef REPR_IS_REAL
visible void add_assign(suNfc *s1, suNfc *s2);
#endif
visible void add_assign(suNg *s1, suNg *s2);
visible void add_assign(suNf_flt *s1, suNf_flt *s2);
visible void add_assign(suNg_flt *s1, suNg_flt *s2);
visible void add_assign(suNf_vector *s1, suNf_vector *s2);
visible void add_assign(suNg_vector *s1, suNg_vector *s2);
visible void add_assign(suNg_algebra_vector *s1, suNg_algebra_vector *s2);
visible void add_assign(double *s1, double *s2);
visible void add_assign(float *s1, float *s2);
visible double prod_re(suNf_spinor *s1, suNf_spinor *s2);
visible double prod_re(suNf_spinor_flt *s1, suNf_spinor_flt *s2);
visible double prod_re(suNf *s1, suNf *s2);
#ifdef REPR_IS_REAL
visible double prod_re(suNfc *s1, suNfc *s2);
#endif
visible double prod_re(suNg *s1, suNg *s2);
visible double prod_re(suNf_flt *s1, suNf_flt *s2);
visible double prod_re(suNg_flt *s1, suNg_flt *s2);
visible double prod_re(suNf_vector *s1, suNf_vector *s2);
visible double prod_re(suNg_vector *s1, suNg_vector *s2);
visible double prod_re(suNg_algebra_vector *s1, suNg_algebra_vector *s2);
visible double prod_re(double *s1, double *s2);
visible double prod_re(float *s1, float *s2);
visible double prod_im(suNf_spinor *s1, suNf_spinor *s2);
visible double prod_im(suNf_spinor_flt *s1, suNf_spinor_flt *s2);
visible double prod_im(suNf *s1, suNf *s2);
#ifdef REPR_IS_REAL
visible double prod_im(suNfc *s1, suNfc *s2);
#endif
visible double prod_im(suNg *s1, suNg *s2);
visible double prod_im(suNf_flt *s1, suNf_flt *s2);
visible double prod_im(suNg_flt *s1, suNg_flt *s2);
visible double prod_im(suNf_vector *s1, suNf_vector *s2);
visible double prod_im(suNg_vector *s1, suNg_vector *s2);
visible double prod_im(suNg_algebra_vector *s1, suNg_algebra_vector *s2);
visible double prod_im(double *s1, double *s2);
visible double prod_im(float *s1, float *s2);

visible hr_complex prod(suNf_spinor *s1, suNf_spinor *s2);
visible hr_complex prod(suNf_spinor_flt *s1, suNf_spinor_flt *s2);
visible hr_complex prod(suNf *s1, suNf *s2);
#ifdef REPR_IS_REAL
visible hr_complex prod(suNfc *s1, suNfc *s2);
#endif
visible hr_complex prod(suNg *s1, suNg *s2);
visible hr_complex prod(suNf_flt *s1, suNf_flt *s2);
visible hr_complex prod(suNg_flt *s1, suNg_flt *s2);
visible hr_complex prod(suNf_vector *s1, suNf_vector *s2);
visible hr_complex prod(suNg_vector *s1, suNg_vector *s2);
visible hr_complex prod(suNg_algebra_vector *s1, suNg_algebra_vector *s2);
visible hr_complex prod(double *s1, double *s2);
visible hr_complex prod(float *s1, float *s2);
visible double sqnorm(suNf_spinor *r);
visible double sqnorm(suNf_spinor_flt *r);
visible double sqnorm(suNf *r);
#ifdef REPR_IS_REAL
visible double sqnorm(suNfc *r);
#endif
visible double sqnorm(suNg *r);
visible double sqnorm(suNf_flt *r);
visible double sqnorm(suNg_flt *r);
visible double sqnorm(suNf_vector *r);
visible double sqnorm(suNg_vector *r);
visible double sqnorm(suNg_algebra_vector *r);
visible double sqnorm(double *r);
visible double sqnorm(float *r);

visible double max(suNf_spinor *r);
visible double max(suNf_spinor_flt *r);
visible double max(suNf *r);
#ifdef REPR_IS_REAL
visible double max(suNfc *r);
#endif
visible double max(suNg *r);
visible double max(suNf_flt *r);
visible double max(suNg_flt *r);
visible double max(suNf_vector *r);
visible double max(suNg_vector *r);
visible double max(suNg_algebra_vector *r);
visible double max(double *r);
visible double max(float *r);

visible double g5_prod_re(suNf_spinor *s1, suNf_spinor *s2);
visible double g5_prod_re(suNf_spinor_flt *s1, suNf_spinor_flt *s2);

visible double g5_prod_im(suNf_spinor *s1, suNf_spinor *s2);
visible double g5_prod_im(suNf_spinor_flt *s1, suNf_spinor_flt *s2);

visible void g5(suNf_spinor *s1, suNf_spinor *s2);
visible void g5(suNf_spinor_flt *s1, suNf_spinor_flt *s2);

visible void g5_assign(suNf_spinor *s);
visible void g5_assign(suNf_spinor_flt *s);

visible void g0(suNf_spinor *s1, suNf_spinor *s2);
visible void g0(suNf_spinor_flt *s1, suNf_spinor_flt *s2);

visible void g1(suNf_spinor *s1, suNf_spinor *s2);
visible void g1(suNf_spinor_flt *s1, suNf_spinor_flt *s2);

visible void g2(suNf_spinor *s1, suNf_spinor *s2);
visible void g2(suNf_spinor_flt *s1, suNf_spinor_flt *s2);

visible void g3(suNf_spinor *s1, suNf_spinor *s2);
visible void g3(suNf_spinor_flt *s1, suNf_spinor_flt *s2);

visible void lc(suNf_spinor *r, double k1, suNf_spinor *s1, double k2, suNf_spinor *s2);
visible void lc(suNf_spinor_flt *r, double k1, suNf_spinor_flt *s1, double k2, suNf_spinor_flt *s2);

visible void lc_add_assign(suNf_spinor *r, double k1, suNf_spinor *s1, double k2, suNf_spinor *s2);
visible void lc_add_assign(suNf_spinor_flt *r, double k1, suNf_spinor_flt *s1, double k2, suNf_spinor_flt *s2);

visible void clc(suNf_spinor *r, hr_complex k1, suNf_spinor *s1, hr_complex k2, suNf_spinor *s2);
visible void clc(suNf_spinor_flt *r, hr_complex_flt k1, suNf_spinor_flt *s1, hr_complex_flt k2, suNf_spinor_flt *s2);

visible void clc_add_assign(suNf_spinor *r, hr_complex k1, suNf_spinor *s1, hr_complex k2, suNf_spinor *s2);
visible void clc_add_assign(suNf_spinor_flt *r, hr_complex_flt k1, suNf_spinor_flt *s1, hr_complex_flt k2, suNf_spinor_flt *s2);

#undef _DECLARE_LINA_HEADER

#endif

#endif
