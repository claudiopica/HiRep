/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *
* All rights reserved.                                                      *
\***************************************************************************/

#ifndef LINEAR_ALGEBRA_H
#define LINEAR_ALGEBRA_H

/*
 * LINEAR ALGEBRA FUNCTIONS ARE DEFINED IN THE TEMPLATE
 *
 * TMPL/linear_algebra.h.sdtmpl
 *
 */

#include "suN.h"
#include "spinor_field.h"
#include "hr_complex.h"

/* double precision */
#define _SPINOR_FIELD_TYPE spinor_field
#define _REAL double
#define _COMPLEX hr_complex
#ifdef WITH_GPU
  #define _FUNC(a) a##_f_gpu
  #include "TMPL/linear_algebra.h.sdtmpl"
  #undef _FUNC
#endif
#define _FUNC(a) a##_f_cpu
#include "TMPL/linear_algebra.h.sdtmpl"
#undef _FUNC
#define _FUNC(a) a##_f
#include "TMPL/linear_algebra.h.sdtmpl"
#undef _FUNC
#undef _SPINOR_FIELD_TYPE
#undef _REAL
#undef _COMPLEX

/* single precision */
#define _SPINOR_FIELD_TYPE spinor_field_flt
#define _REAL float
#define _COMPLEX hr_complex_flt
#ifdef WITH_GPU
  #define _FUNC(a) a##_f_flt_gpu
  #include "TMPL/linear_algebra.h.sdtmpl"
  #undef _FUNC
#endif
#define _FUNC(a) a##_f_flt_cpu
#include "TMPL/linear_algebra.h.sdtmpl"
#undef _FUNC
#define _FUNC(a) a##_f_flt
#include "TMPL/linear_algebra.h.sdtmpl"
#undef _FUNC
#undef _SPINOR_FIELD_TYPE
#undef _REAL
#undef _COMPLEX

/* GPU functions*/
#ifdef WITH_GPU
  #ifdef __cplusplus
  template <class T>
  T global_sum_gpu(T *vector, int size);
  extern "C" {
  #endif
  int global_sum_gpu_int(int *vector, int size);
  float global_sum_gpu_float(float *vector, int size);
  double global_sum_gpu_double(double *vector, int size);
  hr_complex_flt global_sum_gpu_complex_flt(hr_complex_flt *vector, int size);
  hr_complex global_sum_gpu_complex(hr_complex *vector, int size);
  #ifdef __cplusplus
  }
  #endif
#endif

#endif
