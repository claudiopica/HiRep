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
#define _BODY(a) a
#ifdef WITH_GPU
    #define _FUNC(a) a##_f_gpu
    #include "TMPL/linear_algebra.c.sdtmpl"
    #undef _FUNC
#endif
#define _FUNC(a,b,c) a b##_f_cpu c
#include "TMPL/linear_algebra.c.sdtmpl"
#undef _FUNC
#undef _BODY
#ifdef WITH_GPU
    #define _FUNC(a,b,c) b##_f = b##_f_gpu
#else
    #define _FUNC(a,b,c) b##_f = b##_f_cpu
#endif
#define _BODY(a) ;
#include "TMPL/linear_algebra.c.sdtmpl"
#undef _FUNC
#undef _BODY
#undef _SPINOR_FIELD_TYPE
#undef _REAL
#undef _COMPLEX

#if 0

/* single precision */
#define _SPINOR_FIELD_TYPE spinor_field_flt
#define _REAL float
#define _COMPLEX hr_complex_flt
#ifdef WITH_GPU
  #define _FUNC(a) a##_f_flt_gpu
  #include "TMPL/linear_algebra.c.sdtmpl"
  #undef _FUNC
#endif
#define _FUNC(a,b,c) a b##_f_flt_cpu c
#define _BODY(a) ;
#include "TMPL/linear_algebra.c.sdtmpl"
#undef _FUNC
#undef _BODY
#define _FUNC(a,b,c) a (*b##_f_flt) c
#define _BODY(a) ;
#include "TMPL/linear_algebra.c.sdtmpl"
#undef _FUNC
#undef _BODY
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

#endif
