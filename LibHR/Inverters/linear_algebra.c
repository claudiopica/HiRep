/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *
* All rights reserved.                                                      *
\***************************************************************************/

#include "global.h"
#include <string.h>
#include "spinor_field.h"
#include "gamma_spinor.h"
#include "hr_complex.h"
#include "communications.h"
#include "suN.h"
#include "suN_types.h"
#include "error.h"

/*
 * LINEAR ALGEBRA FUNCTIONS ARE DEFINED IN THE TEMPLATE
 *
 * TMPL/linear_algebra.c.sdtmpl
 *
 */

#ifdef WITH_GPU
  #include "TMPL/alloc_tmp_fields_gpu.c"
  #include "TMPL/global_sum_gpu.c"
#endif

/* double precision */
  // Declare types for double precision template parsing
#define _SPINOR_FIELD_TYPE spinor_field
#define _SPINOR_TYPE suNf_spinor
#define _REAL double
#define _COMPLEX hr_complex

  // Use GPU template to declare double precision functions (C++)
#ifdef WITH_GPU
  #define _FUNC(a,b,c) a b##_f_gpu c
  #define _BODY(a) a
  #include "TMPL/linear_algebra_gpu.c.sdtmpl"
  #undef _FUNC
  #undef _BODY
#endif

  // Use CPU templates for double precision functions & aliasing (C)
#ifdef __cplusplus
  extern "C" {
#endif
#define _FUNC(a,b,c) a b##_f_cpu c
#define _BODY(a) a
#include "TMPL/linear_algebra.c.sdtmpl"
#undef _FUNC
#undef _BODY
#ifdef WITH_GPU
  #define _FUNC(a,b,c) a (*b##_f) c = b##_f_gpu
#else
  #define _FUNC(a,b,c) a (*b##_f) c = b##_f_cpu
#endif
#define _BODY(a) ;
#include "TMPL/linear_algebra.c.sdtmpl"
#ifdef __cplusplus
  }
#endif

  //Undefine all definitions to move on to single precision
#undef _FUNC
#undef _BODY
#undef _SPINOR_FIELD_TYPE
#undef _SPINOR_TYPE
#undef _REAL
#undef _COMPLEX

/* single precision */
  //Declare types for single precision
#define _SPINOR_FIELD_TYPE spinor_field_flt
#define _SPINOR_TYPE suNf_spinor_flt
#define _REAL float
#define _COMPLEX hr_complex_flt

  // C++ GPU code
#ifdef WITH_GPU
  #define _FUNC(a,b,c) a b##_f_flt_gpu c
  #define _BODY(a) a
  #include "TMPL/linear_algebra_gpu.c.sdtmpl"
  #undef _FUNC
  #undef _BODY
#endif

  // C CPU code + aliasing
#ifdef __cplusplus
  extern "C" {
#endif 
#define _FUNC(a,b,c) a b##_f_flt_cpu c
#define _BODY(a) a
#include "TMPL/linear_algebra.c.sdtmpl"
#undef _FUNC
#undef _BODY

#ifdef WITH_GPU
  #define _FUNC(a,b,c) a (*b##_f_flt) c = b##_f_flt_gpu
#else
  #define _FUNC(a,b,c) a (*b##_f_flt) c = b##_f_flt_cpu
#endif
#define _BODY(a) ;
#include "TMPL/linear_algebra.c.sdtmpl"
#ifdef __cplusplus
  }
#endif

  //Undefine single precision definitions.
#undef _FUNC
#undef _BODY
#undef _SPINOR_FIELD_TYPE
#undef _SPINOR_TYPE
#undef _REAL
#undef _COMPLEX
