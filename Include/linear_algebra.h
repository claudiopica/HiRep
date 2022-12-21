/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *
* All rights reserved.                                                      *
\***************************************************************************/

#ifndef LINEAR_ALGEBRA_H
#define LINEAR_ALGEBRA_H
#include "hr_complex.h"
#ifdef __cplusplus
  extern "C" {
#endif

/*
 * LINEAR ALGEBRA FUNCTIONS ARE DEFINED IN THE TEMPLATE
 *
 * TMPL/linear_algebra.c.sdtmpl
 *
 */
#include "suN.h"
#include "spinor_field.h"
#include "reduction.h"

/* double precision */
#define _SPINOR_FIELD_TYPE spinor_field
#define _REAL double
#define _COMPLEX hr_complex

// Declare GPU functions
#ifdef WITH_GPU
  #define _FUNC(a,b,c) a b##_f_gpu c
  #define _BODY(a) ;
  #include "TMPL/linear_algebra.c.sdtmpl"
  #undef _FUNC
  #undef _BODY
#endif

// Declare CPU functions
#define _FUNC(a,b,c) a b##_f_cpu c
#define _BODY(a) ;
#include "TMPL/linear_algebra.c.sdtmpl"
#undef _FUNC
#undef _BODY

// Declare alias function pointers
#define _FUNC(a,b,c) extern a (*b##_f) c
#define _BODY(a) ;
#include "TMPL/linear_algebra.c.sdtmpl"
#undef _FUNC
#undef _BODY

#undef _SPINOR_FIELD_TYPE
#undef _REAL
#undef _COMPLEX

/* single precision */
#define _SPINOR_FIELD_TYPE spinor_field_flt
#define _REAL float
#define _COMPLEX hr_complex_flt

// Declare GPU functions
#ifdef WITH_GPU
  #define _FUNC(a,b,c) a b##_f_flt_gpu c
  #define _BODY(a) ;
  #include "TMPL/linear_algebra.c.sdtmpl"
  #undef _FUNC
  #undef _BODY
#endif

// Declare CPU functions
#define _FUNC(a,b,c) a b##_f_flt_cpu c
#define _BODY(a) ;
#include "TMPL/linear_algebra.c.sdtmpl"
#undef _FUNC
#undef _BODY

// Declare alias function pointers
#define _FUNC(a,b,c) extern a (*b##_f_flt) c
#define _BODY(a) ;
#include "TMPL/linear_algebra.c.sdtmpl"
#undef _FUNC
#undef _BODY

#undef _SPINOR_FIELD_TYPE
#undef _REAL
#undef _COMPLEX

#ifdef __cplusplus
  }
#endif
#endif
