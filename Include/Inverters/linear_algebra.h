/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *
* All rights reserved.                                                      *
\***************************************************************************/

#ifndef LINEAR_ALGEBRA_H
#define LINEAR_ALGEBRA_H

#include "hr_complex.h"
#include "suN.h"
#include "spinor_field.h"

#ifdef __cplusplus
  extern "C" {
#endif

/*
 * LINEAR ALGEBRA FUNCTIONS ARE DEFINED IN THE TEMPLATE
 *
 * TMPL/linear_algebra.c.sdtmpl
 *
 */

/* double precision */
#define _SPINOR_FIELD_TYPE spinor_field
#define _REAL double
#define _COMPLEX hr_complex

// GPU functions
#ifdef WITH_GPU
#define _GPU_FUNC(a,b,c) a b##_f_gpu c ;
#else 
#define _GPU_FUNC(a,b,c) 
#endif

// Declare functions and alias function pointers
#define _FUNC(a,b,c) _GPU_FUNC(a,b,c) a b##_f_cpu c; extern a (*b##_f) c
#include "./TMPL/linear_algebra.h.sdtmpl"
#undef _FUNC
#undef _GPU_FUNC

#undef _SPINOR_FIELD_TYPE
#undef _REAL
#undef _COMPLEX

/* single precision */
#define _SPINOR_FIELD_TYPE spinor_field_flt
#define _REAL float
#define _COMPLEX hr_complex_flt

// Declare GPU functions
#ifdef WITH_GPU
#define _GPU_FUNC(a,b,c) a b##_f_flt_gpu c ;
#else 
#define _GPU_FUNC(a,b,c) 
#endif

// Declare CPU functions
#define _FUNC(a,b,c) _GPU_FUNC(a,b,c) a b##_f_flt_cpu c; extern a (*b##_f_flt) c
#include "./TMPL/linear_algebra.h.sdtmpl"
#undef _FUNC
#undef _GPU_FUNC

#undef _SPINOR_FIELD_TYPE
#undef _REAL
#undef _COMPLEX

#ifdef __cplusplus
  }
#endif
#endif
