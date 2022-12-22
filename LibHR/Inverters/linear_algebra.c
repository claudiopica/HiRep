/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *
* All rights reserved.                                                      *
\***************************************************************************/

#include "global.h"
#include <string.h>
#include "spinor_field.h"
#include "gamma_spinor.h"
#include "hr_complex.h"
// #include "communications.h"
#include "suN.h"
#include "suN_types.h"
#include "error.h"
#include "linear_algebra.h"

/*
 * LINEAR ALGEBRA FUNCTIONS ARE DEFINED IN THE TEMPLATE
 *
 * TMPL/linear_algebra.c.sdtmpl
 *
 */

/* double precision */
// Declare types for double precision template parsing
#define _SPINOR_FIELD_TYPE spinor_field
#define _SPINOR_TYPE suNf_spinor
#define _REAL double
#define _COMPLEX hr_complex

  // Use CPU templates for double precision functions & aliasing (C)
#define _FUNC(a,b,c) a b##_f_cpu c
#define _BODY(a) a
#include "TMPL/linear_algebra.c.sdtmpl"
#undef _FUNC
#undef _BODY

//set double precision function aliases
#ifdef WITH_GPU
  #define _FUNC(a,b,c) a (*b##_f) c = b##_f_gpu
#else
  #define _FUNC(a,b,c) a (*b##_f) c = b##_f_cpu
#endif
#define _BODY(a) ;
#include "TMPL/linear_algebra.c.sdtmpl"

 //Undefine all definitions to move on to single precision
#undef _FUNC
#undef _BODY
#undef _SPINOR_FIELD_TYPE
#undef _SPINOR_TYPE
#undef _REAL
#undef _COMPLEX

#define _SPINOR_FIELD_TYPE spinor_field_flt
#define _SPINOR_TYPE suNf_spinor_flt
#define _REAL float
#define _COMPLEX hr_complex_flt

#define _FUNC(a,b,c) a b##_f_flt_cpu c
#define _BODY(a) a
#include "TMPL/linear_algebra.c.sdtmpl"
#undef _FUNC
#undef _BODY

//set single precision function aliases
#ifdef WITH_GPU
  #define _FUNC(a,b,c) a (*b##_f_flt) c = b##_f_flt_gpu
#else
  #define _FUNC(a,b,c) a (*b##_f_flt) c = b##_f_flt_cpu
#endif
#define _BODY(a) ;
#include "TMPL/linear_algebra.c.sdtmpl"

#undef _FUNC
#undef _BODY
#undef _SPINOR_FIELD_TYPE
#undef _SPINOR_TYPE
#undef _REAL
#undef _COMPLEX
