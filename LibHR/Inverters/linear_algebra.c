/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *
* All rights reserved.                                                      *
\***************************************************************************/

#include "linear_algebra.h"
#include "global.h"
#include <string.h>
#include "spinor_field.h"
#include "gamma_spinor.h"
#include "hr_complex.h"
#include "communications.h"
#include "suN.h"
#include "suN_types.h"

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
#define _SPINOR_FIELD_TYPE spinor_field
#define _SPINOR_TYPE suNf_spinor
#define _REAL double
#define _COMPLEX hr_complex
#ifdef WITH_GPU
  #define _FUNC(a2,b2,c2) a2 b2##_f_gpu c2
  #define _BODY(a2) a2
  #include "TMPL/linear_algebra_gpu.c.sdtmpl"
  #undef _FUNC
  #undef _BODY
#endif
#define _FUNC(a2,b2,c2) a2 b2##_f_cpu c2
#define _BODY(a2) a2
#include "TMPL/linear_algebra.c.sdtmpl"
#undef _FUNC
#undef _BODY
#ifdef WITH_GPU
  #define _FUNC(a2,b2,c2) a2 (*b2##_f) c2 = b2##_f_gpu
#else
  #define _FUNC(a2,b2,c2) a2 (*b2##_f) c2 = b2##_f_cpu
#endif
#define _BODY(a2) ;
#include "TMPL/linear_algebra.c.sdtmpl"
#undef _FUNC
#undef _BODY
#undef _SPINOR_FIELD_TYPE
#undef _SPINOR_TYPE
#undef _REAL
#undef _COMPLEX

/* single precision */
#define _SPINOR_FIELD_TYPE spinor_field_flt
#define _SPINOR_TYPE suNf_spinor_flt
#define _REAL float
#define _COMPLEX hr_complex_flt
#ifdef WITH_GPU
  #define _FUNC(a1,b1,c1) a1 b1##_f_flt_gpu c1
  #define _BODY(a1) a1
  #include "TMPL/linear_algebra_gpu.c.sdtmpl"
  #undef _FUNC
  #undef _BODY
#endif
#define _FUNC(a1,b1,c1) a1 b1##_f_flt_cpu c1
#define _BODY(a1) a1
#include "TMPL/linear_algebra.c.sdtmpl"
#undef _FUNC
#undef _BODY
#ifdef WITH_GPU
  #define _FUNC(a1,b1,c1) a1 (*b1##_f_flt) c1 = b1##_f_flt_gpu
#else
  #define _FUNC(a1,b1,c1) a1 (*b1##_f_flt) c1 = b1##_f_flt_cpu
#endif
#define _BODY(a1) ;
#include "TMPL/linear_algebra.c.sdtmpl"
#undef _FUNC
#undef _BODY
#undef _SPINOR_FIELD_TYPE
#undef _SPINOR_TYPE
#undef _REAL
#undef _COMPLEX
