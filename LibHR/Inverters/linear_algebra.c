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
#undef _SPINOR_TYPE
#undef _REAL
#undef _COMPLEX


/* single precision */
#define _SPINOR_FIELD_TYPE spinor_field_flt
#define _SPINOR_TYPE suNf_spinor_flt
#define _REAL float
#define _COMPLEX hr_complex_flt
#define _FUNC(a,b,c) a b##_f_flt c
#define _BODY(a) a
#ifdef WITH_GPU
    #include "TMPL/linear_algebra_gpu.c.sdtmpl"
#else
    #include "TMPL/linear_algebra.c.sdtmpl"
#endif
#undef _FUNC
#undef _BODY
#define _FUNC(a,b,c) a b##_f_flt_cpu c_
#define _BODY(a) a
#include "TMPL/linear_algebra.c.sdtmpl"
#undef _FUNC
#undef _BODY
#undef _SPINOR_FIELD_TYPE
#undef _SPINOR_TYPE
#undef _REAL
#undef _COMPLEX
