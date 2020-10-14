/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *   
* All rights reserved.                                                      * 
\***************************************************************************/

#include "linear_algebra.h"
#include "global.h"
#include <string.h>
#include "spinor_field.h"
#include "gamma_spinor.h"

/* 
 * LINEAR ALGEBRA FUNCTIONS ARE DEFINED IN THE TEMPLATE
 *
 * TMPL/linear_algebra.c.sdtmpl
 *
 */


/* double precision */

#define _SPINOR_FIELD_TYPE spinor_field
#define _SPINOR_TYPE suNf_spinor
#define _REAL double
#define _COMPLEX hr_complex
#define _FUNC(a) a##_f
#ifdef WITH_GPU
#include "TMPL/linear_algebra_gpu.c.sdtmpl"
#undef _FUNC
#define _FUNC(a) a##_f_cpu
#endif
#include "TMPL/linear_algebra.c.sdtmpl"
#undef _FUNC
#undef _SPINOR_FIELD_TYPE
#undef _SPINOR_TYPE
#undef _REAL
#undef _COMPLEX


/* single precision */

#define _SPINOR_FIELD_TYPE spinor_field_flt
#define _SPINOR_TYPE suNf_spinor_flt
#define _REAL float
#define _COMPLEX hr_complex_flt
#define _FUNC(a) a##_f_flt
#ifdef WITH_GPU
#include "TMPL/linear_algebra_gpu.c.sdtmpl"
#undef _FUNC
#define _FUNC(a) a##_f_flt_cpu
#endif
#include "TMPL/linear_algebra.c.sdtmpl"
#undef _FUNC
#undef _SPINOR_FIELD_TYPE
#undef _SPINOR_TYPE
#undef _REAL
#undef _COMPLEX




