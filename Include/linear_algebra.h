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

/* double precision */
#define _SPINOR_FIELD_TYPE spinor_field
#ifdef WITH_GPU
#define _FUNC(a) a##_f_cpu
#include "TMPL/linear_algebra.h.sdtmpl"
#undef _FUNC
#endif

#define _FUNC(a) a##_f
#include "TMPL/linear_algebra.h.sdtmpl"

#undef _SPINOR_FIELD_TYPE
#undef _FUNC

/* single precision */
#define _SPINOR_FIELD_TYPE spinor_field_flt
#ifdef WITH_GPU
#define _FUNC(a) a##_f_flt_cpu
#include "TMPL/linear_algebra.h.sdtmpl"
#undef _FUNC
#endif

#define _FUNC(a) a##_f_flt
#include "TMPL/linear_algebra.h.sdtmpl"

#undef _SPINOR_FIELD_TYPE
#undef _FUNC

#endif
