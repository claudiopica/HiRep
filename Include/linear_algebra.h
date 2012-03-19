/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *   
* All rights reserved.                                                      * 
\***************************************************************************/

#ifndef LINEAR_ALGEBRA_H
#define LINEAR_ALGEBRA_H

#ifdef __cplusplus
extern "C" {
#endif

/* 
 * LINEAR ALGEBRA FUNCTIONS ARE DEFINED IN THE TEMPLATE
 *
 * TMPL/linear_algebra.h.sdtmpl
 *
 */

#include "suN.h"
#include "spinor_field.h"

double global_sum_gpu(double* vector, int n);


/* double precision */
#define _SPINOR_FIELD_TYPE spinor_field

#define _FUNC(a) a##_f_cpu
#include "TMPL/linear_algebra.h.sdtmpl"
#undef _FUNC

#define _FUNC(a) a##_f
#include "TMPL/linear_algebra.h.sdtmpl"

#undef _SPINOR_FIELD_TYPE
#undef _FUNC

/* single precision */
#define _SPINOR_FIELD_TYPE spinor_field_flt

#define _FUNC(a) a##_f_flt_cpu
#include "TMPL/linear_algebra.h.sdtmpl"
#undef _FUNC

#define _FUNC(a) a##_f_flt
#include "TMPL/linear_algebra.h.sdtmpl"

#undef _SPINOR_FIELD_TYPE
#undef _FUNC

#ifdef __cplusplus
}
#endif
	
	
#endif


