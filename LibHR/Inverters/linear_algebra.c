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
#define _FUNC(a) a##_f
#define _REAL double
#define _COMPLEX complex

#include "TMPL/linear_algebra.c.sdtmpl"
void _FUNC(spinor_field_copy)(_SPINOR_FIELD_TYPE *s1, _SPINOR_FIELD_TYPE *s2) {
	_TWO_SPINORS_MATCHING(s1,s2);
	memcpy(s1->ptr,s2->ptr,s1->type->gsize_spinor*sizeof(suNf_spinor));
}

#undef _SPINOR_FIELD_TYPE
#undef _SPINOR_TYPE
#undef _FUNC
#undef _REAL
#undef _COMPLEX


/* single precision */

#define _SPINOR_FIELD_TYPE spinor_field_flt
#define _SPINOR_TYPE suNf_spinor_flt
#define _FUNC(a) a##_f_flt
#define _REAL float
#define _COMPLEX complex_flt

#include "TMPL/linear_algebra.c.sdtmpl"
void _FUNC(spinor_field_copy)(_SPINOR_FIELD_TYPE *s1, _SPINOR_FIELD_TYPE *s2) {
	_TWO_SPINORS_MATCHING(s1,s2);
	memcpy(s1->ptr,s2->ptr,s1->type->gsize_spinor*sizeof(suNf_spinor_flt));
}

#undef _SPINOR_FIELD_TYPE
#undef _SPINOR_TYPE
#undef _FUNC
#undef _REAL
#undef _COMPLEX


