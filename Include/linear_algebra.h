#ifndef LINEAR_ALGEBRA_H
#define LINEAR_ALGEBRA_H

/* 
 * LINEAR ALGEBRA FUNCTIONS ARE DEFINED IN THE TEMPLATE
 *
 * TMPL/linear_algebra.h.sdtmpl
 *
 */

#include "suN.h"


/*
typedef struct {
	enum { ALL=0 };
	suNf_spinor *array;
} spinor_field;

typedef struct {
	enum { ALL=0 };
	suNf_spinor_flt *array;
} spinor_field_flt;
*/


/* double precision */
#define _SPINOR_FIELD_TYPE spinor_field
#define _FUNC(a) a##_f
#include "TMPL/linear_algebra.h.sdtmpl"
#undef _SPINOR_FIELD_TYPE
#undef _FUNC

/* single precision */
#define _SPINOR_FIELD_TYPE spinor_field_flt
#define _FUNC(a) a##_f_flt
#include "TMPL/linear_algebra.h.sdtmpl"
#undef _SPINOR_FIELD_TYPE
#undef _FUNC

#endif
