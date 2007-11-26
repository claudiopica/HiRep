#include "linear_algebra.h"
#include "global.h"
#include <string.h>

/* 
 * LINEAR ALGEBRA FUNCTIONS ARE DEFINED IN THE TEMPLATE
 *
 * TMPL/linear_algebra.c.sdtmpl
 *
 */


static unsigned int _spinor_len=0;

void set_spinor_len(unsigned int len) {
  _spinor_len=len;
}

void get_spinor_len(unsigned int *len) {
  *len=_spinor_len;
}


/* double precision */

#define _SPINOR_FIELD_TYPE spinor_field
#define _SPINOR_TYPE suNf_spinor
#define _FUNC(a) a##_f

#include "TMPL/linear_algebra.c.sdtmpl"
void _FUNC(spinor_field_copy)(_SPINOR_FIELD_TYPE *s1, _SPINOR_FIELD_TYPE *s2) {
	memcpy(s1->ptr,s2->ptr,_spinor_len*sizeof(suNf_spinor));
}

#undef _SPINOR_FIELD_TYPE
#undef _SPINOR_TYPE
#undef _FUNC


/* single precision */

#define _SPINOR_FIELD_TYPE spinor_field_flt
#define _SPINOR_TYPE suNf_spinor_flt
#define _FUNC(a) a##_f_flt

#include "TMPL/linear_algebra.c.sdtmpl"
void _FUNC(spinor_field_copy)(_SPINOR_FIELD_TYPE *s1, _SPINOR_FIELD_TYPE *s2) {
	memcpy(s1->ptr,s2->ptr,_spinor_len*sizeof(suNf_spinor));
}

#undef _SPINOR_FIELD_TYPE
#undef _SPINOR_TYPE
#undef _FUNC
