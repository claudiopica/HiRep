#ifndef LINEAR_ALGEBRA_H
#define LINEAR_ALGEBRA_H

/* 
 * LINEAR ALGEBRA FUNCTIONS ARE DEFINED IN THE TEMPLATE
 *
 * TMPL/linear_algebra.h.sdtmpl
 *
 */

#include "suN.h"

void set_spinor_len(unsigned int len);

/* single precision */
#define _SPINOR_TYPE suNf_spinor
#define _FUNC(a) a##_f
#include "TMPL/linear_algebra.h.sdtmpl"
#undef _SPINOR_TYPE
#undef _FUNC

/* double precision */
#define _SPINOR_TYPE suNf_spinor_dble
#define _FUNC(a) a##_dble_f
#include "TMPL/linear_algebra.h.sdtmpl"
#undef _SPINOR_TYPE
#undef _FUNC

#endif
