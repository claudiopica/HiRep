#include "linear_algebra.h"
#include "global.h"
#include <string.h>

/* 
 * LINEAR ALGEBRA FUNCTIONS ARE DEFINED IN THE TEMPLATE
 *
 * TMPL/linear_algebra.c.sdtmpl
 *
 */


static unsigned int _spinor_len=VOLUME;

void set_spinor_len(unsigned int len) {
  _spinor_len=len;
}

void get_spinor_len(unsigned int *len) {
  *len=_spinor_len;
}

/* single precision */
#define _SPINOR_TYPE suNf_spinor
#define _FUNC(a) a##_f
#include "TMPL/linear_algebra.c.sdtmpl"
#undef _SPINOR_TYPE
#undef _FUNC

/* double precision */
#define _SPINOR_TYPE suNf_spinor_dble
#define _FUNC(a) a##_dble_f
#include "TMPL/linear_algebra.c.sdtmpl"
#undef _SPINOR_TYPE
#undef _FUNC
