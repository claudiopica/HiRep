#include "global.h"
#include "suN.h"
#include "random.h"
#include "linear_algebra.h"
#include "update.h"
#include <math.h>

void gaussian_spinor_field(spinor_field *s) {
	const double c1=1./sqrt(2.);
	unsigned int ix,sx;
/*	get_spinor_len(&len); */
   FOR_LOCAL_SC(ix,sx) {
      suNf_spinor* sptr=_SPINOR_AT(s,sx);
 	  gauss((double*)sptr,sizeof(suNf_spinor)/sizeof(double));
   }
	spinor_field_mul_f(s,c1,s);
}

void gaussian_spinor_field_flt(spinor_field_flt *s) {
	const float c1=1./sqrt(2.);
	unsigned int ix,sx;
/*	get_spinor_len(&len); */
   FOR_LOCAL_SC(ix,sx) {
      suNf_spinor_flt* sptr=_SPINOR_AT(s,sx);
 	  gauss_flt((float*)sptr,sizeof(suNf_spinor)/sizeof(float));
   }
	spinor_field_mul_f_flt(s,c1,s);
}
