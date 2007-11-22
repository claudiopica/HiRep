#include "global.h"
#include "suN.h"
#include "random.h"
#include "linear_algebra.h"
#include "update.h"
#include <math.h>

void gaussian_spinor_field(spinor_field *s) {
	const double c1=1./sqrt(2.);
	unsigned int len;
	get_spinor_len(&len);
	gauss((double*)(*s),(sizeof(suNf_spinor)/sizeof(double))*len);
	spinor_field_mul_f(s,c1,s);
}

void gaussian_spinor_field_flt(spinor_field_flt *s) {
	const float c1=1./sqrt(2.);
	unsigned int len;
	get_spinor_len(&len);
	gauss_flt((float*)(*s),(sizeof(suNf_spinor_flt)/sizeof(float))*len);
	spinor_field_mul_f_flt(s,c1,s);
}
