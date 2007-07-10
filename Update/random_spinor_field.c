#include "global.h"
#include "suN.h"
#include "random.h"
#include "linear_algebra.h"
#include "update.h"
#include <math.h>

void gaussian_spinor_field(suNf_spinor *s) {
   const float c1=1./sqrt(2.);
	 unsigned int len;
	 get_spinor_len(&len);
   gauss((float*)s,(sizeof(suNf_spinor)/sizeof(float))*len);
   spinor_field_mul_f(s,c1,s);
}
