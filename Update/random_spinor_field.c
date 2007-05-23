#include "global.h"
#include "suN.h"
#include "random.h"
#include "linear_algebra.h"
#include "update.h"
#include <math.h>

void gaussian_spinor_field(suNf_spinor *s) {
   const float c1=1./sqrt(2.);
   gauss((float*)s,(sizeof(suNf_spinor)/sizeof(float))*VOLUME);
   spinor_field_mul_f(s,c1,s);
}
