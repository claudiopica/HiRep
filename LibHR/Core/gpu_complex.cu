#include "libhr_core.h"

#define PI 3.141592653589793238462643383279502884197

visible double carg(hr_complex c) {
  double arg = 0;
  if (creal(c) > 0) {
    arg = atan(cimag(c) / creal(c));
  } else if (creal(c) < 0 && cimag(c) >= 0) {
    arg = atan(cimag(c) / creal(c)) + PI;
  } else if (creal(c) < 0 && cimag(c) < 0) {
    arg = atan(cimag(c) / creal(c)) - PI;
  } else if (creal(c) == 0 && cimag(c) > 0) {
    arg = PI / 2.0;
  } else if (creal(c) == 0 && cimag(c) < 0) {
    arg = - PI / 2.0;
  } else if (creal(c) == 0 && creal(c) == 0) {
    // This is technically undefined, so this will
    // give 0 in the hopes that we get 0 by the
    // mod and does not quit with error
    arg = 0.0;
  }
  return arg;
}

visible hr_complex cpow(hr_complex c, double pow) {

    double arg = carg(c);
    double mod = sqrt(creal(_complex_prod(c,c)));
    return hr_complex(powf(mod, pow),0) * hr_complex(cos(pow*arg),sin(pow*arg));
}
