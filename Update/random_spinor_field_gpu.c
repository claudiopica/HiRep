#include "global.h"
#include "dirac.h"
#include "suN.h"
#include "random.h"
#include "linear_algebra.h"
#include "update.h"
#include "memory.h"
#include <math.h>

#ifdef WITH_GPU

void gaussian_spinor_field(spinor_field *s) {
  gaussian_spinor_field_cpu(s);
  spinor_field_copy_to_gpu_f(s);
}


void gaussian_spinor_field_flt(spinor_field_flt *s) {
  gaussian_spinor_field_flt_cpu(s);
  spinor_field_copy_to_gpu_f_flt(s);
}

#endif
