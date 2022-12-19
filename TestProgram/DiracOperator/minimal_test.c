#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "io.h"
#include "random.h"
#include "error.h"
#include "global.h"
#include "geometry.h"
#include "memory.h"
#include "update.h"
#include "suN.h"
#include "suN_types.h"
#include "linear_algebra.h"
#include "inverters.h"
#include "representation.h"
#include "utils.h"
#include "setup.h"
#include "logger.h"

//#include "geometry.h"
//#include "memory.h"
//#include "update.h"
//#include "global.h"
//#include "suN_types.h"
//#include "dirac.h"
//#include "linear_algebra.h"
//#include "inverters.h"
//#include "logger.h"
//#include "setup.h"
//#include "hr_complex.h"
//#include "random.h"
//#include "representation.h"

int main(int argc, char *argv[]) 
{
    logger_map("DEBUG", "debug");
    setup_process(&argc, &argv);

    spinor_field *s1;
    s1 = alloc_spinor_field_f(1, &glattice);
    gaussian_spinor_field(s1);
    copy_to_gpu_spinor_field_f(s1);
    double sqnorm = spinor_field_sqnorm_f(s1);

    lprintf("RESULT", 0, "Spinor field sqnorm: %0.2e\n", sqnorm);
    finalize_process();
    return 0;
}