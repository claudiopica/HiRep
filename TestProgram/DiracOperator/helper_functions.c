#include <stdbool.h>
#include "logger.h"
#include "geometry.h"
#include "setup.h"
#include "global.h"
#include "random.h"
#include "representation.h"
#include "suN_types.h"
#include "memory.h"
#include "update.h"

/*
 * Prints log ouput for failing or passing of test.
 *
 * in: bool test_pass      i.e. a test function that returns a bool value
 *     int global_pass     variable defined in this scope that eventually contains information 
 *                         about whether all tests in this suite have passed or whether there 
 *                         was a fail somewhere. If the test function fails then this value will
 *                         be set to fail.
 * */
void run_test(bool test_pass, int global_pass)
{
    if (test_pass) lprintf("RESULT", 0, "OK");
    else
    {
        global_pass = 1;
        lprintf("RESULT", 0, "FAILED");
    }
}

/*
 * Initializes lattice so that one can test spinor operations on it.
 *
 * in: int argc, char *argv[] are command line options from main
 *
 * */
void init_test(int argc, char *argv[]) 
{
    setup_process(&argc, &argv);
    setup_gauge_fields();
    random_u(u_gauge);
    represent_gauge_field();
}

/*
 * Allocates spinor field, fills it with random values and copies to gpu.
 *
 * out: spinor_field pointer of prepared random spinor field on GPU.
 *
 * */
spinor_field* setup_infield() 
{
    spinor_field *s;
    s = alloc_spinor_field_f(1, &glattice);
    gaussian_spinor_field(s);
    spinor_field_copy_to_gpu_f(s);
    return s;
}
