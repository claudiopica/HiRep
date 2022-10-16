/*******************************************************************************
*
* NOCOMPILE= WITH_MPI
*
* Check that after allocating a field, we can write to and read from it.
* This is supposed to be run without MPI as a baseline test
* Identify problems with MPI allocation using the dedicated memory tests.
*
*******************************************************************************/

#define MAIN_PROGRAM

#include "suN.h"
#include "suN_types.h"
#include "setup.h"
#include "global.h"
#include "linear_algebra.h"
#include "basis_linear_algebra.h"
#include "logger.h"
#include "random.h"
#include "memory.h"
#include "update.h"
#include "geometry.h"
#include <math.h>

int test_gfield_allocation();
int test_spinor_field_allocation();

int main(int argc, char *argv[]) 
{
    // Init
    int return_val = 0;

    // Setup process and communication
    logger_map("DEBUG", "debug");
    setup_process(&argc, &argv);

    return_val += test_gfield_allocation();
    return_val += test_spinor_field_allocation();

    // Finalize and return
    finalize_process();
    return return_val;
}

int test_gfield_allocation() 
{
    lprintf("INFO", 0, " ======= TEST GAUGE FIELD ======= \n");
    int return_val = 0;
    suNg_field *f = alloc_gfield(&glattice);
    
    // Fill with random numbers
    random_u(f);

    // Check that sqnorm is unequal to zero, nan or inf
    double sqnorm = sqnorm_gfield_cpu(f);
    if (!isfinite(sqnorm)) {
        lprintf("RESULT", 0, "FAILED\n");
        return_val = 1;
    }
    else 
    {
        lprintf("RESULT", 0, "OK\n");
        return_val = 0;
    }
    lprintf("RESULT", 0, "[Square norm (should be any finite value) %0.2e]\n", sqnorm);

    free_gfield(f);
    return return_val;
}

int test_spinor_field_allocation() 
{
    lprintf("INFO", 0, " ======= TEST SPINOR FIELD ======= \n");
    int return_val = 0;
    spinor_field *f = alloc_spinor_field_f(1, &glattice);

    // Fill with random numbers
    gaussian_spinor_field(f);

    // Check that sqnorm is unequal to zero, nan or inf
    double sqnorm = spinor_field_sqnorm_f_cpu(f);
    if (!isfinite(sqnorm)) {
        lprintf("RESULT", 0, "FAILED\n");
        return_val = 1;
    } 
    else 
    {
        lprintf("RESULT", 0, "OK\n");
        return_val = 0;
    }
    lprintf("RESULT", 0, "[Square norm (should be any finite value) %0.2e]\n", sqnorm);

    free_spinor_field_f(f);
    return return_val;
}