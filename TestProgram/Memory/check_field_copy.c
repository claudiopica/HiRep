/*******************************************************************************
*
* NOCOMPILE= !WITH_GPU
*
* Check that copy back and forth works for all fields.
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

int test_gfield_bijectivity();

int main(int argc, char *argv[]) 
{
    // Init
    int return_val = 0;

    // Setup process and communication
    logger_map("DEBUG", "debug");
    setup_process(&argc, &argv);

    lprintf("Testing copy_to_gpu and copy_from_gpu functions for all field types.\n");

    return_val += test_gfield_bijectivity();

    // Finalize and return
    finalize_process();
    return return_val;
}

int test_gfield_bijectivity() 
{
    lprintf("INFO", 0, " ====== TEST GAUGE FIELD ======= ");
    suNg_field *in, *in_copy;
    in = alloc_gfield(&glattice);
    in_copy = alloc_gfield(&glattice);

    random_u(in);

    copy_gfield_cpu(in_copy, in);
    lprintf("SANITY CHECK", 0, "CPU sqnorm: %0.2e\n", sqnorm_gfield_cpu(in));
    lprintf("SANITY CHECK", 0, "CPU copy sqnorm (should be the same as CPU sqnorm): %0.2e\n", sqnorm_gfield_cpu(in_copy));

    copy_to_gpu_gfield(in);

    zero_gfield_cpu(in);
    lprintf("SANITY CHECK", 0, "CPU copy should be zero in intermediate step: %0.2e\n", sqnorm_gfield_cpu(in));
    copy_from_gpu_gfield(in);

    sub_assign_gfield_cpu(in, in_copy);
    double diff_norm = sqnorm_gfield_cpu(in);

    if (diff_norm != 0) 
    {
        lprintf("RESULT", 0, "FAILED\n");
        return_val = 1;
    } 
    else 
    {
        lprintf("RESULT", 0, "OK \n");
        return_val = 0;
    }
    lprintf("RESULT", 0, "[Diff norm %0.2]\n", sqnorm);

    free_gfield(in);
    free_gfield(in_copy);
}