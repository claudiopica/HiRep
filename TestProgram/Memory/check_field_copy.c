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

int test_bijectivity_gfield();
int test_bijectivity_gfield_f();
int test_bijectivity_spinor_field();

int main(int argc, char *argv[]) 
{
    // Init
    int return_val = 0;

    // Setup process and communication
    logger_map("DEBUG", "debug");
    setup_process(&argc, &argv);

    // Test block
    return_val += test_bijectivity_gfield();
    return_val += test_bijectivity_gfield_f();
    return_val += test_bijectivity_spinor_field();

    // Finalize and return
    finalize_process();
    return return_val;
}

int test_bijectivity_gfield() 
{
    lprintf("INFO", 0, " ====== TEST GAUGE FIELD ======= ");
    int return_val = 0;
    suNg_field *in, *in_copy;
    in = alloc_gfield(&glattice);
    in_copy = alloc_gfield(&glattice);

    random_gfield_cpu(in);
    //random_gfield_cpu(in_copy);

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
    lprintf("RESULT", 0, "[Diff norm %0.2e]\n", diff_norm);

    free_gfield(in);
    free_gfield(in_copy);
    return return_val;
}

int test_bijectivity_gfield_f() 
{
    lprintf("INFO", 0, " ====== TEST GAUGE FIELD ======= ");
    int return_val = 0;
    suNf_field *in, *in_copy;
    in = alloc_gfield_f(&glattice);
    in_copy = alloc_gfield_f(&glattice);

    random_gfield_f_cpu(in);
    random_gfield_f_cpu(in_copy);

    copy_gfield_f_cpu(in_copy, in);
    lprintf("SANITY CHECK", 0, "CPU sqnorm: %0.2e\n", sqnorm_gfield_f_cpu(in));
    lprintf("SANITY CHECK", 0, "CPU copy sqnorm (should be the same as CPU sqnorm): %0.2e\n", sqnorm_gfield_f_cpu(in_copy));

    copy_to_gpu_gfield_f(in);

    //zero_gfield_f_cpu(in);
    lprintf("SANITY CHECK", 0, "CPU copy should be zero in intermediate step: %0.2e\n", sqnorm_gfield_f_cpu(in));
    copy_from_gpu_gfield_f(in);

    sub_assign_gfield_f_cpu(in, in_copy);
    double diff_norm = sqnorm_gfield_f_cpu(in);

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
    lprintf("RESULT", 0, "[Diff norm %0.2e]\n", diff_norm);

    free_gfield_f(in);
    free_gfield_f(in_copy);
    return return_val;
}

int test_bijectivity_spinor_field() 
{
    lprintf("INFO", 0, " ====== TEST SPINOR FIELD ======= ");
    int return_val = 0;
    spinor_field *in, *in_copy;
    in = alloc_spinor_field_f(1, &glattice);
    in_copy = alloc_spinor_field_f(1, &glattice);

    gaussian_spinor_field(in);

    spinor_field_copy_f_cpu(in_copy, in);
    lprintf("SANITY CHECK", 0, "CPU sqnorm: %0.2e\n", spinor_field_sqnorm_f_cpu(in));
    lprintf("SANITY CHECK", 0, "CPU copy sqnorm (should be the same as CPU sqnorm): %0.2e\n", spinor_field_sqnorm_f_cpu(in_copy));

    copy_to_gpu_spinor_field_f(in);

    spinor_field_zero_f_cpu(in);
    lprintf("SANITY CHECK", 0, "CPU copy should be zero in intermediate step: %0.2e\n", spinor_field_sqnorm_f_cpu(in));
    lprintf("SANITY CHECK", 0, "GPU copy should be equal to ealier in square norms: %0.2e\n", spinor_field_sqnorm_f(in));
    copy_from_gpu_spinor_field_f(in);

    spinor_field_sub_assign_f_cpu(in, in_copy);
    double diff_norm = spinor_field_sqnorm_f_cpu(in);

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
    lprintf("RESULT", 0, "[Diff norm %0.2e]\n", diff_norm);

    free_spinor_field_f(in);
    free_spinor_field_f(in_copy);
    return return_val;
}


