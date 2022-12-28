/*******************************************************************************
*
* NOCOMPILE= !WITH_GPU
*
* Check that after allocating a field, we can write to and read from it.
* This is supposed to be run without MPI as a baseline test
* Identify problems with MPI allocation using the dedicated memory tests.
*
*******************************************************************************/

#include "libhr.h"

// TODO: This test relies on too many high level things.
//       I am not sure, we can run such a test at all.
// TODO: Test for GPU whether CPU allocation works, too. The CPU test is not enough because the function names are different
// TODO: Also do not copy back, but check with GPU norm.

// Double precision
int test_gfield_allocation();
int test_gfield_f_allocation();
int test_spinor_field_allocation();
int test_avfield_allocation();
int test_clover_term_allocation();
int test_clover_force_allocation();
int test_sfield_allocation();

// Single precision
int test_gfield_flt_allocation();
int test_gfield_f_flt_allocation();
int test_spinor_field_flt_allocation();

int main(int argc, char *argv[]) 
{
    // Init
    int return_val = 0;

    // Setup process and communication
    logger_map("DEBUG", "debug");
    setup_process(&argc, &argv);
    test_setup();

    // Double precision test block
    return_val += test_gfield_allocation();
    return_val += test_gfield_f_allocation();
    return_val += test_spinor_field_allocation();
    //return_val += test_avfield_allocation();
    //return_val += test_clover_term_allocation();
    //return_val += test_clover_force_allocation();
    //return_val += test_sfield_allocation(); // FIXME: Bus error

    // Single precision test block
    //return_val += test_gfield_flt_allocation();
    //return_val += test_gfield_f_flt_allocation();

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
    random_gfield_cpu(f);

    // Copy back and forth
    copy_to_gpu_gfield(f);
    copy_from_gpu_gfield(f);

    // Check that sqnorm is unequal to zero, nan or inf
    double sqnorm = sqnorm_gfield_cpu(f);
    if (!isfinite(sqnorm) || fabs(sqnorm)  < 1e-14) {
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

int test_gfield_f_allocation() 
{
    lprintf("INFO", 0, " ======= TEST REPRESENTED GAUGE FIELD ======= \n");
    int return_val = 0;
    suNf_field *f = alloc_gfield_f(&glattice);

    // Fill with random numbers
    random_gfield_f_cpu(f);

    // Copy back and forth
    copy_to_gpu_gfield_f(f);
    copy_from_gpu_gfield_f(f);

    // Check that sqnorm is unequal to zero, nan or inf
    double sqnorm = sqnorm_gfield_f_cpu(f);
    if (!isfinite(sqnorm) || fabs(sqnorm) < 1e-14) {
        lprintf("RESULT", 0, "FAILED\n");
        return_val = 1;
    }
    else 
    {
        lprintf("RESULT", 0, "OK\n");
        return_val = 0;
    }
    lprintf("RESULT", 0, "[Square norm (should be any finite value) %0.2e]\n", sqnorm);

    free_gfield_f(f);
    return return_val;
}

int test_gfield_flt_allocation() 
{
    lprintf("INFO", 0, " ======= TEST GAUGE FIELD SINGLE PRECISION ======= \n");
    int return_val = 0;
    suNg_field_flt *f = alloc_gfield_flt(&glattice);
    
    // Fill with random numbers
    random_gfield_flt_cpu(f);

    // Copy back and forth
    copy_to_gpu_gfield_flt(f);
    copy_from_gpu_gfield_flt(f);

    // Check that sqnorm is unequal to zero, nan or inf
    double sqnorm = sqnorm_gfield_flt_cpu(f);
    if (!isfinite(sqnorm) || fabs(sqnorm) < 1e-14) {
        lprintf("RESULT", 0, "FAILED\n");
        return_val = 1;
    }
    else 
    {
        lprintf("RESULT", 0, "OK\n");
        return_val = 0;
    }
    lprintf("RESULT", 0, "[Square norm (should be any finite value) %0.2e]\n", sqnorm);

    free_gfield_flt(f);
    return return_val;
}

int test_gfield_f_flt_allocation() 
{
    lprintf("INFO", 0, " ======= TEST REPRESENTED GAUGE FIELD SINGLE PRECISION ======= \n");
    int return_val = 0;
    suNf_field_flt *f = alloc_gfield_f_flt(&glattice);
    
    // Fill with random numbers
    random_gfield_f_flt_cpu(f);

    // Copy back and forth
    copy_to_gpu_gfield_f_flt(f);
    copy_from_gpu_gfield_f_flt(f);

    // Check that sqnorm is unequal to zero, nan or inf
    double sqnorm = sqnorm_gfield_f_flt_cpu(f);
    if (!isfinite(sqnorm) || fabs(sqnorm) < 1e-14) {
        lprintf("RESULT", 0, "FAILED\n");
        return_val = 1;
    }
    else 
    {
        lprintf("RESULT", 0, "OK\n");
        return_val = 0;
    }
    lprintf("RESULT", 0, "[Square norm (should be any finite value) %0.2e]\n", sqnorm);

    free_gfield_f_flt(f);
    return return_val;
}

int test_avfield_allocation() 
{
    lprintf("INFO", 0, " ======= TEST SU(N) ALGEBRA VECTOR FIELD ======= \n");
    int return_val = 0;
    suNg_av_field *f = alloc_avfield(&glattice);
    
    // Fill with random numbers
    random_avfield_cpu(f);

    // Copy back and forth
    copy_to_gpu_avfield(f);
    copy_from_gpu_avfield(f);

    // Check that sqnorm is unequal to zero, nan or inf
    double sqnorm = sqnorm_avfield_cpu(f);
    if (!isfinite(sqnorm) || fabs(sqnorm) < 1e-14) {
        lprintf("RESULT", 0, "FAILED\n");
        return_val = 1;
    }
    else 
    {
        lprintf("RESULT", 0, "OK\n");
        return_val = 0;
    }
    lprintf("RESULT", 0, "[Square norm (should be any finite value) %0.2e]\n", sqnorm);

    free_avfield(f);
    return return_val;
}

int test_gtransf_allocation() 
{
    lprintf("INFO", 0, " ======= GAUGE TRANSFORMATION ======= \n");
    int return_val = 0;
    suNg_field *f = alloc_gtransf(&glattice);
    
    // Fill with random numbers
    random_gtransf_cpu(f);

    // Copy back and forth
    copy_to_gpu_gtransf(f);
    copy_from_gpu_gtransf(f);

    // Check that sqnorm is unequal to zero, nan or inf
    double sqnorm = sqnorm_gtransf_cpu(f);
    if (!isfinite(sqnorm) || fabs(sqnorm) < 1e-14) {
        lprintf("RESULT", 0, "FAILED\n");
        return_val = 1;
    }
    else 
    {
        lprintf("RESULT", 0, "OK\n");
        return_val = 0;
    }
    lprintf("RESULT", 0, "[Square norm (should be any finite value) %0.2e]\n", sqnorm);

    free_gtransf(f);
    return return_val;
}

int test_clover_term_allocation() 
{
    lprintf("INFO", 0, " ======= CLOVER TERM ======= \n");
    int return_val = 0;
    suNfc_field *f = alloc_clover_term(&glattice);
    
    // Fill with random numbers
    random_clover_term_cpu(f);

    // Copy back and forth
    copy_to_gpu_clover_term(f);
    copy_from_gpu_clover_term(f);

    // Check that sqnorm is unequal to zero, nan or inf
    double sqnorm = sqnorm_clover_term_cpu(f);
    if (!isfinite(sqnorm) || fabs(sqnorm) < 1e-14) {
        lprintf("RESULT", 0, "FAILED\n");
        return_val = 1;
    }
    else 
    {
        lprintf("RESULT", 0, "OK\n");
        return_val = 0;
    }
    lprintf("RESULT", 0, "[Square norm (should be any finite value) %0.2e]\n", sqnorm);

    free_clover_term(f);
    return return_val;
}

int test_clover_force_allocation()
{
    lprintf("INFO", 0, " ======= CLOVER FORCE ======= \n");
    int return_val = 0;
    suNf_field *f = alloc_clover_force(&glattice);
    
    // Fill with random numbers
    random_clover_force_cpu(f);

    // Copy back and forth
    copy_to_gpu_clover_force(f);
    copy_from_gpu_clover_force(f);

    // Check that sqnorm is unequal to zero, nan or inf
    double sqnorm = sqnorm_clover_force_cpu(f);
    if (!isfinite(sqnorm) || fabs(sqnorm) < 1e-14) {
        lprintf("RESULT", 0, "FAILED\n");
        return_val = 1;
    }
    else 
    {
        lprintf("RESULT", 0, "OK\n");
        return_val = 0;
    }
    lprintf("RESULT", 0, "[Square norm (should be any finite value) %0.2e]\n", sqnorm);

    free_clover_force(f);
    return return_val;
}

int test_spinor_field_allocation() 
{
    lprintf("INFO", 0, " ======= TEST SPINOR FIELD ======= \n");
    int return_val = 0;
    spinor_field *f = alloc_spinor_field_f(1, &glattice);

    // Fill with random numbers
    gaussian_spinor_field(f);

    // Copy back and forth
    /*copy_to_gpu_spinor_field_f(f);
    copy_from_gpu_spinor_field_f(f);

    // Check that sqnorm is unequal to zero, nan or inf
    double sqnorm = spinor_field_sqnorm_f_cpu(f);
    if (!isfinite(sqnorm) || fabs(sqnorm) < 1e-14) {
        lprintf("RESULT", 0, "FAILED\n");
        return_val = 1;
    } 
    else 
    {
        lprintf("RESULT", 0, "OK\n");
        return_val = 0;
    }
    lprintf("RESULT", 0, "[Square norm (should be any finite value) %0.2e]\n", sqnorm);

    free_spinor_field_f(f); */
    return return_val;
}

int test_spinor_field_flt_allocation() 
{
    lprintf("INFO", 0, " ======= TEST SPINOR FIELD ======= \n");
    int return_val = 0;
    spinor_field_flt *f = alloc_spinor_field_f_flt(1, &glattice);

    // Fill with random numbers
    gaussian_spinor_field_flt(f);

    // Copy back and forth
    copy_to_gpu_spinor_field_f_flt(f);
    copy_from_gpu_spinor_field_f_flt(f);

    // Check that sqnorm is unequal to zero, nan or inf
    double sqnorm = spinor_field_sqnorm_f_flt_cpu(f);
    if (!isfinite(sqnorm) || fabs(sqnorm) < 1e-14) {
        lprintf("RESULT", 0, "FAILED\n");
        return_val = 1;
    } 
    else 
    {
        lprintf("RESULT", 0, "OK\n");
        return_val = 0;
    }
    lprintf("RESULT", 0, "[Square norm (should be any finite value) %0.2e]\n", sqnorm);

    free_spinor_field_f_flt(f);
    return return_val;
}

int test_sfield_allocation() 
{
    lprintf("INFO", 0, " ======= TEST SPINOR FIELD ======= \n");
    int return_val = 0;
    scalar_field *f = alloc_sfield(1, &glattice);

    // Fill with random numbers
    random_sfield_cpu(f);

    // Copy back and forth
    copy_to_gpu_sfield(f);
    copy_from_gpu_sfield(f);

    // Check that sqnorm is unequal to zero, nan or inf
    double sqnorm = sqnorm_sfield_cpu(f);
    if (!isfinite(sqnorm) || fabs(sqnorm) < 1e-14) {
        lprintf("RESULT", 0, "FAILED\n");
        return_val = 1;
    } 
    else 
    {
        lprintf("RESULT", 0, "OK\n");
        return_val = 0;
    }
    lprintf("RESULT", 0, "[Square norm (should be any finite value) %0.2e]\n", sqnorm);

    free_sfield(f);
    return return_val;
}




