/*******************************************************************************
*
* NOCOMPILE= !WITH_GPU
* NOCOMPILE= !WITH_MPI
*
* Check that copy sync and buffer comms execute without errors
* and don't change the sqnorm of the field
*
*******************************************************************************/

#include "libhr.h"

int test_comms_spinor_field_f();
int test_comms_spinor_field_f_flt();
int test_comms_sfield();
int test_comms_gfield();
int test_comms_gfield_flt();
int test_comms_gfield_f();
int test_comms_gfield_f_flt();
int test_comms_suNg_scalar_field();
int test_comms_avfield();
int test_comms_gtransf();
int test_comms_clover_ldl();
int test_comms_clover_term();
int test_comms_clover_force();

int main(int argc, char *argv[]) 
{
    // Init
    int return_val = 0;

    // Setup process and communication
    logger_map("DEBUG", "debug");
    setup_process(&argc, &argv);
    test_setup();

    // Run tests
      /* Double precision */
    return_val += test_comms_spinor_field_f();
    return_val += test_comms_sfield();
    return_val += test_comms_gfield();
    return_val += test_comms_gfield_f();
    return_val += test_comms_suNg_scalar_field();
    return_val += test_comms_avfield();
    return_val += test_comms_gtransf();
    return_val += test_comms_clover_ldl();
    return_val += test_comms_clover_term();
    return_val += test_comms_clover_force();

      /* Single precision */
    return_val += test_comms_spinor_field_f_flt();
    return_val += test_comms_gfield_flt();
    return_val += test_comms_gfield_f_flt();

    // Finalize and return
    finalize_process();
    return return_val;
}

int test_comms_spinor_field_f() 
{
    lprintf("INFO", 0, " ======= TEST SPINOR FIELD ======= \n");
    lprintf("INFO", 0, "Sqnorm testing\n");

    // Setup fields on GPU
    int return_val = 0;
    spinor_field *f = alloc_spinor_field_f(1, &glattice);
    gaussian_spinor_field(f);
    copy_to_gpu_spinor_field_f(f);

    // Evaluate sqnorm in the beginning
    double sqnorm_start = spinor_field_sqnorm_f(f);
    lprintf("SANITY CHECK", 0, "[In field GPU copy norm unequal zero: %0.2e]\n", sqnorm_start);
    
    // Execute communications
    start_sendrecv_gpu_spinor_field_f(f);
    complete_sendrecv_gpu_spinor_field_f(f);

    // Evaluate sqnorm after comms
    double sqnorm_end = spinor_field_sqnorm_f(f);
    lprintf("SANITY CHECK", 0, "[Out field CPU copy norm unequal zero: %0.2e]\n", sqnorm_end);

    return_val += check_finiteness(sqnorm_start);
    return_val += check_finiteness(sqnorm_end);
    return_val += check_diff_norm_zero(sqnorm_start-sqnorm_end);// TODO: Check diff not diff norm (SAM)

    free_spinor_field_f(f);

    return return_val;
}

int test_comms_spinor_field_f_flt() 
{
    lprintf("INFO", 0, " ======= TEST SPINOR FIELD SINGLE PRECISION ======= \n");

    // Setup fields on GPU
    int return_val = 0;
    spinor_field_flt *f = alloc_spinor_field_f_flt(1, &glattice);
    gaussian_spinor_field_flt(f);
    copy_to_gpu_spinor_field_f_flt(f);

    // Evaluate sqnorm in the beginning
    double sqnorm_start = spinor_field_sqnorm_f_flt(f);
    lprintf("SANITY CHECK", 0, "[In field GPU copy norm unequal zero: %0.2e]\n", sqnorm_start);
    
    // Execute communications
    start_sendrecv_gpu_spinor_field_f_flt(f);
    complete_sendrecv_gpu_spinor_field_f_flt(f);

    // Evaluate sqnorm after comms
    double sqnorm_end = spinor_field_sqnorm_f_flt(f);
    lprintf("SANITY CHECK", 0, "[Out field GPU copy norm unequal zero: %0.2e]\n", sqnorm_end);

    return_val += check_finiteness(sqnorm_start);
    return_val += check_finiteness(sqnorm_end);
    return_val += check_diff_norm_zero(sqnorm_start-sqnorm_end);// TODO: Check diff not diff norm (SAM)

    free_spinor_field_f_flt(f);

    return return_val;
}

int test_comms_sfield() 
{
    lprintf("INFO", 0, " ======= TEST SFIELD ======= \n");

    // Setup fields on GPU
    int return_val = 0;
    scalar_field *f = alloc_sfield(1, &glattice);
    random_sfield_cpu(f);

    // Evaluate sqnorm in the beginning
    double sqnorm_start = sqnorm_sfield_cpu(f);
    lprintf("SANITY CHECK", 0, "[In field CPU copy norm unequal zero: %0.2e]\n", sqnorm_start);
    copy_to_gpu_sfield(f);
    
    // Execute communications
    start_sendrecv_gpu_sfield(f);
    complete_sendrecv_gpu_sfield(f);

    // Evaluate sqnorm after comms
    copy_from_gpu_sfield(f);
    double sqnorm_end = sqnorm_sfield_cpu(f);
    lprintf("SANITY CHECK", 0, "[Out field CPU copy norm unequal zero: %0.2e]\n", sqnorm_end);

    return_val += check_finiteness(sqnorm_start);
    return_val += check_finiteness(sqnorm_end);
    return_val += check_diff_norm_zero(sqnorm_start-sqnorm_end);

    free_sfield(f);

    return return_val;
}

int test_comms_gfield() 
{
    lprintf("INFO", 0, " ======= TEST GFIELD ======= \n");

    // Setup fields on GPU
    int return_val = 0;
    suNg_field *f = alloc_gfield(&glattice);
    random_gfield_cpu(f);

    // Evaluate sqnorm in the beginning
    double sqnorm_start = sqnorm_gfield_cpu(f);
    lprintf("SANITY CHECK", 0, "[In field CPU copy norm unequal zero: %0.2e]\n", sqnorm_start);
    copy_to_gpu_gfield(f);
    
    // Execute communications
    start_sendrecv_gpu_gfield(f);
    complete_sendrecv_gpu_gfield(f);

    // Evaluate sqnorm after comms
    copy_from_gpu_gfield(f);
    double sqnorm_end = sqnorm_gfield_cpu(f);
    lprintf("SANITY CHECK", 0, "[Out field CPU copy norm unequal zero: %0.2e]\n", sqnorm_end);

    return_val += check_finiteness(sqnorm_start);
    return_val += check_finiteness(sqnorm_end);
    return_val += check_diff_norm_zero(sqnorm_start-sqnorm_end);

    free_gfield(f);

    return return_val;
}

int test_comms_gfield_flt() 
{
    lprintf("INFO", 0, " ======= TEST GFIELD SINGLE PRECISION ======= \n");

    // Setup fields on GPU
    int return_val = 0;
    suNg_field_flt *f = alloc_gfield_flt(&glattice);
    random_gfield_flt_cpu(f);

    // Evaluate sqnorm in the beginning
    double sqnorm_start = sqnorm_gfield_flt_cpu(f);
    lprintf("SANITY CHECK", 0, "[In field CPU copy norm unequal zero: %0.2e]\n", sqnorm_start);
    copy_to_gpu_gfield_flt(f);
    
    // Execute communications
    start_sendrecv_gpu_gfield_flt(f);
    complete_sendrecv_gpu_gfield_flt(f);

    // Evaluate sqnorm after comms
    copy_from_gpu_gfield_flt(f);
    double sqnorm_end = sqnorm_gfield_flt_cpu(f);
    lprintf("SANITY CHECK", 0, "[Out field CPU copy norm unequal zero: %0.2e]\n", sqnorm_end);

    return_val += check_finiteness(sqnorm_start);
    return_val += check_finiteness(sqnorm_end);
    return_val += check_diff_norm_zero(sqnorm_start-sqnorm_end);

    free_gfield_flt(f);

    return return_val;
}

int test_comms_gfield_f() 
{
    lprintf("INFO", 0, " ======= TEST GFIELD FERMION REP ======= \n");

    // Setup fields on GPU
    int return_val = 0;
    suNf_field *f = alloc_gfield_f(&glattice);
    random_gfield_f_cpu(f);

    // Evaluate sqnorm in the beginning
    double sqnorm_start = sqnorm_gfield_f_cpu(f);
    lprintf("SANITY CHECK", 0, "[In field CPU copy norm unequal zero: %0.2e]\n", sqnorm_start);
    copy_to_gpu_gfield_f(f);
    
    // Execute communications
    start_sendrecv_gpu_gfield_f(f);
    complete_sendrecv_gpu_gfield_f(f);

    // Evaluate sqnorm after comms
    copy_from_gpu_gfield_f(f);
    double sqnorm_end = sqnorm_gfield_f_cpu(f);
    lprintf("SANITY CHECK", 0, "[Out field CPU copy norm unequal zero: %0.2e]\n", sqnorm_end);

    return_val += check_finiteness(sqnorm_start);
    return_val += check_finiteness(sqnorm_end);
    return_val += check_diff_norm_zero(sqnorm_start-sqnorm_end);

    free_gfield_f(f);

    return return_val;
}

int test_comms_gfield_f_flt() 
{
    lprintf("INFO", 0, " ======= TEST GFIELD FERMION REP SINGLE PRECISION ======= \n");

    // Setup fields on GPU
    int return_val = 0;
    suNf_field_flt *f = alloc_gfield_f_flt(&glattice);
    random_gfield_f_flt_cpu(f);

    // Evaluate sqnorm in the beginning
    double sqnorm_start = sqnorm_gfield_f_flt_cpu(f);
    lprintf("SANITY CHECK", 0, "[In field CPU copy norm unequal zero: %0.2e]\n", sqnorm_start);
    copy_to_gpu_gfield_f_flt(f);
    
    // Execute communications
    start_sendrecv_gpu_gfield_f_flt(f);
    complete_sendrecv_gpu_gfield_f_flt(f);

    // Evaluate sqnorm after comms
    copy_from_gpu_gfield_f_flt(f);
    double sqnorm_end = sqnorm_gfield_f_flt_cpu(f);
    lprintf("SANITY CHECK", 0, "[Out field CPU copy norm unequal zero: %0.2e]\n", sqnorm_end);

    return_val += check_finiteness(sqnorm_start);
    return_val += check_finiteness(sqnorm_end);
    return_val += check_diff_norm_zero(sqnorm_start-sqnorm_end);

    free_gfield_f_flt(f);

    return return_val;
}

int test_comms_suNg_scalar_field() 
{
    lprintf("INFO", 0, " ======= TEST SU(NG) SCALAR FIELD ======= \n");

    // Setup fields on GPU
    int return_val = 0;
    suNg_scalar_field *f = alloc_suNg_scalar_field(&glattice);
    random_suNg_scalar_field_cpu(f);

    // Evaluate sqnorm in the beginning
    double sqnorm_start = sqnorm_suNg_scalar_field_cpu(f);
    lprintf("SANITY CHECK", 0, "[In field CPU copy norm unequal zero: %0.2e]\n", sqnorm_start);
    copy_to_gpu_suNg_scalar_field(f);
    
    // Execute communications
    start_sendrecv_gpu_suNg_scalar_field(f);
    complete_sendrecv_gpu_suNg_scalar_field(f);

    // Evaluate sqnorm after comms
    copy_from_gpu_suNg_scalar_field(f);
    double sqnorm_end = sqnorm_suNg_scalar_field_cpu(f);
    lprintf("SANITY CHECK", 0, "[Out field CPU copy norm unequal zero: %0.2e]\n", sqnorm_end);

    return_val += check_finiteness(sqnorm_start);
    return_val += check_finiteness(sqnorm_end);
    return_val += check_diff_norm_zero(sqnorm_start-sqnorm_end);

    free_suNg_scalar_field(f);

    return return_val;
}

int test_comms_avfield() 
{
    lprintf("INFO", 0, " ======= TEST ALGEBRA VECTOR FIELD ======= \n");

    // Setup fields on GPU
    int return_val = 0;
    suNg_av_field *f = alloc_avfield(&glattice);
    random_avfield_cpu(f);

    // Evaluate sqnorm in the beginning
    double sqnorm_start = sqnorm_avfield_cpu(f);
    lprintf("SANITY CHECK", 0, "[In field CPU copy norm unequal zero: %0.2e]\n", sqnorm_start);
    copy_to_gpu_avfield(f);
    
    // Execute communications
    start_sendrecv_gpu_avfield(f);
    complete_sendrecv_gpu_avfield(f);

    // Evaluate sqnorm after comms
    copy_from_gpu_avfield(f);
    double sqnorm_end = sqnorm_avfield_cpu(f);
    lprintf("SANITY CHECK", 0, "[Out field CPU copy norm unequal zero: %0.2e]\n", sqnorm_end);

    return_val += check_finiteness(sqnorm_start);
    return_val += check_finiteness(sqnorm_end);
    return_val += check_diff_norm_zero(sqnorm_start-sqnorm_end);

    free_avfield(f);

    return return_val;
}

int test_comms_gtransf() 
{
    lprintf("INFO", 0, " ======= TEST GAUGE TRANSFORMATION ======= \n");

    // Setup fields on GPU
    int return_val = 0;
    suNg_field *f = alloc_gtransf(&glattice);
    random_gtransf_cpu(f);

    // Evaluate sqnorm in the beginning
    double sqnorm_start = sqnorm_gtransf_cpu(f);
    lprintf("SANITY CHECK", 0, "[In field CPU copy norm unequal zero: %0.2e]\n", sqnorm_start);
    copy_to_gpu_gtransf(f);
    
    // Execute communications
    start_sendrecv_gpu_gtransf(f);
    complete_sendrecv_gpu_gtransf(f);

    // Evaluate sqnorm after comms
    copy_from_gpu_gtransf(f);
    double sqnorm_end = sqnorm_gtransf_cpu(f);
    lprintf("SANITY CHECK", 0, "[Out field CPU copy norm unequal zero: %0.2e]\n", sqnorm_end);

    return_val += check_finiteness(sqnorm_start);
    return_val += check_finiteness(sqnorm_end);
    return_val += check_diff_norm_zero(sqnorm_start-sqnorm_end);

    free_gtransf(f);

    return return_val;
}

int test_comms_clover_ldl() 
{
    lprintf("INFO", 0, " ======= TEST CLOVER LDL ======= \n");

    // Setup fields on GPU
    int return_val = 0;
    ldl_field *f = alloc_clover_ldl(&glattice);
    random_clover_ldl_cpu(f);

    // Evaluate sqnorm in the beginning
    double sqnorm_start = sqnorm_clover_ldl_cpu(f);
    lprintf("SANITY CHECK", 0, "[In field CPU copy norm unequal zero: %0.2e]\n", sqnorm_start);
    copy_to_gpu_clover_ldl(f);
    
    // Execute communications
    start_sendrecv_gpu_clover_ldl(f);
    complete_sendrecv_gpu_clover_ldl(f);

    // Evaluate sqnorm after comms
    copy_from_gpu_clover_ldl(f);
    double sqnorm_end = sqnorm_clover_ldl_cpu(f);
    lprintf("SANITY CHECK", 0, "[Out field CPU copy norm unequal zero: %0.2e]\n", sqnorm_end);

    return_val += check_finiteness(sqnorm_start);
    return_val += check_finiteness(sqnorm_end);
    return_val += check_diff_norm_zero(sqnorm_start-sqnorm_end);

    free_clover_ldl(f);

    return return_val;
}

int test_comms_clover_term() 
{
    lprintf("INFO", 0, " ======= TEST CLOVER TERM ======= \n");

    // Setup fields on GPU
    int return_val = 0;
    suNfc_field *f = alloc_clover_term(&glattice);
    random_clover_term_cpu(f);

    // Evaluate sqnorm in the beginning
    double sqnorm_start = sqnorm_clover_term_cpu(f);
    lprintf("SANITY CHECK", 0, "[In field CPU copy norm unequal zero: %0.2e]\n", sqnorm_start);
    copy_to_gpu_clover_term(f);
    
    // Execute communications
    start_sendrecv_gpu_clover_term(f);
    complete_sendrecv_gpu_clover_term(f);

    // Evaluate sqnorm after comms
    copy_from_gpu_clover_term(f);
    double sqnorm_end = sqnorm_clover_term_cpu(f);
    lprintf("SANITY CHECK", 0, "[Out field CPU copy norm unequal zero: %0.2e]\n", sqnorm_end);

    return_val += check_finiteness(sqnorm_start);
    return_val += check_finiteness(sqnorm_end);
    return_val += check_diff_norm_zero(sqnorm_start-sqnorm_end);

    free_clover_term(f);

    return return_val;
}

int test_comms_clover_force() 
{
    lprintf("INFO", 0, " ======= TEST CLOVER FORCE ======= \n");

    // Setup fields on GPU
    int return_val = 0;
    suNf_field *f = alloc_clover_force(&glattice);
    random_clover_force_cpu(f);

    // Evaluate sqnorm in the beginning
    double sqnorm_start = sqnorm_clover_force_cpu(f);
    lprintf("SANITY CHECK", 0, "[In field CPU copy norm unequal zero: %0.2e]\n", sqnorm_start);
    copy_to_gpu_clover_force(f);
    
    // Execute communications
    start_sendrecv_gpu_clover_force(f);
    complete_sendrecv_gpu_clover_force(f);

    // Evaluate sqnorm after comms
    copy_from_gpu_clover_force(f);
    double sqnorm_end = sqnorm_clover_force_cpu(f);
    lprintf("SANITY CHECK", 0, "[Out field CPU copy norm unequal zero: %0.2e]\n", sqnorm_end);

    return_val += check_finiteness(sqnorm_start);
    return_val += check_finiteness(sqnorm_end);
    return_val += check_diff_norm_zero(sqnorm_start-sqnorm_end);

    free_clover_force(f);

    return return_val;
}

#endif