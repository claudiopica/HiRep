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
#include <string.h>

int test_comms_spinor_field(geometry_descriptor *);
int test_comms_spinor_field_flt(geometry_descriptor *);
int test_comms_scalar_field(geometry_descriptor *);
int test_comms_suNg_field();
int test_comms_suNg_field_flt();
int test_comms_suNf_field();
int test_comms_suNf_field_flt();
int test_comms_suNg_scalar_field();
int test_comms_suNg_av_field();
int test_comms_gtransf();
int test_comms_ldl_field();
int test_comms_clover_term();
int test_comms_clover_force();
int test_comms_staple_field();

int main(int argc, char *argv[]) {
    // Init
    int return_val = 0;

    // Setup process and communication
    logger_map("DEBUG", "debug");
    setup_process(&argc, &argv);
    test_setup();

    // Run tests
    /* Double precision */
    lprintf("INFO", 0, "\n\nFull lattice tests\n\n");
    return_val += test_comms_spinor_field(&glattice);
    return_val += test_comms_scalar_field(&glattice);
    return_val += test_comms_suNg_field();
    return_val += test_comms_suNf_field();
    return_val += test_comms_suNg_scalar_field();
    return_val += test_comms_suNg_av_field();
    return_val += test_comms_gtransf();
    return_val += test_comms_ldl_field();
    return_val += test_comms_clover_term();
    return_val += test_comms_clover_force();
    return_val += test_comms_staple_field();

    /* Single precision */
    return_val += test_comms_spinor_field_flt(&glattice);
    return_val += test_comms_suNg_field_flt();
    return_val += test_comms_suNf_field_flt();

    lprintf("INFO", 0, "\n\nSpinor tests on even lattice\n\n");
    return_val += test_comms_spinor_field(&glat_even);
    return_val += test_comms_spinor_field_flt(&glat_even);
    return_val += test_comms_scalar_field(&glat_even);

    lprintf("INFO", 0, "\n\nSpinor tests on odd lattice\n\n");
    return_val += test_comms_spinor_field(&glat_odd);
    return_val += test_comms_spinor_field_flt(&glat_odd);
    return_val += test_comms_scalar_field(&glat_odd);

    // Finalize and return
    finalize_process();
    return return_val;
}

int test_comms_spinor_field(geometry_descriptor *gd) {
    lprintf("INFO", 0, " ======= TEST SPINOR FIELD ======= \n");
    lprintf("INFO", 0, "Sqnorm testing\n");

    // Setup fields on GPU
    int return_val = 0;
    spinor_field *f = alloc_spinor_field(1, gd);
    gaussian_spinor_field(f);
    copy_to_gpu_spinor_field(f);

    // Evaluate sqnorm in the beginning
    double sqnorm_start = spinor_field_sqnorm_f(f);
    lprintf("SANITY CHECK", 0, "[In field GPU copy norm unequal zero: %0.2e]\n", sqnorm_start);

    // Execute communications
    start_sendrecv_spinor_field(f);
    complete_sendrecv_spinor_field(f);

    // Evaluate sqnorm after comms
    double sqnorm_end = spinor_field_sqnorm_f(f);
    lprintf("SANITY CHECK", 0, "[Out field CPU copy norm unequal zero: %0.2e]\n", sqnorm_end);

    return_val += check_finiteness(sqnorm_start);
    return_val += check_finiteness(sqnorm_end);
    return_val += check_diff_norm_zero(sqnorm_start - sqnorm_end); // TODO: Check diff not diff norm (SAM)

    free_spinor_field(f);

    return return_val;
}

int test_comms_spinor_field_flt(geometry_descriptor *gd) {
    lprintf("INFO", 0, " ======= TEST SPINOR FIELD SINGLE PRECISION ======= \n");

    // Setup fields on GPU
    int return_val = 0;
    spinor_field_flt *f = alloc_spinor_field_flt(1, gd);
    gaussian_spinor_field_flt(f);
    copy_to_gpu_spinor_field_flt(f);

    // Evaluate sqnorm in the beginning
    double sqnorm_start = spinor_field_sqnorm_f_flt(f);
    lprintf("SANITY CHECK", 0, "[In field GPU copy norm unequal zero: %0.2e]\n", sqnorm_start);

    // Execute communications
    start_sendrecv_spinor_field_flt(f);
    complete_sendrecv_spinor_field_flt(f);

    // Evaluate sqnorm after comms
    double sqnorm_end = spinor_field_sqnorm_f_flt(f);
    lprintf("SANITY CHECK", 0, "[Out field GPU copy norm unequal zero: %0.2e]\n", sqnorm_end);

    return_val += check_finiteness(sqnorm_start);
    return_val += check_finiteness(sqnorm_end);
    return_val += check_diff_norm_zero(sqnorm_start - sqnorm_end); // TODO: Check diff not diff norm (SAM)

    free_spinor_field_flt(f);

    return return_val;
}

int test_comms_scalar_field(geometry_descriptor *gd) {
    lprintf("INFO", 0, " ======= TEST SFIELD ======= \n");

    // Setup fields on GPU
    int return_val = 0;
    scalar_field *f = alloc_scalar_field(1, gd);
    random_scalar_field_cpu(f);

    // Evaluate sqnorm in the beginning
    double sqnorm_start = sqnorm_scalar_field_cpu(f);
    lprintf("SANITY CHECK", 0, "[In field CPU copy norm unequal zero: %0.2e]\n", sqnorm_start);
    copy_to_gpu_scalar_field(f);

    // Execute communications
    start_sendrecv_scalar_field(f);
    complete_sendrecv_scalar_field(f);

    // Evaluate sqnorm after comms
    copy_from_gpu_scalar_field(f);
    double sqnorm_end = sqnorm_scalar_field_cpu(f);
    lprintf("SANITY CHECK", 0, "[Out field CPU copy norm unequal zero: %0.2e]\n", sqnorm_end);

    return_val += check_finiteness(sqnorm_start);
    return_val += check_finiteness(sqnorm_end);
    return_val += check_diff_norm_zero(sqnorm_start - sqnorm_end);

    free_scalar_field(f);

    return return_val;
}

int test_comms_suNg_field() {
    lprintf("INFO", 0, " ======= TEST GFIELD ======= \n");

    // Setup fields on GPU
    int return_val = 0;
    suNg_field *f = alloc_suNg_field(&glattice);
    random_suNg_field_cpu(f);

    // Evaluate sqnorm in the beginning
    double sqnorm_start = sqnorm_suNg_field_cpu(f);
    lprintf("SANITY CHECK", 0, "[In field CPU copy norm unequal zero: %0.2e]\n", sqnorm_start);
    copy_to_gpu_suNg_field(f);

    // Execute communications
    start_sendrecv_suNg_field(f);
    complete_sendrecv_suNg_field(f);

    // Evaluate sqnorm after comms
    copy_from_gpu_suNg_field(f);
    double sqnorm_end = sqnorm_suNg_field_cpu(f);
    lprintf("SANITY CHECK", 0, "[Out field CPU copy norm unequal zero: %0.2e]\n", sqnorm_end);

    return_val += check_finiteness(sqnorm_start);
    return_val += check_finiteness(sqnorm_end);
    return_val += check_diff_norm_zero(sqnorm_start - sqnorm_end);

    free_suNg_field(f);

    return return_val;
}

int test_comms_suNg_field_flt() {
    lprintf("INFO", 0, " ======= TEST GFIELD SINGLE PRECISION ======= \n");

    // Setup fields on GPU
    int return_val = 0;
    suNg_field_flt *f = alloc_suNg_field_flt(&glattice);
    random_suNg_field_flt_cpu(f);

    // Evaluate sqnorm in the beginning
    double sqnorm_start = sqnorm_suNg_field_flt_cpu(f);
    lprintf("SANITY CHECK", 0, "[In field CPU copy norm unequal zero: %0.2e]\n", sqnorm_start);
    copy_to_gpu_suNg_field_flt(f);

    // Execute communications
    start_sendrecv_suNg_field_flt(f);
    complete_sendrecv_suNg_field_flt(f);

    // Evaluate sqnorm after comms
    copy_from_gpu_suNg_field_flt(f);
    double sqnorm_end = sqnorm_suNg_field_flt_cpu(f);
    lprintf("SANITY CHECK", 0, "[Out field CPU copy norm unequal zero: %0.2e]\n", sqnorm_end);

    return_val += check_finiteness(sqnorm_start);
    return_val += check_finiteness(sqnorm_end);
    return_val += check_diff_norm_zero(sqnorm_start - sqnorm_end);

    free_suNg_field_flt(f);

    return return_val;
}

int test_comms_suNf_field() {
    lprintf("INFO", 0, " ======= TEST GFIELD FERMION REP ======= \n");

    // Setup fields on GPU
    int return_val = 0;
    suNf_field *f = alloc_suNf_field(&glattice);
    random_suNf_field_cpu(f);

    // Evaluate sqnorm in the beginning
    double sqnorm_start = sqnorm_suNf_field_cpu(f);
    lprintf("SANITY CHECK", 0, "[In field CPU copy norm unequal zero: %0.2e]\n", sqnorm_start);
    copy_to_gpu_suNf_field(f);

    // Execute communications
    start_sendrecv_suNf_field(f);
    complete_sendrecv_suNf_field(f);

    // Evaluate sqnorm after comms
    copy_from_gpu_suNf_field(f);
    double sqnorm_end = sqnorm_suNf_field_cpu(f);
    lprintf("SANITY CHECK", 0, "[Out field CPU copy norm unequal zero: %0.2e]\n", sqnorm_end);

    return_val += check_finiteness(sqnorm_start);
    return_val += check_finiteness(sqnorm_end);
    return_val += check_diff_norm_zero(sqnorm_start - sqnorm_end);

    free_suNf_field(f);

    return return_val;
}

int test_comms_suNf_field_flt() {
    lprintf("INFO", 0, " ======= TEST GFIELD FERMION REP SINGLE PRECISION ======= \n");

    // Setup fields on GPU
    int return_val = 0;
    suNf_field_flt *f = alloc_suNf_field_flt(&glattice);
    random_suNf_field_flt_cpu(f);

    // Evaluate sqnorm in the beginning
    double sqnorm_start = sqnorm_suNf_field_flt_cpu(f);
    lprintf("SANITY CHECK", 0, "[In field CPU copy norm unequal zero: %0.2e]\n", sqnorm_start);
    copy_to_gpu_suNf_field_flt(f);

    // Execute communications
    start_sendrecv_suNf_field_flt(f);
    complete_sendrecv_suNf_field_flt(f);

    // Evaluate sqnorm after comms
    copy_from_gpu_suNf_field_flt(f);
    double sqnorm_end = sqnorm_suNf_field_flt_cpu(f);
    lprintf("SANITY CHECK", 0, "[Out field CPU copy norm unequal zero: %0.2e]\n", sqnorm_end);

    return_val += check_finiteness(sqnorm_start);
    return_val += check_finiteness(sqnorm_end);
    return_val += check_diff_norm_zero(sqnorm_start - sqnorm_end);

    free_suNf_field_flt(f);

    return return_val;
}

int test_comms_suNg_scalar_field() {
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
    start_sendrecv_suNg_scalar_field(f);
    complete_sendrecv_suNg_scalar_field(f);

    // Evaluate sqnorm after comms
    copy_from_gpu_suNg_scalar_field(f);
    double sqnorm_end = sqnorm_suNg_scalar_field_cpu(f);
    lprintf("SANITY CHECK", 0, "[Out field CPU copy norm unequal zero: %0.2e]\n", sqnorm_end);

    return_val += check_finiteness(sqnorm_start);
    return_val += check_finiteness(sqnorm_end);
    return_val += check_diff_norm_zero(sqnorm_start - sqnorm_end);

    free_suNg_scalar_field(f);

    return return_val;
}

int test_comms_suNg_av_field() {
    lprintf("INFO", 0, " ======= TEST ALGEBRA VECTOR FIELD ======= \n");

    // Setup fields on GPU
    int return_val = 0;
    suNg_av_field *f = alloc_suNg_av_field(&glattice);
    random_suNg_av_field_cpu(f);

    // Evaluate sqnorm in the beginning
    double sqnorm_start = sqnorm_suNg_av_field_cpu(f);
    lprintf("SANITY CHECK", 0, "[In field CPU copy norm unequal zero: %0.2e]\n", sqnorm_start);
    copy_to_gpu_suNg_av_field(f);

    // Execute communications
    start_sendrecv_suNg_av_field(f);
    complete_sendrecv_suNg_av_field(f);

    // Evaluate sqnorm after comms
    copy_from_gpu_suNg_av_field(f);
    double sqnorm_end = sqnorm_suNg_av_field_cpu(f);
    lprintf("SANITY CHECK", 0, "[Out field CPU copy norm unequal zero: %0.2e]\n", sqnorm_end);

    return_val += check_finiteness(sqnorm_start);
    return_val += check_finiteness(sqnorm_end);
    return_val += check_diff_norm_zero(sqnorm_start - sqnorm_end);

    free_suNg_av_field(f);

    return return_val;
}

int test_comms_gtransf() {
    lprintf("INFO", 0, " ======= TEST GAUGE TRANSFORMATION ======= \n");

    // Setup fields on GPU
    int return_val = 0;
    gtransf *f = alloc_gtransf(&glattice);
    random_gtransf_cpu(f);

    // Evaluate sqnorm in the beginning
    double sqnorm_start = sqnorm_gtransf_cpu(f);
    lprintf("SANITY CHECK", 0, "[In field CPU copy norm unequal zero: %0.2e]\n", sqnorm_start);
    copy_to_gpu_gtransf(f);

    // Execute communications
    start_sendrecv_gtransf(f);
    complete_sendrecv_gtransf(f);

    // Evaluate sqnorm after comms
    copy_from_gpu_gtransf(f);
    double sqnorm_end = sqnorm_gtransf_cpu(f);
    lprintf("SANITY CHECK", 0, "[Out field CPU copy norm unequal zero: %0.2e]\n", sqnorm_end);

    return_val += check_finiteness(sqnorm_start);
    return_val += check_finiteness(sqnorm_end);
    return_val += check_diff_norm_zero(sqnorm_start - sqnorm_end);

    free_gtransf(f);

    return return_val;
}

int test_comms_ldl_field() {
    lprintf("INFO", 0, " ======= TEST CLOVER LDL ======= \n");

    // Setup fields on GPU
    int return_val = 0;
    ldl_field *f = alloc_ldl_field(&glattice);
    random_ldl_field_cpu(f);

    // Evaluate sqnorm in the beginning
    double sqnorm_start = sqnorm_ldl_field_cpu(f);
    lprintf("SANITY CHECK", 0, "[In field CPU copy norm unequal zero: %0.2e]\n", sqnorm_start);
    copy_to_gpu_ldl_field(f);

    // Execute communications
    start_sendrecv_ldl_field(f);
    complete_sendrecv_ldl_field(f);

    // Evaluate sqnorm after comms
    copy_from_gpu_ldl_field(f);
    double sqnorm_end = sqnorm_ldl_field_cpu(f);
    lprintf("SANITY CHECK", 0, "[Out field CPU copy norm unequal zero: %0.2e]\n", sqnorm_end);

    return_val += check_finiteness(sqnorm_start);
    return_val += check_finiteness(sqnorm_end);
    return_val += check_diff_norm_zero(sqnorm_start - sqnorm_end);

    free_ldl_field(f);

    return return_val;
}

int test_comms_clover_term() {
    lprintf("INFO", 0, " ======= TEST CLOVER TERM ======= \n");

    // Setup fields on GPU
    int return_val = 0;
    clover_term *f = alloc_clover_term(&glattice);
    random_clover_term_cpu(f);

    // Evaluate sqnorm in the beginning
    double sqnorm_start = sqnorm_clover_term_cpu(f);
    lprintf("SANITY CHECK", 0, "[In field CPU copy norm unequal zero: %0.2e]\n", sqnorm_start);
    copy_to_gpu_clover_term(f);

    // Execute communications
    start_sendrecv_clover_term(f);
    complete_sendrecv_clover_term(f);

    // Evaluate sqnorm after comms
    copy_from_gpu_clover_term(f);
    double sqnorm_end = sqnorm_clover_term_cpu(f);
    lprintf("SANITY CHECK", 0, "[Out field CPU copy norm unequal zero: %0.2e]\n", sqnorm_end);

    return_val += check_finiteness(sqnorm_start);
    return_val += check_finiteness(sqnorm_end);
    return_val += check_diff_norm_zero(sqnorm_start - sqnorm_end);

    free_clover_term(f);

    return return_val;
}

int test_comms_clover_force() {
    lprintf("INFO", 0, " ======= TEST CLOVER FORCE ======= \n");

    // Setup fields on GPU
    int return_val = 0;
    clover_force *f = alloc_clover_force(&glattice);
    random_clover_force_cpu(f);

    // Evaluate sqnorm in the beginning
    double sqnorm_start = sqnorm_clover_force_cpu(f);
    lprintf("SANITY CHECK", 0, "[In field CPU copy norm unequal zero: %0.2e]\n", sqnorm_start);
    copy_to_gpu_clover_force(f);

    // Execute communications
    start_sendrecv_clover_force(f);
    complete_sendrecv_clover_force(f);

    // Evaluate sqnorm after comms
    copy_from_gpu_clover_force(f);
    double sqnorm_end = sqnorm_clover_force_cpu(f);
    lprintf("SANITY CHECK", 0, "[Out field CPU copy norm unequal zero: %0.2e]\n", sqnorm_end);

    return_val += check_finiteness(sqnorm_start);
    return_val += check_finiteness(sqnorm_end);
    return_val += check_diff_norm_zero(sqnorm_start - sqnorm_end);

    free_clover_force(f);

    return return_val;
}

int test_comms_staple_field() {
    lprintf("INFO", 0, " ======= TEST STAPLE FIELD ======= \n");

    // Setup fields on GPU
    int return_val = 0;
    staple_field *f = alloc_staple_field(&glattice);
    random_staple_field_cpu(f);

    // Evaluate sqnorm in the beginning
    double sqnorm_start = sqnorm_staple_field_cpu(f);
    lprintf("SANITY CHECK", 0, "[In field CPU copy norm unequal zero: %0.2e]\n", sqnorm_start);
    copy_to_gpu_staple_field(f);

    // Execute communications
    start_sendrecv_staple_field(f);
    complete_sendrecv_staple_field(f);

    // Evaluate sqnorm after comms
    copy_from_gpu_staple_field(f);
    double sqnorm_end = sqnorm_staple_field_cpu(f);
    lprintf("SANITY CHECK", 0, "[Out field CPU copy norm unequal zero: %0.2e]\n", sqnorm_end);

    return_val += check_finiteness(sqnorm_start);
    return_val += check_finiteness(sqnorm_end);
    return_val += check_diff_norm_zero(sqnorm_start - sqnorm_end);

    free_staple_field(f);

    return return_val;
}
