/*******************************************************************************
*
* Check that after allocating a field, we can write to and read from it on 
* the CPU. (Test should work both WITH_GPU and !WITH_GPU, but this is
* important to check, whether the CPU part of the WITH_GPU code works)
*
*******************************************************************************/

#include "libhr.h"

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

int main(int argc, char *argv[]) {
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
    return_val += test_avfield_allocation();
    return_val += test_clover_term_allocation();
    return_val += test_clover_force_allocation();
    return_val += test_sfield_allocation();

    // Single precision test block
    return_val += test_gfield_flt_allocation();
    return_val += test_gfield_f_flt_allocation();
    return_val += test_spinor_field_flt_allocation();

    // Finalize and return
    finalize_process();
    return return_val;
}

int test_gfield_allocation() {
    lprintf("INFO", 0, " ======= TEST GAUGE FIELD ======= \n");
    int return_val = 0;
    suNg_field *f = alloc_gfield(&glattice);
    random_gfield_cpu(f);
    double sqnorm = sqnorm_gfield_cpu(f);
    return_val += check_finiteness(sqnorm);
    free_gfield(f);
    return return_val;
}

int test_gfield_f_allocation() {
    lprintf("INFO", 0, " ======= TEST REPRESENTED GAUGE FIELD ======= \n");
    int return_val = 0;
    suNf_field *f = alloc_gfield_f(&glattice);
    random_gfield_f_cpu(f);
    double sqnorm = sqnorm_gfield_f_cpu(f);
    return_val += check_finiteness(sqnorm);
    free_gfield_f(f);
    return return_val;
}

int test_gfield_flt_allocation() {
    lprintf("INFO", 0, " ======= TEST GAUGE FIELD SINGLE PRECISION ======= \n");
    int return_val = 0;
    suNg_field_flt *f = alloc_gfield_flt(&glattice);
    random_gfield_flt_cpu(f);
    float sqnorm = sqnorm_gfield_flt_cpu(f);
    return_val += check_finiteness(sqnorm);
    free_gfield_flt(f);
    return return_val;
}

int test_gfield_f_flt_allocation() {
    lprintf("INFO", 0, " ======= TEST REPRESENTED GAUGE FIELD SINGLE PRECISION ======= \n");
    int return_val = 0;
    suNf_field_flt *f = alloc_gfield_f_flt(&glattice);
    random_gfield_f_flt_cpu(f);
    float sqnorm = sqnorm_gfield_f_flt_cpu(f);
    return_val += check_finiteness(sqnorm);
    free_gfield_f_flt(f);
    return return_val;
}

int test_avfield_allocation() {
    lprintf("INFO", 0, " ======= TEST SU(N) ALGEBRA VECTOR FIELD ======= \n");
    int return_val = 0;
    suNg_av_field *f = alloc_avfield(&glattice);
    random_avfield_cpu(f);
    double sqnorm = sqnorm_avfield_cpu(f);
    return_val += check_finiteness(sqnorm);
    free_avfield(f);
    return return_val;
}

int test_gtransf_allocation() {
    lprintf("INFO", 0, " ======= GAUGE TRANSFORMATION ======= \n");
    int return_val = 0;
    suNg_field *f = alloc_gtransf(&glattice);
    random_gtransf_cpu(f);
    double sqnorm = sqnorm_gtransf_cpu(f);
    return_val += check_finiteness(sqnorm);
    free_gtransf(f);
    return return_val;
}

int test_clover_term_allocation() {
    lprintf("INFO", 0, " ======= CLOVER TERM ======= \n");
    int return_val = 0;
    suNfc_field *f = alloc_clover_term(&glattice);
    random_clover_term_cpu(f);
    double sqnorm = sqnorm_clover_term_cpu(f);
    return_val += check_finiteness(sqnorm);
    free_clover_term(f);
    return return_val;
}

int test_clover_force_allocation() {
    lprintf("INFO", 0, " ======= CLOVER FORCE ======= \n");
    int return_val = 0;
    suNf_field *f = alloc_clover_force(&glattice);
    random_clover_force_cpu(f);
    double sqnorm = sqnorm_clover_force_cpu(f);
    return_val += check_finiteness(sqnorm);
    free_clover_force(f);
    return return_val;
}

int test_spinor_field_allocation() {
    lprintf("INFO", 0, " ======= TEST SPINOR FIELD ======= \n");
    int return_val = 0;
    spinor_field *f = alloc_spinor_field_f(1, &glattice);
    gaussian_spinor_field(f);
    double sqnorm = spinor_field_sqnorm_f_cpu(f);
    return_val += check_finiteness(sqnorm);
    free_spinor_field_f(f);
    return return_val;
}

int test_spinor_field_flt_allocation() {
    lprintf("INFO", 0, " ======= TEST SPINOR FIELD SINGLE PRECISION ======= \n");
    int return_val = 0;
    spinor_field_flt *f = alloc_spinor_field_f_flt(1, &glattice);
    gaussian_spinor_field_flt(f);
    double sqnorm = spinor_field_sqnorm_f_flt_cpu(f);
    return_val += check_finiteness(sqnorm);
    free_spinor_field_f_flt(f);
    return return_val;
}

int test_sfield_allocation() {
    lprintf("INFO", 0, " ======= TEST SFIELD ======= \n");
    int return_val = 0;
    scalar_field *f = alloc_sfield(1, &glattice);
    random_sfield_cpu(f);
    double sqnorm = sqnorm_sfield_cpu(f);
    return_val += check_finiteness(sqnorm);
    free_sfield(f);
    return return_val;
}
