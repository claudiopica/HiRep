/*******************************************************************************
*
* Check that after allocating a field, we can write to and read from it on 
* the CPU. (Test should work both WITH_GPU and !WITH_GPU, but this is
* important to check, whether the CPU part of the WITH_GPU code works)
*
*******************************************************************************/

#include "libhr.h"

// Double precision
int test_suNg_field_allocation();
int test_suNf_field_allocation();
int test_suNg_scalar_field_allocation();
int test_suNg_av_field_allocation();
int test_gtransf_allocation();
int test_clover_term_allocation();
int test_clover_force_allocation();
int test_spinor_field_allocation();
int test_scalar_field_allocation();
int test_ldl_field_allocation();
int test_staple_field_allocation();

// Single precision
int test_suNg_field_flt_allocation();
int test_suNf_field_flt_allocation();
int test_spinor_field_flt_allocation();

int main(int argc, char *argv[]) {
    // Init
    int return_val = 0;

    // Setup process and communication
    logger_map("DEBUG", "debug");
    setup_process(&argc, &argv);
    test_setup();

    // Double precision test block
    return_val += test_suNg_field_allocation();
    return_val += test_suNf_field_allocation();
    return_val += test_spinor_field_allocation();
    return_val += test_suNg_av_field_allocation();
    return_val += test_clover_term_allocation();
    return_val += test_clover_force_allocation();
    return_val += test_scalar_field_allocation();
    return_val += test_suNg_scalar_field_allocation();
    return_val += test_gtransf_allocation();
    return_val += test_ldl_field_allocation();
    return_val += test_staple_field_allocation();

    // Single precision test block
    return_val += test_suNg_field_flt_allocation();
    return_val += test_suNf_field_flt_allocation();
    return_val += test_spinor_field_flt_allocation();

    // Finalize and return
    finalize_process();
    return return_val;
}

int test_suNg_field_allocation() {
    lprintf("INFO", 0, " ======= TEST GAUGE FIELD ======= \n");
    int return_val = 0;
    suNg_field *f = alloc_suNg_field(&glattice);
    random_suNg_field_cpu(f);
    double sqnorm = sqnorm_suNg_field_cpu(f);
    return_val += check_finiteness(sqnorm);
    free_suNg_field(f);
    return return_val;
}

int test_suNf_field_allocation() {
    lprintf("INFO", 0, " ======= TEST REPRESENTED GAUGE FIELD ======= \n");
    int return_val = 0;
    suNf_field *f = alloc_suNf_field(&glattice);
    random_suNf_field_cpu(f);
    double sqnorm = sqnorm_suNf_field_cpu(f);
    return_val += check_finiteness(sqnorm);
    free_suNf_field(f);
    return return_val;
}

int test_suNg_field_flt_allocation() {
    lprintf("INFO", 0, " ======= TEST GAUGE FIELD SINGLE PRECISION ======= \n");
    int return_val = 0;
    suNg_field_flt *f = alloc_suNg_field_flt(&glattice);
    random_suNg_field_flt_cpu(f);
    float sqnorm = sqnorm_suNg_field_flt_cpu(f);
    return_val += check_finiteness(sqnorm);
    free_suNg_field_flt(f);
    return return_val;
}

int test_suNf_field_flt_allocation() {
    lprintf("INFO", 0, " ======= TEST REPRESENTED GAUGE FIELD SINGLE PRECISION ======= \n");
    int return_val = 0;
    suNf_field_flt *f = alloc_suNf_field_flt(&glattice);
    random_suNf_field_flt_cpu(f);
    float sqnorm = sqnorm_suNf_field_flt_cpu(f);
    return_val += check_finiteness(sqnorm);
    free_suNf_field_flt(f);
    return return_val;
}

int test_suNg_av_field_allocation() {
    lprintf("INFO", 0, " ======= TEST SU(N) ALGEBRA VECTOR FIELD ======= \n");
    int return_val = 0;
    suNg_av_field *f = alloc_suNg_av_field(&glattice);
    random_suNg_av_field_cpu(f);
    double sqnorm = sqnorm_suNg_av_field_cpu(f);
    return_val += check_finiteness(sqnorm);
    free_suNg_av_field(f);
    return return_val;
}

int test_gtransf_allocation() {
    lprintf("INFO", 0, " ======= GAUGE TRANSFORMATION ======= \n");
    int return_val = 0;
    gtransf *f = alloc_gtransf(&glattice);
    random_gtransf_cpu(f);
    double sqnorm = sqnorm_gtransf_cpu(f);
    return_val += check_finiteness(sqnorm);
    free_gtransf(f);
    return return_val;
}

int test_clover_term_allocation() {
    lprintf("INFO", 0, " ======= CLOVER TERM ======= \n");
    int return_val = 0;
    clover_term *f = alloc_clover_term(&glattice);
    random_clover_term_cpu(f);
    double sqnorm = sqnorm_clover_term_cpu(f);
    return_val += check_finiteness(sqnorm);
    free_clover_term(f);
    return return_val;
}

int test_clover_force_allocation() {
    lprintf("INFO", 0, " ======= CLOVER FORCE ======= \n");
    int return_val = 0;
    clover_force *f = alloc_clover_force(&glattice);
    random_clover_force_cpu(f);
    double sqnorm = sqnorm_clover_force_cpu(f);
    return_val += check_finiteness(sqnorm);
    free_clover_force(f);
    return return_val;
}

int test_spinor_field_allocation() {
    lprintf("INFO", 0, " ======= TEST SPINOR FIELD ======= \n");
    int return_val = 0;
    spinor_field *f = alloc_spinor_field(1, &glattice);
    gaussian_spinor_field(f);
    double sqnorm = spinor_field_sqnorm_f_cpu(f);
    return_val += check_finiteness(sqnorm);
    free_spinor_field(f);
    return return_val;
}

int test_spinor_field_flt_allocation() {
    lprintf("INFO", 0, " ======= TEST SPINOR FIELD SINGLE PRECISION ======= \n");
    int return_val = 0;
    spinor_field_flt *f = alloc_spinor_field_flt(1, &glattice);
    gaussian_spinor_field_flt(f);
    double sqnorm = spinor_field_sqnorm_f_flt_cpu(f);
    return_val += check_finiteness(sqnorm);
    free_spinor_field_flt(f);
    return return_val;
}

int test_scalar_field_allocation() {
    lprintf("INFO", 0, " ======= TEST SFIELD ======= \n");
    int return_val = 0;
    scalar_field *f = alloc_scalar_field(1, &glattice);
    random_scalar_field_cpu(f);
    double sqnorm = sqnorm_scalar_field_cpu(f);
    return_val += check_finiteness(sqnorm);
    free_scalar_field(f);
    return return_val;
}

int test_suNg_scalar_field_allocation() {
    lprintf("INFO", 0, " ======= TEST SU(NG) SCALAR FIELD ======= \n");
    int return_val = 0;
    suNg_scalar_field *f = alloc_suNg_scalar_field(&glattice);
    random_suNg_scalar_field_cpu(f);
    double sqnorm = sqnorm_suNg_scalar_field_cpu(f);
    return_val += check_finiteness(sqnorm);
    free_suNg_scalar_field(f);
    return return_val;
}

int test_ldl_field_allocation() {
    lprintf("INFO", 0, " ======= TEST CLOVER LDL ======= \n");
    int return_val = 0;
    ldl_field *f = alloc_ldl_field(&glattice);
    random_ldl_field_cpu(f);
    double sqnorm = sqnorm_ldl_field_cpu(f);
    return_val += check_finiteness(sqnorm);
    free_ldl_field(f);
    return return_val;
}

int test_staple_field_allocation() {
    lprintf("INFO", 0, " ======= TEST STAPLE FIELD ======= \n");
    int return_val = 0;
    staple_field *f = alloc_staple_field(&glattice);
    random_staple_field_cpu(f);
    double sqnorm = sqnorm_staple_field_cpu(f);
    return_val += check_finiteness(sqnorm);
    free_staple_field(f);
    return return_val;
}
