/*******************************************************************************
*
* NOCOMPILE= !WITH_GPU
*
* Check that after allocating a field, we can write to and read from it.
*
*******************************************************************************/

#include "libhr.h"

// Double precision
int test_gfield_allocation();
int test_gfield_f_allocation();
int test_spinor_field_allocation(geometry_descriptor *);
int test_avfield_allocation();
int test_clover_term_allocation();
int test_clover_force_allocation();
int test_sfield_allocation(geometry_descriptor *);
int test_suNg_scalar_field_allocation();
int test_gtransf_allocation();
int test_clover_ldl_allocation();
int test_staple_field_allocation();

// Single precision
int test_gfield_flt_allocation();
int test_gfield_f_flt_allocation();
int test_spinor_field_flt_allocation(geometry_descriptor *);

int main(int argc, char *argv[]) {
    // Init
    int return_val = 0;

    // Setup process and communication
    logger_map("DEBUG", "debug");
    setup_process(&argc, &argv);
    test_setup();

    // Double precision test block
    lprintf("INFO", 0, "\n\nFull lattice tests\n\n");
    return_val += test_gfield_allocation();
    return_val += test_gfield_f_allocation();
    return_val += test_spinor_field_allocation(&glattice);
    return_val += test_avfield_allocation();
    return_val += test_clover_term_allocation();
    return_val += test_clover_force_allocation();
    return_val += test_sfield_allocation(&glattice);
    return_val += test_suNg_scalar_field_allocation();
    return_val += test_gtransf_allocation();
    return_val += test_clover_ldl_allocation();
    return_val += test_staple_field_allocation();

    // Single precision test block
    return_val += test_gfield_flt_allocation();
    return_val += test_gfield_f_flt_allocation();
    return_val += test_spinor_field_flt_allocation(&glattice);

    // Spinor-like with even parity
    lprintf("INFO", 0, "\n\n Even parity\n\n");
    return_val += test_spinor_field_allocation(&glat_even);
    return_val += test_spinor_field_flt_allocation(&glat_even);
    return_val += test_sfield_allocation(&glat_even);

    // Spinor-like with odd parity
    lprintf("INFO", 0, "\n\n Odd parity\n\n");
    return_val += test_spinor_field_allocation(&glat_odd);
    return_val += test_spinor_field_flt_allocation(&glat_odd);
    return_val += test_sfield_allocation(&glat_odd);

    // Finalize and return
    finalize_process();
    return return_val;
}

int test_gfield_allocation() {
    lprintf("INFO", 0, " ======= TEST GAUGE FIELD ======= \n");
    int return_val = 0;
    suNg_field *f = alloc_gfield(&glattice);
    random_gfield_cpu(f);
    copy_to_gpu_gfield(f);
    zero_gfield_cpu(f);
    copy_from_gpu_gfield(f);
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
    copy_to_gpu_gfield_f(f);
    zero_gfield_f_cpu(f);
    copy_from_gpu_gfield_f(f);
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
    copy_to_gpu_gfield_flt(f);
    zero_gfield_flt_cpu(f);
    copy_from_gpu_gfield_flt(f);
    double sqnorm = sqnorm_gfield_flt_cpu(f);
    return_val += check_finiteness(sqnorm);
    free_gfield_flt(f);
    return return_val;
}

int test_gfield_f_flt_allocation() {
    lprintf("INFO", 0, " ======= TEST REPRESENTED GAUGE FIELD SINGLE PRECISION ======= \n");
    int return_val = 0;
    suNf_field_flt *f = alloc_gfield_f_flt(&glattice);
    random_gfield_f_flt_cpu(f);
    copy_to_gpu_gfield_f_flt(f);
    zero_gfield_f_flt_cpu(f);
    copy_from_gpu_gfield_f_flt(f);
    double sqnorm = sqnorm_gfield_f_flt_cpu(f);
    return_val += check_finiteness(sqnorm);
    free_gfield_f_flt(f);
    return return_val;
}

int test_avfield_allocation() {
    lprintf("INFO", 0, " ======= TEST SU(N) ALGEBRA VECTOR FIELD ======= \n");
    int return_val = 0;
    suNg_av_field *f = alloc_avfield(&glattice);
    random_avfield_cpu(f);
    copy_to_gpu_avfield(f);
    zero_avfield_cpu(f);
    copy_from_gpu_avfield(f);
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
    copy_to_gpu_gtransf(f);
    zero_gtransf_cpu(f);
    copy_from_gpu_gtransf(f);
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
    copy_to_gpu_clover_term(f);
    zero_clover_term_cpu(f);
    copy_from_gpu_clover_term(f);
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
    copy_to_gpu_clover_force(f);
    zero_clover_force_cpu(f);
    copy_from_gpu_clover_force(f);
    double sqnorm = sqnorm_clover_force_cpu(f);
    return_val += check_finiteness(sqnorm);
    free_clover_force(f);
    return return_val;
}

int test_spinor_field_allocation(geometry_descriptor *gd) {
    lprintf("INFO", 0, " ======= TEST SPINOR FIELD ======= \n");
    int return_val = 0;
    spinor_field *f = alloc_spinor_field_f(1, gd);
    gaussian_spinor_field(f);
    copy_to_gpu_spinor_field_f(f);
    spinor_field_zero_f_cpu(f);
    copy_from_gpu_spinor_field_f(f);
    double sqnorm = spinor_field_sqnorm_f_cpu(f);
    return_val += check_finiteness(sqnorm);
    free_spinor_field_f(f);
    return return_val;
}

int test_spinor_field_flt_allocation(geometry_descriptor *gd) {
    lprintf("INFO", 0, " ======= TEST SPINOR FIELD SINGLE PRECISION ======= \n");
    int return_val = 0;
    spinor_field_flt *f = alloc_spinor_field_f_flt(1, gd);
    gaussian_spinor_field_flt(f);
    copy_to_gpu_spinor_field_f_flt(f);
    spinor_field_zero_f_flt_cpu(f);
    copy_from_gpu_spinor_field_f_flt(f);
    double sqnorm = spinor_field_sqnorm_f_flt_cpu(f);
    return_val += check_finiteness(sqnorm);
    free_spinor_field_f_flt(f);
    return return_val;
}

int test_sfield_allocation(geometry_descriptor *gd) {
    lprintf("INFO", 0, " ======= TEST SCALAR FIELD ======= \n");
    int return_val = 0;
    scalar_field *f = alloc_sfield(1, gd);
    random_sfield_cpu(f);
    copy_to_gpu_sfield(f);
    zero_sfield_cpu(f);
    copy_from_gpu_sfield(f);
    double sqnorm = sqnorm_sfield_cpu(f);
    return_val += check_finiteness(sqnorm);
    free_sfield(f);
    return return_val;
}

int test_suNg_scalar_field_allocation() {
    lprintf("INFO", 0, " ======= SU(NG) SCALAR FIELD ======= \n");
    int return_val = 0;
    suNg_scalar_field *f = alloc_suNg_scalar_field(&glattice);
    random_suNg_scalar_field_cpu(f);
    copy_to_gpu_suNg_scalar_field(f);
    zero_suNg_scalar_field_cpu(f);
    copy_from_gpu_suNg_scalar_field(f);
    double sqnorm = sqnorm_suNg_scalar_field_cpu(f);
    return_val += check_finiteness(sqnorm);
    free_suNg_scalar_field(f);
    return return_val;
}

int test_clover_ldl_allocation() {
    lprintf("INFO", 0, " ======= CLOVER LDL ======= \n");
    int return_val = 0;
    ldl_field *f = alloc_clover_ldl(&glattice);
    random_clover_ldl_cpu(f);
    copy_to_gpu_clover_ldl(f);
    zero_clover_ldl_cpu(f);
    copy_from_gpu_clover_ldl(f);
    double sqnorm = sqnorm_clover_ldl_cpu(f);
    return_val += check_finiteness(sqnorm);
    free_clover_ldl(f);
    return return_val;
}

int test_staple_field_allocation() {
    lprintf("INFO", 0, " ======= STAPLE FIELD ======= \n");
    int return_val = 0;
    suNg_field *f = alloc_staple_field(&glattice);
    random_staple_field_cpu(f);
    copy_to_gpu_staple_field(f);
    zero_staple_field_cpu(f);
    copy_from_gpu_staple_field(f);
    double sqnorm = sqnorm_staple_field_cpu(f);
    return_val += check_finiteness(sqnorm);
    free_staple_field(f);
    return return_val;
}
