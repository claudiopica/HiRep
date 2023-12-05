/*******************************************************************************
*
* NOCOMPILE= !WITH_GPU
*
* Check that after allocating a field, we can write to and read from it.
*
*******************************************************************************/

#include "libhr.h"

// Double precision
int test_suNg_field_allocation();
int test_suNf_field_allocation();
int test_spinor_field_allocation(geometry_descriptor *);
int test_suNg_av_field_allocation();
int test_clover_term_allocation();
int test_clover_force_allocation();
int test_scalar_field_allocation(geometry_descriptor *);
int test_suNg_scalar_field_allocation();
int test_gtransf_allocation();
//int test_clover_ldl_allocation();
int test_staple_field_allocation();

// Single precision
int test_suNg_field_flt_allocation();
int test_suNf_field_flt_allocation();
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
    return_val += test_suNg_field_allocation();
    return_val += test_suNf_field_allocation();
    return_val += test_spinor_field_allocation(&glattice);
    return_val += test_suNg_av_field_allocation();
    return_val += test_clover_term_allocation();
    return_val += test_clover_force_allocation();
    return_val += test_scalar_field_allocation(&glattice);
    return_val += test_suNg_scalar_field_allocation();
    return_val += test_gtransf_allocation();
    //return_val += test_clover_ldl_allocation();
    return_val += test_staple_field_allocation();

    // Single precision test block
    return_val += test_suNg_field_flt_allocation();
    return_val += test_suNf_field_flt_allocation();
    return_val += test_spinor_field_flt_allocation(&glattice);

    // Spinor-like with even parity
    lprintf("INFO", 0, "\n\n Even parity\n\n");
    return_val += test_spinor_field_allocation(&glat_even);
    return_val += test_spinor_field_flt_allocation(&glat_even);
    return_val += test_scalar_field_allocation(&glat_even);

    // Spinor-like with odd parity
    lprintf("INFO", 0, "\n\n Odd parity\n\n");
    return_val += test_spinor_field_allocation(&glat_odd);
    return_val += test_spinor_field_flt_allocation(&glat_odd);
    return_val += test_scalar_field_allocation(&glat_odd);

    // Finalize and return
    finalize_process();
    return return_val;
}

int test_suNg_field_allocation() {
    lprintf("INFO", 0, " ======= TEST GAUGE FIELD ======= \n");
    int return_val = 0;
    suNg_field *f = alloc_suNg_field(&glattice);
    random_suNg_field_cpu(f);
    copy_to_gpu_suNg_field(f);
    zero_suNg_field_cpu(f);
    copy_from_gpu_suNg_field(f);
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
    copy_to_gpu_suNf_field(f);
    zero_suNf_field_cpu(f);
    copy_from_gpu_suNf_field(f);
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
    copy_to_gpu_suNg_field_flt(f);
    zero_suNg_field_flt_cpu(f);
    copy_from_gpu_suNg_field_flt(f);
    double sqnorm = sqnorm_suNg_field_flt_cpu(f);
    return_val += check_finiteness(sqnorm);
    free_suNg_field_flt(f);
    return return_val;
}

int test_suNf_field_flt_allocation() {
    lprintf("INFO", 0, " ======= TEST REPRESENTED GAUGE FIELD SINGLE PRECISION ======= \n");
    int return_val = 0;
    suNf_field_flt *f = alloc_suNf_field_flt(&glattice);
    random_suNf_field_flt_cpu(f);
    copy_to_gpu_suNf_field_flt(f);
    zero_suNf_field_flt_cpu(f);
    copy_from_gpu_suNf_field_flt(f);
    double sqnorm = sqnorm_suNf_field_flt_cpu(f);
    return_val += check_finiteness(sqnorm);
    free_suNf_field_flt(f);
    return return_val;
}

int test_suNg_av_field_allocation() {
    lprintf("INFO", 0, " ======= TEST SU(N) ALGEBRA VECTOR FIELD ======= \n");
    int return_val = 0;
    suNg_av_field *f = alloc_suNg_av_field(&glattice);
    random_suNg_av_field_cpu(f);
    copy_to_gpu_suNg_av_field(f);
    zero_suNg_av_field_cpu(f);
    copy_from_gpu_suNg_av_field(f);
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
    clover_term *f = alloc_clover_term(&glattice);
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
    clover_force *f = alloc_clover_force(&glattice);
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
    spinor_field *f = alloc_spinor_field(1, gd);
    gaussian_spinor_field(f);
    copy_to_gpu_spinor_field(f);
    zero_spinor_field_cpu(f);
    copy_from_gpu_spinor_field(f);
    double sqnorm = sqnorm_spinor_field_cpu(f);
    return_val += check_finiteness(sqnorm);
    free_spinor_field(f);
    return return_val;
}

int test_spinor_field_flt_allocation(geometry_descriptor *gd) {
    lprintf("INFO", 0, " ======= TEST SPINOR FIELD SINGLE PRECISION ======= \n");
    int return_val = 0;
    spinor_field_flt *f = alloc_spinor_field_flt(1, gd);
    gaussian_spinor_field_flt(f);
    copy_to_gpu_spinor_field_flt(f);
    zero_spinor_field_flt_cpu(f);
    copy_from_gpu_spinor_field_flt(f);
    double sqnorm = sqnorm_spinor_field_flt_cpu(f);
    return_val += check_finiteness(sqnorm);
    free_spinor_field_flt(f);
    return return_val;
}

int test_scalar_field_allocation(geometry_descriptor *gd) {
    lprintf("INFO", 0, " ======= TEST SCALAR FIELD ======= \n");
    int return_val = 0;
    scalar_field *f = alloc_scalar_field(1, gd);
    random_scalar_field_cpu(f);
    copy_to_gpu_scalar_field(f);
    zero_scalar_field_cpu(f);
    copy_from_gpu_scalar_field(f);
    double sqnorm = sqnorm_scalar_field_cpu(f);
    return_val += check_finiteness(sqnorm);
    free_scalar_field(f);
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

#if 0
int test_clover_ldl_allocation() {
    lprintf("INFO", 0, " ======= CLOVER LDL ======= \n");
    int return_val = 0;
    ldl_field *f = alloc_clover_ldl(&glattice);
    random_clover_ldl_cpu(f);
    copy_to_gpu_clover_ldl(f);
    zero_clover_ldl_cpu(f);
    copy_from_gpu_clover_ldl(f);
    double sqnorm = sqnorm_ldl_field_cpu(f);
    return_val += check_finiteness(sqnorm);
    free_ldl_field(f);
    return return_val;
}
#endif

int test_staple_field_allocation() {
    lprintf("INFO", 0, " ======= STAPLE FIELD ======= \n");
    int return_val = 0;
    staple_field *f = alloc_staple_field(&glattice);
    random_staple_field_cpu(f);
    copy_to_gpu_staple_field(f);
    zero_staple_field_cpu(f);
    copy_from_gpu_staple_field(f);
    double sqnorm = sqnorm_staple_field_cpu(f);
    return_val += check_finiteness(sqnorm);
    free_staple_field(f);
    return return_val;
}
