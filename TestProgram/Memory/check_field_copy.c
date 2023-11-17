/*******************************************************************************
*
* NOCOMPILE= !WITH_GPU
*
* Check that copy back and forth works for all fields.
*
*******************************************************************************/

#include "libhr.h"

// Double precision
int test_bijectivity_suNg_field();
int test_bijectivity_suNf_field();
int test_bijectivity_suNg_scalar_field();
int test_bijectivity_suNg_av_field();
int test_bijectivity_gtransf();
int test_bijectivity_ldl_field();
int test_bijectivity_clover_term();
int test_bijectivity_clover_force();
int test_bijectivity_spinor_field(geometry_descriptor *);
int test_bijectivity_scalar_field(geometry_descriptor *);
int test_bijectivity_staple_field();

// Single precision
int test_bijectivity_suNg_field_flt();
int test_bijectivity_suNf_field_flt();
int test_bijectivity_spinor_field_flt(geometry_descriptor *);

int main(int argc, char *argv[]) {
    // Init
    int return_val = 0;

    // Setup process and communication
    logger_map("DEBUG", "debug");
    setup_process(&argc, &argv);
    test_setup();

    // Test block

    lprintf("INFO", 0, "\n\n Testing Full Lattice \n\n");
    /* Double precision */
    return_val += test_bijectivity_suNg_field();
    return_val += test_bijectivity_suNf_field();
    return_val += test_bijectivity_suNg_scalar_field();
    return_val += test_bijectivity_suNg_av_field();
    return_val += test_bijectivity_gtransf();
    return_val += test_bijectivity_ldl_field();
    return_val += test_bijectivity_clover_term();
    return_val += test_bijectivity_clover_force();
    return_val += test_bijectivity_spinor_field(&glattice);
    return_val += test_bijectivity_scalar_field(&glattice);
    return_val += test_bijectivity_staple_field();

    /* Single precision */
    return_val += test_bijectivity_suNg_field_flt();
    return_val += test_bijectivity_suNf_field_flt();
    return_val += test_bijectivity_spinor_field_flt(&glattice);

    lprintf("INFO", 0, "\n\n Testing Even Lattice \n\n");
    return_val += test_bijectivity_spinor_field(&glat_even);
    return_val += test_bijectivity_spinor_field_flt(&glat_even);
    return_val += test_bijectivity_scalar_field(&glat_even);

    lprintf("INFO", 0, "\n\n Testing Odd Lattice \n\n");
    return_val += test_bijectivity_spinor_field(&glat_odd);
    return_val += test_bijectivity_spinor_field_flt(&glat_odd);
    return_val += test_bijectivity_scalar_field(&glat_odd);

    // Finalize and return
    finalize_process();
    return return_val;
}

int test_bijectivity_suNg_field() {
    lprintf("INFO", 0, " ====== TEST GAUGE FIELD ======= ");
    int return_val = 0;
    suNg_field *in, *in_copy;
    in = alloc_suNg_field(&glattice);
    in_copy = alloc_suNg_field(&glattice);

    random_suNg_field_cpu(in);

    copy_suNg_field_cpu(in_copy, in);
    lprintf("SANITY CHECK", 0, "CPU sqnorm: %0.2e\n", sqnorm_suNg_field_cpu(in));
    lprintf("SANITY CHECK", 0, "CPU copy sqnorm (should be the same as CPU sqnorm): %0.2e\n", sqnorm_suNg_field_cpu(in_copy));

    copy_to_gpu_suNg_field(in);
    zero_suNg_field_cpu(in);
    fill_buffers_with_zeroes_suNg_field(in);
    lprintf("SANITY CHECK", 0, "CPU copy should be zero in intermediate step: %0.2e\n", sqnorm_suNg_field_cpu(in));
    copy_from_gpu_suNg_field(in);

    sub_assign_suNg_field_cpu(in, in_copy);
    double diff_norm = sqnorm_suNg_field_cpu(in);
    return_val += check_diff_norm_zero(diff_norm);

    free_suNg_field(in);
    free_suNg_field(in_copy);
    return return_val;
}

int test_bijectivity_suNf_field() {
    lprintf("INFO", 0, " ====== TEST GAUGE FIELD ======= ");
    int return_val = 0;
    suNf_field *in, *in_copy;
    in = alloc_suNf_field(&glattice);
    in_copy = alloc_suNf_field(&glattice);

    random_suNf_field_cpu(in);
    random_suNf_field_cpu(in_copy);

    copy_suNf_field_cpu(in_copy, in);
    lprintf("SANITY CHECK", 0, "CPU sqnorm: %0.2e\n", sqnorm_suNf_field_cpu(in));
    lprintf("SANITY CHECK", 0, "CPU copy sqnorm (should be the same as CPU sqnorm): %0.2e\n", sqnorm_suNf_field_cpu(in_copy));

    copy_to_gpu_suNf_field(in);
    zero_suNf_field_cpu(in);
    fill_buffers_with_zeroes_suNf_field(in);
    lprintf("SANITY CHECK", 0, "CPU copy should be zero in intermediate step: %0.2e\n", sqnorm_suNf_field_cpu(in));
    copy_from_gpu_suNf_field(in);

    sub_assign_suNf_field_cpu(in, in_copy);
    double diff_norm = sqnorm_suNf_field_cpu(in);
    return_val += check_diff_norm_zero(diff_norm);

    free_suNf_field(in);
    free_suNf_field(in_copy);
    return return_val;
}

int test_bijectivity_spinor_field(geometry_descriptor *gd) {
    lprintf("INFO", 0, " ====== TEST SPINOR FIELD ======= ");
    int return_val = 0;
    spinor_field *in, *in_copy;
    in = alloc_spinor_field(1, gd);
    in_copy = alloc_spinor_field(1, gd);

    gaussian_spinor_field(in);

    spinor_field_copy_f_cpu(in_copy, in);
    lprintf("SANITY CHECK", 0, "CPU sqnorm: %0.2e\n", spinor_field_sqnorm_f_cpu(in));
    lprintf("SANITY CHECK", 0, "CPU copy sqnorm (should be the same as CPU sqnorm): %0.2e\n",
            spinor_field_sqnorm_f_cpu(in_copy));

    copy_to_gpu_spinor_field(in);
    fill_buffers_with_zeroes_spinor_field(in);
    spinor_field_zero_f_cpu(in);
    lprintf("SANITY CHECK", 0, "GPU copy non-zero in intermediate step: %0.2e\n", spinor_field_sqnorm_f(in));
    lprintf("SANITY CHECK", 0, "CPU copy should be zero in intermediate step: %0.2e\n", spinor_field_sqnorm_f_cpu(in));
    copy_from_gpu_spinor_field(in);

    lprintf("SANITY CHECK", 0, "CPU after copying back: %0.2e\n", spinor_field_sqnorm_f_cpu(in));

    spinor_field_sub_assign_f_cpu(in, in_copy);
    double diff_norm = spinor_field_sqnorm_f_cpu(in);
    return_val += check_diff_norm_zero(diff_norm);

    free_spinor_field(in);
    free_spinor_field(in_copy);
    return return_val;
}

int test_bijectivity_spinor_field_flt(geometry_descriptor *gd) {
    lprintf("INFO", 0, " ====== TEST SPINOR FIELD SINGLE PRECISION ======= ");
    int return_val = 0;
    spinor_field_flt *in, *in_copy;
    in = alloc_spinor_field_flt(1, gd);
    in_copy = alloc_spinor_field_flt(1, gd);

    gaussian_spinor_field_flt(in);

    spinor_field_copy_f_flt_cpu(in_copy, in);
    lprintf("SANITY CHECK", 0, "CPU sqnorm: %0.2e\n", spinor_field_sqnorm_f_flt_cpu(in));
    lprintf("SANITY CHECK", 0, "CPU copy sqnorm (should be the same as CPU sqnorm): %0.2e\n",
            spinor_field_sqnorm_f_flt_cpu(in_copy));

    copy_to_gpu_spinor_field_flt(in);
    spinor_field_zero_f_flt_cpu(in);
    fill_buffers_with_zeroes_spinor_field_flt(in);
    lprintf("SANITY CHECK", 0, "GPU copy non-zero in intermediate step: %0.2e\n", spinor_field_sqnorm_f_flt(in));
    lprintf("SANITY CHECK", 0, "CPU copy should be zero in intermediate step: %0.2e\n", spinor_field_sqnorm_f_flt_cpu(in));
    copy_from_gpu_spinor_field_flt(in);

    spinor_field_sub_assign_f_flt_cpu(in, in_copy);
    double diff_norm = spinor_field_sqnorm_f_flt_cpu(in);
    return_val += check_diff_norm_zero(diff_norm);

    free_spinor_field_flt(in);
    free_spinor_field_flt(in_copy);
    return return_val;
}

int test_bijectivity_scalar_field(geometry_descriptor *gd) {
    lprintf("INFO", 0, " ====== TEST SCALAR FIELD ======= ");
    int return_val = 0;
    scalar_field *in, *in_copy;
    in = alloc_scalar_field(1, gd);
    in_copy = alloc_scalar_field(1, gd);

    random_scalar_field_cpu(in);
    copy_scalar_field_cpu(in_copy, in);
    lprintf("SANITY CHECK", 0, "CPU sqnorm: %0.2e\n", sqnorm_scalar_field_cpu(in));
    lprintf("SANITY CHECK", 0, "CPU copy sqnorm (should be the same as CPU sqnorm): %0.2e\n", sqnorm_scalar_field_cpu(in_copy));

    copy_to_gpu_scalar_field(in);
    fill_buffers_with_zeroes_scalar_field(in);
    zero_scalar_field_cpu(in);
    lprintf("SANITY CHECK", 0, "CPU copy should be zero in intermediate step: %0.2e\n", sqnorm_scalar_field_cpu(in));
    copy_from_gpu_scalar_field(in);

    sub_assign_scalar_field_cpu(in, in_copy);
    double diff_norm = sqnorm_scalar_field_cpu(in);
    return_val += check_diff_norm_zero(diff_norm);

    free_scalar_field(in);
    free_scalar_field(in_copy);
    return return_val;
}

int test_bijectivity_suNg_field_flt() {
    lprintf("INFO", 0, " ====== TEST GAUGE FIELD SINGLE PRECISION ======= ");
    int return_val = 0;
    suNg_field_flt *in, *in_copy;
    in = alloc_suNg_field_flt(&glattice);
    in_copy = alloc_suNg_field_flt(&glattice);

    random_suNg_field_flt_cpu(in);

    copy_suNg_field_flt_cpu(in_copy, in);
    lprintf("SANITY CHECK", 0, "CPU sqnorm: %0.2e\n", sqnorm_suNg_field_flt_cpu(in));
    lprintf("SANITY CHECK", 0, "CPU copy sqnorm (should be the same as CPU sqnorm): %0.2e\n",
            sqnorm_suNg_field_flt_cpu(in_copy));

    copy_to_gpu_suNg_field_flt(in);
    zero_suNg_field_flt_cpu(in);
    fill_buffers_with_zeroes_suNg_field_flt(in);
    lprintf("SANITY CHECK", 0, "CPU copy should be zero in intermediate step: %0.2e\n", sqnorm_suNg_field_flt_cpu(in));
    copy_from_gpu_suNg_field_flt(in);

    sub_assign_suNg_field_flt_cpu(in, in_copy);
    double diff_norm = sqnorm_suNg_field_flt_cpu(in);
    return_val += check_diff_norm_zero(diff_norm);

    free_suNg_field_flt(in);
    free_suNg_field_flt(in_copy);
    return return_val;
}

int test_bijectivity_suNf_field_flt() {
    lprintf("INFO", 0, " ====== TEST GAUGE FIELD REPRESENTED SINGLE PRECISION ======= ");
    int return_val = 0;
    suNf_field_flt *in, *in_copy;
    in = alloc_suNf_field_flt(&glattice);
    in_copy = alloc_suNf_field_flt(&glattice);

    random_suNf_field_flt_cpu(in);

    copy_suNf_field_flt_cpu(in_copy, in);
    lprintf("SANITY CHECK", 0, "CPU sqnorm: %0.2e\n", sqnorm_suNf_field_flt_cpu(in));
    lprintf("SANITY CHECK", 0, "CPU copy sqnorm (should be the same as CPU sqnorm): %0.2e\n",
            sqnorm_suNf_field_flt_cpu(in_copy));

    copy_to_gpu_suNf_field_flt(in);
    zero_suNf_field_flt_cpu(in);
    fill_buffers_with_zeroes_suNf_field_flt(in);
    lprintf("SANITY CHECK", 0, "CPU copy should be zero in intermediate step: %0.2e\n", sqnorm_suNf_field_flt_cpu(in));
    copy_from_gpu_suNf_field_flt(in);

    sub_assign_suNf_field_flt_cpu(in, in_copy);
    double diff_norm = sqnorm_suNf_field_flt_cpu(in);
    return_val += check_diff_norm_zero(diff_norm);

    free_suNf_field_flt(in);
    free_suNf_field_flt(in_copy);
    return return_val;
}

int test_bijectivity_suNg_scalar_field() {
    lprintf("INFO", 0, " ====== TEST SU(NG) SCALAR FIELD ======= ");
    int return_val = 0;
    suNg_scalar_field *in, *in_copy;
    in = alloc_suNg_scalar_field(&glattice);
    in_copy = alloc_suNg_scalar_field(&glattice);

    random_suNg_scalar_field_cpu(in);

    copy_suNg_scalar_field_cpu(in_copy, in);
    lprintf("SANITY CHECK", 0, "CPU sqnorm: %0.2e\n", sqnorm_suNg_scalar_field_cpu(in));
    lprintf("SANITY CHECK", 0, "CPU copy sqnorm (should be the same as CPU sqnorm): %0.2e\n",
            sqnorm_suNg_scalar_field_cpu(in_copy));

    copy_to_gpu_suNg_scalar_field(in);
    zero_suNg_scalar_field_cpu(in);
    fill_buffers_with_zeroes_suNg_scalar_field(in);
    lprintf("SANITY CHECK", 0, "CPU copy should be zero in intermediate step: %0.2e\n", sqnorm_suNg_scalar_field_cpu(in));
    copy_from_gpu_suNg_scalar_field(in);

    sub_assign_suNg_scalar_field_cpu(in, in_copy);
    double diff_norm = sqnorm_suNg_scalar_field_cpu(in);
    return_val += check_diff_norm_zero(diff_norm);

    free_suNg_scalar_field(in);
    free_suNg_scalar_field(in_copy);
    return return_val;
}

int test_bijectivity_suNg_av_field() {
    lprintf("INFO", 0, " ====== TEST ALGEBRA VECTOR FIELD ======= ");
    int return_val = 0;
    suNg_av_field *in, *in_copy;
    in = alloc_suNg_av_field(&glattice);
    in_copy = alloc_suNg_av_field(&glattice);

    random_suNg_av_field_cpu(in);

    copy_suNg_av_field_cpu(in_copy, in);
    lprintf("SANITY CHECK", 0, "CPU sqnorm: %0.2e\n", sqnorm_suNg_av_field_cpu(in));
    lprintf("SANITY CHECK", 0, "CPU copy sqnorm (should be the same as CPU sqnorm): %0.2e\n",
            sqnorm_suNg_av_field_cpu(in_copy));

    copy_to_gpu_suNg_av_field(in);
    zero_suNg_av_field_cpu(in);
    fill_buffers_with_zeroes_suNg_av_field(in);
    lprintf("SANITY CHECK", 0, "CPU copy should be zero in intermediate step: %0.2e\n", sqnorm_suNg_av_field_cpu(in));
    copy_from_gpu_suNg_av_field(in);

    sub_assign_suNg_av_field_cpu(in, in_copy);
    double diff_norm = sqnorm_suNg_av_field_cpu(in);
    return_val += check_diff_norm_zero(diff_norm);

    free_suNg_av_field(in);
    free_suNg_av_field(in_copy);
    return return_val;
}

int test_bijectivity_gtransf() {
    lprintf("INFO", 0, " ====== TEST GAUGE TRANSFORMATION ======= ");
    int return_val = 0;
    gtransf *in, *in_copy;
    in = alloc_gtransf(&glattice);
    in_copy = alloc_gtransf(&glattice);

    random_gtransf_cpu(in);

    copy_gtransf_cpu(in_copy, in);
    lprintf("SANITY CHECK", 0, "CPU sqnorm: %0.2e\n", sqnorm_gtransf_cpu(in));
    lprintf("SANITY CHECK", 0, "CPU copy sqnorm (should be the same as CPU sqnorm): %0.2e\n", sqnorm_gtransf_cpu(in_copy));

    copy_to_gpu_gtransf(in);
    zero_gtransf_cpu(in);
    fill_buffers_with_zeroes_gtransf(in);
    lprintf("SANITY CHECK", 0, "CPU copy should be zero in intermediate step: %0.2e\n", sqnorm_gtransf_cpu(in));
    copy_from_gpu_gtransf(in);

    sub_assign_gtransf_cpu(in, in_copy);
    double diff_norm = sqnorm_gtransf_cpu(in);
    return_val += check_diff_norm_zero(diff_norm);

    free_gtransf(in);
    free_gtransf(in_copy);
    return return_val;
}

int test_bijectivity_clover_term() {
    lprintf("INFO", 0, " ====== TEST CLOVER TERM ======= ");
    int return_val = 0;
    clover_term *in, *in_copy;
    in = alloc_clover_term(&glattice);
    in_copy = alloc_clover_term(&glattice);

    random_clover_term_cpu(in);

    copy_clover_term_cpu(in_copy, in);
    lprintf("SANITY CHECK", 0, "CPU sqnorm: %0.2e\n", sqnorm_clover_term_cpu(in));
    lprintf("SANITY CHECK", 0, "CPU copy sqnorm (should be the same as CPU sqnorm): %0.2e\n", sqnorm_clover_term_cpu(in_copy));

    copy_to_gpu_clover_term(in);
    zero_clover_term_cpu(in);
    fill_buffers_with_zeroes_clover_term(in);
    lprintf("SANITY CHECK", 0, "CPU copy should be zero in intermediate step: %0.2e\n", sqnorm_clover_term_cpu(in));
    copy_from_gpu_clover_term(in);

    sub_assign_clover_term_cpu(in, in_copy);
    double diff_norm = sqnorm_clover_term_cpu(in);
    return_val += check_diff_norm_zero(diff_norm);

    free_clover_term(in);
    free_clover_term(in_copy);
    return return_val;
}

int test_bijectivity_clover_force() {
    lprintf("INFO", 0, " ====== TEST CLOVER FORCE ======= ");
    int return_val = 0;
    clover_force *in, *in_copy;
    in = alloc_clover_force(&glattice);
    in_copy = alloc_clover_force(&glattice);

    random_clover_force_cpu(in);

    copy_clover_force_cpu(in_copy, in);
    lprintf("SANITY CHECK", 0, "CPU sqnorm: %0.2e\n", sqnorm_clover_force_cpu(in));
    lprintf("SANITY CHECK", 0, "CPU copy sqnorm (should be the same as CPU sqnorm): %0.2e\n", sqnorm_clover_force_cpu(in_copy));

    copy_to_gpu_clover_force(in);
    zero_clover_force_cpu(in);
    fill_buffers_with_zeroes_clover_force(in);
    lprintf("SANITY CHECK", 0, "CPU copy should be zero in intermediate step: %0.2e\n", sqnorm_clover_force_cpu(in));
    copy_from_gpu_clover_force(in);

    sub_assign_clover_force_cpu(in, in_copy);
    double diff_norm = sqnorm_clover_force_cpu(in);
    return_val += check_diff_norm_zero(diff_norm);

    free_clover_force(in);
    free_clover_force(in_copy);
    return return_val;
}

int test_bijectivity_ldl_field() {
    lprintf("INFO", 0, " ====== TEST CLOVER LDL ======= ");
    int return_val = 0;
    ldl_field *in, *in_copy;
    in = alloc_ldl_field(&glattice);
    in_copy = alloc_ldl_field(&glattice);

    random_ldl_field_cpu(in);

    copy_ldl_field_cpu(in_copy, in);
    lprintf("SANITY CHECK", 0, "CPU sqnorm: %0.2e\n", sqnorm_ldl_field_cpu(in));
    lprintf("SANITY CHECK", 0, "CPU copy sqnorm (should be the same as CPU sqnorm): %0.2e\n", sqnorm_ldl_field_cpu(in_copy));

    copy_to_gpu_ldl_field(in);
    zero_ldl_field_cpu(in);
    fill_buffers_with_zeroes_ldl_field(in);
    lprintf("SANITY CHECK", 0, "CPU copy should be zero in intermediate step: %0.2e\n", sqnorm_ldl_field_cpu(in));
    copy_from_gpu_ldl_field(in);

    sub_assign_ldl_field_cpu(in, in_copy);
    double diff_norm = sqnorm_ldl_field_cpu(in);
    return_val += check_diff_norm_zero(diff_norm);

    free_ldl_field(in);
    free_ldl_field(in_copy);
    return return_val;
}

int test_bijectivity_staple_field() {
    lprintf("INFO", 0, " ====== TEST STAPLE FIELD ======= ");
    int return_val = 0;
    staple_field *in, *in_copy;
    in = alloc_staple_field(&glattice);
    in_copy = alloc_staple_field(&glattice);

    random_staple_field_cpu(in);

    copy_staple_field_cpu(in_copy, in);
    lprintf("SANITY CHECK", 0, "CPU sqnorm: %0.2e\n", sqnorm_staple_field_cpu(in));
    lprintf("SANITY CHECK", 0, "CPU copy sqnorm (should be the same as CPU sqnorm): %0.2e\n", sqnorm_staple_field_cpu(in_copy));

    copy_to_gpu_staple_field(in);
    zero_staple_field_cpu(in);
    fill_buffers_with_zeroes_staple_field(in);
    lprintf("SANITY CHECK", 0, "CPU copy should be zero in intermediate step: %0.2e\n", sqnorm_staple_field_cpu(in));
    copy_from_gpu_staple_field(in);

    sub_assign_staple_field_cpu(in, in_copy);
    double diff_norm = sqnorm_staple_field_cpu(in);
    return_val += check_diff_norm_zero(diff_norm);

    free_staple_field(in);
    free_staple_field(in_copy);
    return return_val;
}
