/*******************************************************************************
*
* NOCOMPILE= !WITH_GPU
*
* Check that copy back and forth works for all fields.
*
*******************************************************************************/

#include "libhr.h"

// Double precision
int test_bijectivity_gfield();
int test_bijectivity_gfield_f();
int test_bijectivity_suNg_scalar_field();
int test_bijectivity_avfield();
int test_bijectivity_gtransf();
int test_bijectivity_clover_ldl();
int test_bijectivity_clover_term();
int test_bijectivity_clover_force();
int test_bijectivity_spinor_field_f(geometry_descriptor *);
int test_bijectivity_sfield(geometry_descriptor *);
int test_bijectivity_staple_field();

// Single precision
int test_bijectivity_gfield_flt();
int test_bijectivity_gfield_f_flt();
int test_bijectivity_spinor_field_f_flt(geometry_descriptor *);

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
    return_val += test_bijectivity_gfield();
    return_val += test_bijectivity_gfield_f();
    return_val += test_bijectivity_suNg_scalar_field();
    return_val += test_bijectivity_avfield();
    return_val += test_bijectivity_gtransf();
    return_val += test_bijectivity_clover_ldl();
    return_val += test_bijectivity_clover_term();
    return_val += test_bijectivity_clover_force();
    return_val += test_bijectivity_spinor_field_f(&glattice);
    return_val += test_bijectivity_sfield(&glattice);
    return_val += test_bijectivity_staple_field();

    /* Single precision */
    return_val += test_bijectivity_gfield_flt();
    return_val += test_bijectivity_gfield_f_flt();
    return_val += test_bijectivity_spinor_field_f_flt(&glattice);

    lprintf("INFO", 0, "\n\n Testing Even Lattice \n\n");
    return_val += test_bijectivity_spinor_field_f(&glat_even);
    return_val += test_bijectivity_spinor_field_f_flt(&glat_even);
    return_val += test_bijectivity_sfield(&glat_even);

    lprintf("INFO", 0, "\n\n Testing Odd Lattice \n\n");
    return_val += test_bijectivity_spinor_field_f(&glat_odd);
    return_val += test_bijectivity_spinor_field_f_flt(&glat_odd);
    return_val += test_bijectivity_sfield(&glat_odd);

    // Finalize and return
    finalize_process();
    return return_val;
}

int test_bijectivity_gfield() {
    lprintf("INFO", 0, " ====== TEST GAUGE FIELD ======= ");
    int return_val = 0;
    suNg_field *in, *in_copy;
    in = alloc_gfield(&glattice);
    in_copy = alloc_gfield(&glattice);

    random_gfield_cpu(in);

    copy_gfield_cpu(in_copy, in);
    lprintf("SANITY CHECK", 0, "CPU sqnorm: %0.2e\n", sqnorm_gfield_cpu(in));
    lprintf("SANITY CHECK", 0, "CPU copy sqnorm (should be the same as CPU sqnorm): %0.2e\n", sqnorm_gfield_cpu(in_copy));

    copy_to_gpu_gfield(in);
    zero_gfield_cpu(in);
    fill_buffers_with_zeroes_gfield(in);
    lprintf("SANITY CHECK", 0, "CPU copy should be zero in intermediate step: %0.2e\n", sqnorm_gfield_cpu(in));
    copy_from_gpu_gfield(in);

    sub_assign_gfield_cpu(in, in_copy);
    double diff_norm = sqnorm_gfield_cpu(in);
    return_val += check_diff_norm_zero(diff_norm);

    free_gfield(in);
    free_gfield(in_copy);
    return return_val;
}

int test_bijectivity_gfield_f() {
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
    zero_gfield_f_cpu(in);
    fill_buffers_with_zeroes_gfield_f(in);
    lprintf("SANITY CHECK", 0, "CPU copy should be zero in intermediate step: %0.2e\n", sqnorm_gfield_f_cpu(in));
    copy_from_gpu_gfield_f(in);

    sub_assign_gfield_f_cpu(in, in_copy);
    double diff_norm = sqnorm_gfield_f_cpu(in);
    return_val += check_diff_norm_zero(diff_norm);

    free_gfield_f(in);
    free_gfield_f(in_copy);
    return return_val;
}

int test_bijectivity_spinor_field_f(geometry_descriptor *gd) {
    lprintf("INFO", 0, " ====== TEST SPINOR FIELD ======= ");
    int return_val = 0;
    spinor_field *in, *in_copy;
    in = alloc_spinor_field_f(1, gd);
    in_copy = alloc_spinor_field_f(1, gd);

    gaussian_spinor_field(in);

    spinor_field_copy_f_cpu(in_copy, in);
    lprintf("SANITY CHECK", 0, "CPU sqnorm: %0.2e\n", spinor_field_sqnorm_f_cpu(in));
    lprintf("SANITY CHECK", 0, "CPU copy sqnorm (should be the same as CPU sqnorm): %0.2e\n",
            spinor_field_sqnorm_f_cpu(in_copy));

    copy_to_gpu_spinor_field_f(in);
    fill_buffers_with_zeroes_spinor_field_f(in);
    spinor_field_zero_f_cpu(in);
    lprintf("SANITY CHECK", 0, "GPU copy non-zero in intermediate step: %0.2e\n", spinor_field_sqnorm_f(in));
    lprintf("SANITY CHECK", 0, "CPU copy should be zero in intermediate step: %0.2e\n", spinor_field_sqnorm_f_cpu(in));
    copy_from_gpu_spinor_field_f(in);

    lprintf("SANITY CHECK", 0, "CPU after copying back: %0.2e\n", spinor_field_sqnorm_f_cpu(in));

    spinor_field_sub_assign_f_cpu(in, in_copy);
    double diff_norm = spinor_field_sqnorm_f_cpu(in);
    return_val += check_diff_norm_zero(diff_norm);

    free_spinor_field_f(in);
    free_spinor_field_f(in_copy);
    return return_val;
}

int test_bijectivity_spinor_field_f_flt(geometry_descriptor *gd) {
    lprintf("INFO", 0, " ====== TEST SPINOR FIELD SINGLE PRECISION ======= ");
    int return_val = 0;
    spinor_field_flt *in, *in_copy;
    in = alloc_spinor_field_f_flt(1, gd);
    in_copy = alloc_spinor_field_f_flt(1, gd);

    gaussian_spinor_field_flt(in);

    spinor_field_copy_f_flt_cpu(in_copy, in);
    lprintf("SANITY CHECK", 0, "CPU sqnorm: %0.2e\n", spinor_field_sqnorm_f_flt_cpu(in));
    lprintf("SANITY CHECK", 0, "CPU copy sqnorm (should be the same as CPU sqnorm): %0.2e\n",
            spinor_field_sqnorm_f_flt_cpu(in_copy));

    copy_to_gpu_spinor_field_f_flt(in);
    spinor_field_zero_f_flt_cpu(in);
    fill_buffers_with_zeroes_spinor_field_f_flt(in);
    lprintf("SANITY CHECK", 0, "GPU copy non-zero in intermediate step: %0.2e\n", spinor_field_sqnorm_f_flt(in));
    lprintf("SANITY CHECK", 0, "CPU copy should be zero in intermediate step: %0.2e\n", spinor_field_sqnorm_f_flt_cpu(in));
    copy_from_gpu_spinor_field_f_flt(in);

    spinor_field_sub_assign_f_flt_cpu(in, in_copy);
    double diff_norm = spinor_field_sqnorm_f_flt_cpu(in);
    return_val += check_diff_norm_zero(diff_norm);

    free_spinor_field_f_flt(in);
    free_spinor_field_f_flt(in_copy);
    return return_val;
}

int test_bijectivity_sfield(geometry_descriptor *gd) {
    lprintf("INFO", 0, " ====== TEST SCALAR FIELD ======= ");
    int return_val = 0;
    scalar_field *in, *in_copy;
    in = alloc_sfield(1, gd);
    in_copy = alloc_sfield(1, gd);

    random_sfield_cpu(in);
    copy_sfield_cpu(in_copy, in);
    lprintf("SANITY CHECK", 0, "CPU sqnorm: %0.2e\n", sqnorm_sfield_cpu(in));
    lprintf("SANITY CHECK", 0, "CPU copy sqnorm (should be the same as CPU sqnorm): %0.2e\n", sqnorm_sfield_cpu(in_copy));

    copy_to_gpu_sfield(in);
    fill_buffers_with_zeroes_sfield(in);
    zero_sfield_cpu(in);
    lprintf("SANITY CHECK", 0, "CPU copy should be zero in intermediate step: %0.2e\n", sqnorm_sfield_cpu(in));
    copy_from_gpu_sfield(in);

    sub_assign_sfield_cpu(in, in_copy);
    double diff_norm = sqnorm_sfield_cpu(in);
    return_val += check_diff_norm_zero(diff_norm);

    free_sfield(in);
    free_sfield(in_copy);
    return return_val;
}

int test_bijectivity_gfield_flt() {
    lprintf("INFO", 0, " ====== TEST GAUGE FIELD SINGLE PRECISION ======= ");
    int return_val = 0;
    suNg_field_flt *in, *in_copy;
    in = alloc_gfield_flt(&glattice);
    in_copy = alloc_gfield_flt(&glattice);

    random_gfield_flt_cpu(in);

    copy_gfield_flt_cpu(in_copy, in);
    lprintf("SANITY CHECK", 0, "CPU sqnorm: %0.2e\n", sqnorm_gfield_flt_cpu(in));
    lprintf("SANITY CHECK", 0, "CPU copy sqnorm (should be the same as CPU sqnorm): %0.2e\n", sqnorm_gfield_flt_cpu(in_copy));

    copy_to_gpu_gfield_flt(in);
    zero_gfield_flt_cpu(in);
    fill_buffers_with_zeroes_gfield_flt(in);
    lprintf("SANITY CHECK", 0, "CPU copy should be zero in intermediate step: %0.2e\n", sqnorm_gfield_flt_cpu(in));
    copy_from_gpu_gfield_flt(in);

    sub_assign_gfield_flt_cpu(in, in_copy);
    double diff_norm = sqnorm_gfield_flt_cpu(in);
    return_val += check_diff_norm_zero(diff_norm);

    free_gfield_flt(in);
    free_gfield_flt(in_copy);
    return return_val;
}

int test_bijectivity_gfield_f_flt() {
    lprintf("INFO", 0, " ====== TEST GAUGE FIELD REPRESENTED SINGLE PRECISION ======= ");
    int return_val = 0;
    suNf_field_flt *in, *in_copy;
    in = alloc_gfield_f_flt(&glattice);
    in_copy = alloc_gfield_f_flt(&glattice);

    random_gfield_f_flt_cpu(in);

    copy_gfield_f_flt_cpu(in_copy, in);
    lprintf("SANITY CHECK", 0, "CPU sqnorm: %0.2e\n", sqnorm_gfield_f_flt_cpu(in));
    lprintf("SANITY CHECK", 0, "CPU copy sqnorm (should be the same as CPU sqnorm): %0.2e\n", sqnorm_gfield_f_flt_cpu(in_copy));

    copy_to_gpu_gfield_f_flt(in);
    zero_gfield_f_flt_cpu(in);
    fill_buffers_with_zeroes_gfield_f_flt(in);
    lprintf("SANITY CHECK", 0, "CPU copy should be zero in intermediate step: %0.2e\n", sqnorm_gfield_f_flt_cpu(in));
    copy_from_gpu_gfield_f_flt(in);

    sub_assign_gfield_f_flt_cpu(in, in_copy);
    double diff_norm = sqnorm_gfield_f_flt_cpu(in);
    return_val += check_diff_norm_zero(diff_norm);

    free_gfield_f_flt(in);
    free_gfield_f_flt(in_copy);
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

int test_bijectivity_avfield() {
    lprintf("INFO", 0, " ====== TEST ALGEBRA VECTOR FIELD ======= ");
    int return_val = 0;
    suNg_av_field *in, *in_copy;
    in = alloc_avfield(&glattice);
    in_copy = alloc_avfield(&glattice);

    random_avfield_cpu(in);

    copy_avfield_cpu(in_copy, in);
    lprintf("SANITY CHECK", 0, "CPU sqnorm: %0.2e\n", sqnorm_avfield_cpu(in));
    lprintf("SANITY CHECK", 0, "CPU copy sqnorm (should be the same as CPU sqnorm): %0.2e\n", sqnorm_avfield_cpu(in_copy));

    copy_to_gpu_avfield(in);
    zero_avfield_cpu(in);
    fill_buffers_with_zeroes_avfield(in);
    lprintf("SANITY CHECK", 0, "CPU copy should be zero in intermediate step: %0.2e\n", sqnorm_avfield_cpu(in));
    copy_from_gpu_avfield(in);

    sub_assign_avfield_cpu(in, in_copy);
    double diff_norm = sqnorm_avfield_cpu(in);
    return_val += check_diff_norm_zero(diff_norm);

    free_avfield(in);
    free_avfield(in_copy);
    return return_val;
}

int test_bijectivity_gtransf() {
    lprintf("INFO", 0, " ====== TEST GAUGE TRANSFORMATION ======= ");
    int return_val = 0;
    suNg_field *in, *in_copy;
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
    suNfc_field *in, *in_copy;
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
    suNf_field *in, *in_copy;
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

int test_bijectivity_clover_ldl() {
    lprintf("INFO", 0, " ====== TEST CLOVER LDL ======= ");
    int return_val = 0;
    ldl_field *in, *in_copy;
    in = alloc_clover_ldl(&glattice);
    in_copy = alloc_clover_ldl(&glattice);

    random_clover_ldl_cpu(in);

    copy_clover_ldl_cpu(in_copy, in);
    lprintf("SANITY CHECK", 0, "CPU sqnorm: %0.2e\n", sqnorm_clover_ldl_cpu(in));
    lprintf("SANITY CHECK", 0, "CPU copy sqnorm (should be the same as CPU sqnorm): %0.2e\n", sqnorm_clover_ldl_cpu(in_copy));

    copy_to_gpu_clover_ldl(in);
    zero_clover_ldl_cpu(in);
    fill_buffers_with_zeroes_clover_ldl(in);
    lprintf("SANITY CHECK", 0, "CPU copy should be zero in intermediate step: %0.2e\n", sqnorm_clover_ldl_cpu(in));
    copy_from_gpu_clover_ldl(in);

    sub_assign_clover_ldl_cpu(in, in_copy);
    double diff_norm = sqnorm_clover_ldl_cpu(in);
    return_val += check_diff_norm_zero(diff_norm);

    free_clover_ldl(in);
    free_clover_ldl(in_copy);
    return return_val;
}

int test_bijectivity_staple_field() {
    lprintf("INFO", 0, " ====== TEST STAPLE FIELD ======= ");
    int return_val = 0;
    suNg_field *in, *in_copy;
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
