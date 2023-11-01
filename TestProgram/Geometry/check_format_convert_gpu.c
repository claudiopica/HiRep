/*******************************************************************************
*
* NOCOMPILE= !WITH_GPU
*
* Check that the conversion functions from GPU to CPU format and back 
* are bijective
*
*******************************************************************************/

#include "libhr.h"

/* Double precision tests */
int test_convert_back_forth_suNg_field();
int test_convert_back_forth_suNf_field();
int test_convert_back_forth_suNg_scalar_field();
int test_convert_back_forth_suNg_av_field();
int test_convert_back_forth_gtransf();
int test_convert_back_forth_ldl_field();
int test_convert_back_forth_clover_term();
int test_convert_back_forth_clover_force();
int test_convert_back_forth_spinor_field();
int test_convert_back_forth_scalar_field();
int test_convert_back_forth_staple_field();

/* Single precision tests */
int test_convert_back_forth_suNg_field_flt();
int test_convert_back_forth_suNf_field_flt();
int test_convert_back_forth_spinor_field_flt();

int main(int argc, char *argv[]) {
    // Init
    int return_val = 0;

    // Setup process and communication
    logger_map("DEBUG", "debug");
    setup_process(&argc, &argv);
    test_setup();

    // Run tests
    /* Double precision */
    return_val += test_convert_back_forth_suNg_field();
    return_val += test_convert_back_forth_suNf_field();
    return_val += test_convert_back_forth_suNg_scalar_field();
    return_val += test_convert_back_forth_suNg_av_field();
    return_val += test_convert_back_forth_gtransf();
    return_val += test_convert_back_forth_ldl_field();
    return_val += test_convert_back_forth_clover_term();
    return_val += test_convert_back_forth_clover_force();
    return_val += test_convert_back_forth_spinor_field();
    return_val += test_convert_back_forth_scalar_field();
    return_val += test_convert_back_forth_staple_field();

    /* Single precision */
    return_val += test_convert_back_forth_suNg_field_flt();
    return_val += test_convert_back_forth_suNf_field_flt();
    return_val += test_convert_back_forth_spinor_field_flt();

    // Finalize and return
    finalize_process();
    return return_val;
}

int test_convert_back_forth_spinor_field() {
    lprintf("INFO", 0, " ======= TEST SPINOR FIELD ======= \n");

    // Setup spinor fields
    spinor_field *in, *tmp, *out;
    in = alloc_spinor_field(1, &glattice);
    tmp = alloc_spinor_field(1, &glattice);
    out = alloc_spinor_field(1, &glattice);
    gaussian_spinor_field(in);

    lprintf("SANITY CHECK", 0, "[In field CPU copy norm unequal zero: %0.2e]\n", spinor_field_sqnorm_f_cpu(in));

    // Convert twice
    to_gpu_format_spinor_field(tmp, in);
    fill_buffers_spinor_field(tmp);
    to_cpu_format_spinor_field(out, tmp);

    lprintf("SANITY CHECK", 0, "[In and outfield sqnorms unequal zero and equal to each other: in %0.2e out %0.2e]\n",
            spinor_field_sqnorm_f_cpu(in), spinor_field_sqnorm_f_cpu(out));

    // Assert fields are equal over sqnorm
    spinor_field_sub_assign_f_cpu(out, in);
    double diff_norm = spinor_field_sqnorm_f_cpu(out);

    // Free and return
    free_spinor_field(in);
    free_spinor_field(tmp);
    free_spinor_field(out);
    return check_diff_norm_zero(diff_norm);
}

int test_convert_back_forth_spinor_field_flt() {
    lprintf("INFO", 0, " ======= TEST SPINOR FIELD SINGLE PRECISION ======= \n");

    // Setup spinor fields
    spinor_field_flt *in, *tmp, *out;
    in = alloc_spinor_field_flt(1, &glattice);
    tmp = alloc_spinor_field_flt(1, &glattice);
    out = alloc_spinor_field_flt(1, &glattice);
    gaussian_spinor_field_flt(in);

    double sqnorm = spinor_field_sqnorm_f_flt_cpu(in);
    lprintf("SANITY CHECK", 0, "[In field CPU copy norm unequal zero: %0.2e]\n", sqnorm);

    // Convert twice
    to_gpu_format_spinor_field_flt(tmp, in);
    fill_buffers_spinor_field_flt(tmp);
    to_cpu_format_spinor_field_flt(out, tmp);
    lprintf("SANITY CHECK", 0, "[In and outfield sqnorms unequal zero and equal to each other: in: %0.2e out: %0.2e]\n",
            spinor_field_sqnorm_f_flt_cpu(in), spinor_field_sqnorm_f_flt_cpu(out));

    // Assert fields are equal over sqnorm
    spinor_field_sub_assign_f_flt_cpu(out, in);
    double diff_norm = spinor_field_sqnorm_f_flt_cpu(out);

    // Free and return
    free_spinor_field_flt(in);
    free_spinor_field_flt(tmp);
    free_spinor_field_flt(out);
    return check_diff_norm_zero(diff_norm);
}

int test_convert_back_forth_suNf_field() {
    lprintf("INFO", 0, " ======= TEST GAUGE FIELD REPRESENTED ======= \n");

    // Setup suNg_fields
    suNf_field *in, *tmp, *out;
    in = alloc_suNf_field(&glattice);

    tmp = alloc_suNf_field(&glattice);
    out = alloc_suNf_field(&glattice);

    random_u_f(in);
    lprintf("SANITY CHECK", 0, "[In field CPU copy norm unequal zero: %0.2e]\n", sqnorm_suNf_field_cpu(in));

    // Convert twice
    to_gpu_format_suNf_field(tmp, in);
    fill_buffers_suNf_field(tmp);
    to_cpu_format_suNf_field(out, tmp);

    lprintf("SANITY CHECK", 0, "[In and outfield sqnorms unequal zero and equal to each other: in %0.2e out %0.2e]\n",
            sqnorm_suNf_field_cpu(in), sqnorm_suNf_field_cpu(out));

    // Assert fields are equal over sqnorm
    sub_assign_suNf_field_cpu(out, in);
    double diff_norm = sqnorm_suNf_field_cpu(out);
    return check_diff_norm_zero(diff_norm);

    free_suNf_field(in);
    free_suNf_field(tmp);
    free_suNf_field(out);
    return 0;
}

int test_convert_back_forth_suNf_field_flt() {
    lprintf("INFO", 0, " ======= TEST GAUGE FIELD REPRESENTED SINGLE PRECISION ======= \n");

    // Setup suNg_fields
    suNf_field_flt *in, *tmp, *out;
    in = alloc_suNf_field_flt(&glattice);
    tmp = alloc_suNf_field_flt(&glattice);
    out = alloc_suNf_field_flt(&glattice);

    random_suNf_field_flt_cpu(in);
    lprintf("SANITY CHECK", 0, "[In field CPU copy norm unequal zero: %0.2e]\n", sqnorm_suNf_field_flt_cpu(in));

    // Convert twice
    to_gpu_format_suNf_field_flt(tmp, in);
    fill_buffers_suNf_field_flt(tmp);
    to_cpu_format_suNf_field_flt(out, tmp);
    lprintf("SANITY CHECK", 0, "[In and outfield sqnorms unequal zero and equal to each other: in %0.2e out %0.2e]\n",
            sqnorm_suNf_field_flt_cpu(in), sqnorm_suNf_field_flt_cpu(out));

    // Assert fields are equal over sqnorm
    sub_assign_suNf_field_flt_cpu(out, in);
    double diff_norm = sqnorm_suNf_field_flt_cpu(out);

    // Free and return
    free_suNf_field_flt(in);
    free_suNf_field_flt(tmp);
    free_suNf_field_flt(out);
    return check_diff_norm_zero(diff_norm);
}

int test_convert_back_forth_suNg_field() {
    lprintf("INFO", 0, " ======= TEST GAUGE FIELD ======= \n");

    // Setup suNg_fields
    suNg_field *in, *tmp, *out;
    in = alloc_suNg_field(&glattice);
    tmp = alloc_suNg_field(&glattice);
    out = alloc_suNg_field(&glattice);

    random_u(in);
    lprintf("SANITY CHECK", 0, "[In field CPU copy norm unequal zero: %0.2e]\n", sqnorm_suNg_field_cpu(in));

    // Convert twice
    to_gpu_format_suNg_field(tmp, in);
    fill_buffers_suNg_field(tmp);
    to_cpu_format_suNg_field(out, tmp);
    lprintf("SANITY CHECK", 0, "[In and outfield sqnorms unequal zero and equal to each other: in %0.2e out %0.2e]\n",
            sqnorm_suNg_field_cpu(in), sqnorm_suNg_field_cpu(out));

    // Assert fields are equal over sqnorm
    sub_assign_suNg_field_cpu(out, in);
    double diff_norm = sqnorm_suNg_field_cpu(out);

    // Free and return
    free_suNg_field(in);
    free_suNg_field(tmp);
    free_suNg_field(out);
    return check_diff_norm_zero(diff_norm);
}

int test_convert_back_forth_suNg_field_flt() {
    lprintf("INFO", 0, " ======= TEST GAUGE FIELD SINGLE PRECISION ======= \n");

    // Setup suNg_fields
    suNg_field_flt *in, *tmp, *out;
    in = alloc_suNg_field_flt(&glattice);
    tmp = alloc_suNg_field_flt(&glattice);
    out = alloc_suNg_field_flt(&glattice);

    random_suNg_field_flt_cpu(in);
    lprintf("SANITY CHECK", 0, "[In field CPU copy norm unequal zero: %0.2e]\n", sqnorm_suNg_field_flt_cpu(in));

    // Convert twice
    to_gpu_format_suNg_field_flt(tmp, in);
    fill_buffers_suNg_field_flt(tmp);
    to_cpu_format_suNg_field_flt(out, tmp);
    lprintf("SANITY CHECK", 0, "[In and outfield sqnorms unequal zero and equal to each other: in %0.2e out %0.2e]\n",
            sqnorm_suNg_field_flt_cpu(in), sqnorm_suNg_field_flt_cpu(out));

    // Assert fields are equal over sqnorm
    sub_assign_suNg_field_flt_cpu(out, in);
    double diff_norm = sqnorm_suNg_field_flt_cpu(out);

    free_suNg_field_flt(in);
    free_suNg_field_flt(tmp);
    free_suNg_field_flt(out);
    return check_diff_norm_zero(diff_norm);
}

int test_convert_back_forth_suNg_scalar_field() {
    lprintf("INFO", 0, " ======= TEST SU(NG) SCALAR FIELD ======= \n");

    // Setup scalar fields
    suNg_scalar_field *in, *tmp, *out;
    in = alloc_suNg_scalar_field(&glattice);
    tmp = alloc_suNg_scalar_field(&glattice);
    out = alloc_suNg_scalar_field(&glattice);

    random_suNg_scalar_field_cpu(in);
    lprintf("SANITY CHECK", 0, "[In field CPU copy norm unequal zero: %0.2e]\n", sqnorm_suNg_scalar_field_cpu(in));

    // Convert twice
    to_gpu_format_suNg_scalar_field(tmp, in);
    fill_buffers_suNg_scalar_field(tmp);
    to_cpu_format_suNg_scalar_field(out, tmp);
    lprintf("SANITY CHECK", 0, "[In and outfield sqnorm unequal zero and equal to each other: in %0.2e out %0.2e]\n",
            sqnorm_suNg_scalar_field_cpu(in), sqnorm_suNg_scalar_field_cpu(out));

    // Assert field are equal over sqnorm
    sub_assign_suNg_scalar_field_cpu(out, in);
    double diff_norm = sqnorm_suNg_scalar_field_cpu(out);

    free_suNg_scalar_field(in);
    free_suNg_scalar_field(tmp);
    free_suNg_scalar_field(out);
    return check_diff_norm_zero(diff_norm);
}

int test_convert_back_forth_suNg_av_field() {
    lprintf("INFO", 0, " ======= TEST suNg_av_field ======= \n");

    // Setup fields
    suNg_av_field *in, *tmp, *out;
    in = alloc_suNg_av_field(&glattice);
    tmp = alloc_suNg_av_field(&glattice);
    out = alloc_suNg_av_field(&glattice);

    random_suNg_av_field_cpu(in);
    lprintf("SANITY CHECK", 0, "[In field CPU copy norm unequal zero: %0.2e]\n", sqnorm_suNg_av_field_cpu(in));

    // Convert twice
    to_gpu_format_suNg_av_field(tmp, in);
    fill_buffers_suNg_av_field(tmp);
    to_cpu_format_suNg_av_field(out, tmp);
    lprintf("SANITY CHECK", 0, "[In and outfield sqnorm unequal zero and equal to each other: in %0.2e out %0.2e]\n",
            sqnorm_suNg_av_field_cpu(in), sqnorm_suNg_av_field_cpu(out));

    // Assert field are equal over sqnorm
    sub_assign_suNg_av_field_cpu(out, in);
    double diff_norm = sqnorm_suNg_av_field_cpu(out);

    free_suNg_av_field(in);
    free_suNg_av_field(tmp);
    free_suNg_av_field(out);
    return check_diff_norm_zero(diff_norm);
}

int test_convert_back_forth_gtransf() {
    lprintf("INFO", 0, " ======= TEST GTRANSF ======= \n");

    // Setup fields
    suNg_field *in, *tmp, *out;
    in = alloc_gtransf(&glattice);
    tmp = alloc_gtransf(&glattice);
    out = alloc_gtransf(&glattice);

    random_gtransf_cpu(in);
    lprintf("SANITY CHECK", 0, "[In field CPU copy norm unequal zero: %0.2e]\n", sqnorm_gtransf_cpu(in));

    // Convert twice
    to_gpu_format_gtransf(tmp, in);
    fill_buffers_gtransf(tmp);
    to_cpu_format_gtransf(out, tmp);
    lprintf("SANITY CHECK", 0, "[In and outfield sqnorm unequal zero and equal to each other: in %0.2e out %0.2e]\n",
            sqnorm_gtransf_cpu(in), sqnorm_gtransf_cpu(out));

    // Assert fields are equal over sqnorm
    sub_assign_gtransf_cpu(out, in);
    double diff_norm = sqnorm_gtransf_cpu(out);

    free_gtransf(in);
    free_gtransf(tmp);
    free_gtransf(out);
    return check_diff_norm_zero(diff_norm);
}

int test_convert_back_forth_clover_term() {
    lprintf("INFO", 0, " ======= TEST CLOVER TERM ======= \n");

    // Setup fields
    suNfc_field *in, *tmp, *out;
    in = alloc_clover_term(&glattice);
    tmp = alloc_clover_term(&glattice);
    out = alloc_clover_term(&glattice);

    random_clover_term_cpu(in);
    lprintf("SANITY CHECK", 0, "[In field CPU copy norm unequal zero: %0.2e]\n", sqnorm_clover_term_cpu(in));

    // Convert twice
    to_gpu_format_clover_term(tmp, in);
    fill_buffers_clover_term(tmp);
    to_cpu_format_clover_term(out, tmp);
    lprintf("SANITY CHECK", 0, "[In and outfield sqnorm unequal zero and equal to each other: in %0.2e out %0.2e]\n",
            sqnorm_clover_term_cpu(in), sqnorm_clover_term_cpu(out));

    // Assert fields are equal over sqnorm
    sub_assign_clover_term_cpu(out, in);
    double diff_norm = sqnorm_clover_term_cpu(out);

    free_clover_term(in);
    free_clover_term(tmp);
    free_clover_term(out);
    return check_diff_norm_zero(diff_norm);
}

int test_convert_back_forth_clover_force() {
    lprintf("INFO", 0, " ======= TEST CLOVER FORCE ======= \n");

    // Setup fields
    suNf_field *in, *tmp, *out;
    in = alloc_clover_force(&glattice);
    tmp = alloc_clover_force(&glattice);
    out = alloc_clover_force(&glattice);

    random_clover_force_cpu(in);
    lprintf("SANITY CHECK", 0, "[In field CPU copy norm unequal zero: %0.2e]\n", sqnorm_clover_force_cpu(in));

    // Convert twice
    to_gpu_format_clover_force(tmp, in);
    fill_buffers_clover_force(tmp);
    to_cpu_format_clover_force(out, tmp);
    lprintf("SANITY CHECK", 0, "[In and outfield sqnorm unequal zero and equal to each other: in %0.2e out %0.2e]\n",
            sqnorm_clover_force_cpu(in), sqnorm_clover_force_cpu(out));

    // Assert fields are equal over sqnorm
    sub_assign_clover_force_cpu(out, in);
    double diff_norm = sqnorm_clover_force_cpu(out);

    free_clover_force(in);
    free_clover_force(tmp);
    free_clover_force(out);
    return check_diff_norm_zero(diff_norm);
}

int test_convert_back_forth_scalar_field() {
    lprintf("INFO", 0, " ======= TEST SFIELD ======= \n");

    // Setup fields
    scalar_field *in, *tmp, *out;
    in = alloc_scalar_field(1, &glattice);
    tmp = alloc_scalar_field(1, &glattice);
    out = alloc_scalar_field(1, &glattice);

    random_scalar_field_cpu(in);
    lprintf("SANITY CHECK", 0, "[In field CPU copy norm unequal zero: %0.2e]\n", sqnorm_scalar_field_cpu(in));

    // Convert twice
    to_gpu_format_scalar_field(tmp, in);
    fill_buffers_scalar_field(tmp);
    to_cpu_format_scalar_field(out, tmp);
    lprintf("SANITY CHECK", 0, "[In and outfield sqnorm unequal zero and equal to each other: in %0.2e out %0.2e]\n",
            sqnorm_scalar_field_cpu(in), sqnorm_scalar_field_cpu(out));

    // Assert fields are equal over sqnorm
    sub_assign_scalar_field_cpu(out, in);
    double diff_norm = sqnorm_scalar_field_cpu(out);

    free_scalar_field(in);
    free_scalar_field(tmp);
    free_scalar_field(out);
    return check_diff_norm_zero(diff_norm);
}

int test_convert_back_forth_ldl_field() {
    lprintf("INFO", 0, " ======= TEST CLOVER LDL ======= \n");

    // Setup fields
    ldl_field *in, *tmp, *out;
    in = alloc_ldl_field(&glattice);
    tmp = alloc_ldl_field(&glattice);
    out = alloc_ldl_field(&glattice);

    random_ldl_field_cpu(in);
    lprintf("SANITY CHECK", 0, "[In field CPU copy norm unequal zero: %0.2e]\n", sqnorm_ldl_field_cpu(in));

    // Convert twice
    to_gpu_format_ldl_field(tmp, in);
    fill_buffers_ldl_field(tmp);
    to_cpu_format_ldl_field(out, tmp);

    lprintf("SANITY CHECK", 0, "[In and outfield sqnorm unequal zero and equal to each other: in %0.2e out %0.2e]\n",
            sqnorm_ldl_field_cpu(in), sqnorm_ldl_field_cpu(out));

    // Assert fields are equal over sqnorm
    sub_assign_ldl_field_cpu(out, in);
    double diff_norm = sqnorm_ldl_field_cpu(out);

    // Free and return
    free_ldl_field(in);
    free_ldl_field(tmp);
    free_ldl_field(out);
    return check_diff_norm_zero(diff_norm);
}

int test_convert_back_forth_staple_field() {
    lprintf("INFO", 0, " ======= TEST STAPLE FIELD ======= \n");

    // Setup suNg_fields
    suNg_field *in, *tmp, *out;
    in = alloc_staple_field(&glattice);
    tmp = alloc_staple_field(&glattice);
    out = alloc_staple_field(&glattice);

    random_staple_field_cpu(in);
    lprintf("SANITY CHECK", 0, "[In field CPU copy norm unequal zero: %0.2e]\n", sqnorm_staple_field_cpu(in));

    // Convert twice
    to_gpu_format_staple_field(tmp, in);
    fill_buffers_staple_field(tmp);
    to_cpu_format_staple_field(out, tmp);
    lprintf("SANITY CHECK", 0, "[In and outfield sqnorms unequal zero and equal to each other: in %0.2e out %0.2e]\n",
            sqnorm_staple_field_cpu(in), sqnorm_staple_field_cpu(out));

    // Assert fields are equal over sqnorm
    sub_assign_staple_field_cpu(out, in);
    double diff_norm = sqnorm_staple_field_cpu(out);

    // Free and return
    free_staple_field(in);
    free_staple_field(tmp);
    free_staple_field(out);
    return check_diff_norm_zero(diff_norm);
}