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
int test_convert_back_forth_gfield();
int test_convert_back_forth_gfield_f();
int test_convert_back_forth_suNg_scalar_field();
int test_convert_back_forth_avfield();
int test_convert_back_forth_gtransf();
int test_convert_back_forth_clover_ldl();
int test_convert_back_forth_clover_term();
int test_convert_back_forth_clover_force();
int test_convert_back_forth_spinor_field();
int test_convert_back_forth_sfield();

/* Single precision tests */
int test_convert_back_forth_gfield_flt();
int test_convert_back_forth_gfield_f_flt();
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
    return_val += test_convert_back_forth_gfield();
    return_val += test_convert_back_forth_gfield_f();
    return_val += test_convert_back_forth_suNg_scalar_field();
    return_val += test_convert_back_forth_avfield();
    return_val += test_convert_back_forth_gtransf();
    return_val += test_convert_back_forth_clover_ldl();
    return_val += test_convert_back_forth_clover_term();
    return_val += test_convert_back_forth_clover_force();
    return_val += test_convert_back_forth_spinor_field();
    return_val += test_convert_back_forth_sfield();

    /* Single precision */
    return_val += test_convert_back_forth_gfield_flt();
    return_val += test_convert_back_forth_gfield_f_flt();
    return_val += test_convert_back_forth_spinor_field_flt();

    // Finalize and return
    finalize_process();
    return return_val;
}

int test_convert_back_forth_spinor_field() {
    lprintf("INFO", 0, " ======= TEST SPINOR FIELD ======= \n");

    // Setup spinor fields
    spinor_field *in, *tmp, *out;
    in = alloc_spinor_field_f(1, &glattice);
    tmp = alloc_spinor_field_f(1, &glattice);
    out = alloc_spinor_field_f(1, &glattice);
    gaussian_spinor_field(in);

    lprintf("SANITY CHECK", 0, "[In field CPU copy norm unequal zero: %0.2e]\n", spinor_field_sqnorm_f_cpu(in));

    // Convert twice
    to_gpu_format_spinor_field_f(tmp, in);
    fill_buffers_spinor_field_f(tmp);
    to_cpu_format_spinor_field_f(out, tmp);

    lprintf("SANITY CHECK", 0, "[In and outfield sqnorms unequal zero and equal to each other: in %0.2e out %0.2e]\n",
            spinor_field_sqnorm_f_cpu(in), spinor_field_sqnorm_f_cpu(out));

    // Assert fields are equal over sqnorm
    spinor_field_sub_assign_f_cpu(out, in);
    double diff_norm = spinor_field_sqnorm_f_cpu(out);

    // Free and return
    free_spinor_field_f(in);
    free_spinor_field_f(tmp);
    free_spinor_field_f(out);
    return check_diff_norm_zero(diff_norm);
}

int test_convert_back_forth_spinor_field_flt() {
    lprintf("INFO", 0, " ======= TEST SPINOR FIELD SINGLE PRECISION ======= \n");

    // Setup spinor fields
    spinor_field_flt *in, *tmp, *out;
    in = alloc_spinor_field_f_flt(1, &glattice);
    tmp = alloc_spinor_field_f_flt(1, &glattice);
    out = alloc_spinor_field_f_flt(1, &glattice);
    gaussian_spinor_field_flt(in);

    double sqnorm = spinor_field_sqnorm_f_flt_cpu(in);
    lprintf("SANITY CHECK", 0, "[In field CPU copy norm unequal zero: %0.2e]\n", sqnorm);

    // Convert twice
    to_gpu_format_spinor_field_f_flt(tmp, in);
    fill_buffers_spinor_field_f_flt(tmp);
    to_cpu_format_spinor_field_f_flt(out, tmp);
    lprintf("SANITY CHECK", 0, "[In and outfield sqnorms unequal zero and equal to each other: in: %0.2e out: %0.2e]\n",
            spinor_field_sqnorm_f_flt_cpu(in), spinor_field_sqnorm_f_flt_cpu(out));

    // Assert fields are equal over sqnorm
    spinor_field_sub_assign_f_flt_cpu(out, in);
    double diff_norm = spinor_field_sqnorm_f_flt_cpu(out);

    // Free and return
    free_spinor_field_f_flt(in);
    free_spinor_field_f_flt(tmp);
    free_spinor_field_f_flt(out);
    return check_diff_norm_zero(diff_norm);
}

int test_convert_back_forth_gfield_f() {
    lprintf("INFO", 0, " ======= TEST GAUGE FIELD REPRESENTED ======= \n");

    // Setup gfields
    suNf_field *in, *tmp, *out;
    in = alloc_gfield_f(&glattice);

    tmp = alloc_gfield_f(&glattice);
    out = alloc_gfield_f(&glattice);

    random_u_f(in);
    lprintf("SANITY CHECK", 0, "[In field CPU copy norm unequal zero: %0.2e]\n", sqnorm_gfield_f_cpu(in));

    // Convert twice
    to_gpu_format_gfield_f(tmp, in);
    fill_buffers_gfield_f(tmp);
    to_cpu_format_gfield_f(out, tmp);

    lprintf("SANITY CHECK", 0, "[In and outfield sqnorms unequal zero and equal to each other: in %0.2e out %0.2e]\n",
            sqnorm_gfield_f_cpu(in), sqnorm_gfield_f_cpu(out));

    // Assert fields are equal over sqnorm
    sub_assign_gfield_f_cpu(out, in);
    double diff_norm = sqnorm_gfield_f_cpu(out);
    return check_diff_norm_zero(diff_norm);

    free_gfield_f(in);
    free_gfield_f(tmp);
    free_gfield_f(out);
    return 0;
}

int test_convert_back_forth_gfield_f_flt() {
    lprintf("INFO", 0, " ======= TEST GAUGE FIELD REPRESENTED SINGLE PRECISION ======= \n");

    // Setup gfields
    suNf_field_flt *in, *tmp, *out;
    in = alloc_gfield_f_flt(&glattice);
    tmp = alloc_gfield_f_flt(&glattice);
    out = alloc_gfield_f_flt(&glattice);

    random_gfield_f_flt_cpu(in);
    lprintf("SANITY CHECK", 0, "[In field CPU copy norm unequal zero: %0.2e]\n", sqnorm_gfield_f_flt_cpu(in));

    // Convert twice
    to_gpu_format_gfield_f_flt(tmp, in);
    fill_buffers_gfield_f_flt(tmp);
    to_cpu_format_gfield_f_flt(out, tmp);
    lprintf("SANITY CHECK", 0, "[In and outfield sqnorms unequal zero and equal to each other: in %0.2e out %0.2e]\n",
            sqnorm_gfield_f_flt_cpu(in), sqnorm_gfield_f_flt_cpu(out));

    // Assert fields are equal over sqnorm
    sub_assign_gfield_f_flt_cpu(out, in);
    double diff_norm = sqnorm_gfield_f_flt_cpu(out);

    // Free and return
    free_gfield_f_flt(in);
    free_gfield_f_flt(tmp);
    free_gfield_f_flt(out);
    return check_diff_norm_zero(diff_norm);
}

int test_convert_back_forth_gfield() {
    lprintf("INFO", 0, " ======= TEST GAUGE FIELD ======= \n");

    // Setup gfields
    suNg_field *in, *tmp, *out;
    in = alloc_gfield(&glattice);
    tmp = alloc_gfield(&glattice);
    out = alloc_gfield(&glattice);

    random_u(in);
    lprintf("SANITY CHECK", 0, "[In field CPU copy norm unequal zero: %0.2e]\n", sqnorm_gfield_cpu(in));

    // Convert twice
    to_gpu_format_gfield(tmp, in);
    fill_buffers_gfield(tmp);
    to_cpu_format_gfield(out, tmp);
    lprintf("SANITY CHECK", 0, "[In and outfield sqnorms unequal zero and equal to each other: in %0.2e out %0.2e]\n",
            sqnorm_gfield_cpu(in), sqnorm_gfield_cpu(out));

    // Assert fields are equal over sqnorm
    sub_assign_gfield_cpu(out, in);
    double diff_norm = sqnorm_gfield_cpu(out);

    // Free and return
    free_gfield(in);
    free_gfield(tmp);
    free_gfield(out);
    return check_diff_norm_zero(diff_norm);
}

int test_convert_back_forth_gfield_flt() {
    lprintf("INFO", 0, " ======= TEST GAUGE FIELD SINGLE PRECISION ======= \n");

    // Setup gfields
    suNg_field_flt *in, *tmp, *out;
    in = alloc_gfield_flt(&glattice);
    tmp = alloc_gfield_flt(&glattice);
    out = alloc_gfield_flt(&glattice);

    random_gfield_flt_cpu(in);
    lprintf("SANITY CHECK", 0, "[In field CPU copy norm unequal zero: %0.2e]\n", sqnorm_gfield_flt_cpu(in));

    // Convert twice
    to_gpu_format_gfield_flt(tmp, in);
    fill_buffers_gfield_flt(tmp);
    to_cpu_format_gfield_flt(out, tmp);
    lprintf("SANITY CHECK", 0, "[In and outfield sqnorms unequal zero and equal to each other: in %0.2e out %0.2e]\n",
            sqnorm_gfield_flt_cpu(in), sqnorm_gfield_flt_cpu(out));

    // Assert fields are equal over sqnorm
    sub_assign_gfield_flt_cpu(out, in);
    double diff_norm = sqnorm_gfield_flt_cpu(out);

    free_gfield_flt(in);
    free_gfield_flt(tmp);
    free_gfield_flt(out);
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

int test_convert_back_forth_avfield() {
    lprintf("INFO", 0, " ======= TEST AVFIELD ======= \n");

    // Setup fields
    suNg_av_field *in, *tmp, *out;
    in = alloc_avfield(&glattice);
    tmp = alloc_avfield(&glattice);
    out = alloc_avfield(&glattice);

    random_avfield_cpu(in);
    lprintf("SANITY CHECK", 0, "[In field CPU copy norm unequal zero: %0.2e]\n", sqnorm_avfield_cpu(in));

    // Convert twice
    to_gpu_format_avfield(tmp, in);
    fill_buffers_avfield(tmp);
    to_cpu_format_avfield(out, tmp);
    lprintf("SANITY CHECK", 0, "[In and outfield sqnorm unequal zero and equal to each other: in %0.2e out %0.2e]\n",
            sqnorm_avfield_cpu(in), sqnorm_avfield_cpu(out));

    // Assert field are equal over sqnorm
    sub_assign_avfield_cpu(out, in);
    double diff_norm = sqnorm_avfield_cpu(out);

    free_avfield(in);
    free_avfield(tmp);
    free_avfield(out);
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

int test_convert_back_forth_sfield() {
    lprintf("INFO", 0, " ======= TEST SFIELD ======= \n");

    // Setup fields
    scalar_field *in, *tmp, *out;
    in = alloc_sfield(1, &glattice);
    tmp = alloc_sfield(1, &glattice);
    out = alloc_sfield(1, &glattice);

    random_sfield_cpu(in);
    lprintf("SANITY CHECK", 0, "[In field CPU copy norm unequal zero: %0.2e]\n", sqnorm_sfield_cpu(in));

    // Convert twice
    to_gpu_format_sfield(tmp, in);
    fill_buffers_sfield(tmp);
    to_cpu_format_sfield(out, tmp);
    lprintf("SANITY CHECK", 0, "[In and outfield sqnorm unequal zero and equal to each other: in %0.2e out %0.2e]\n",
            sqnorm_sfield_cpu(in), sqnorm_sfield_cpu(out));

    // Assert fields are equal over sqnorm
    sub_assign_sfield_cpu(out, in);
    double diff_norm = sqnorm_sfield_cpu(out);

    free_sfield(in);
    free_sfield(tmp);
    free_sfield(out);
    return check_diff_norm_zero(diff_norm);
}

int test_convert_back_forth_clover_ldl() {
    lprintf("INFO", 0, " ======= TEST CLOVER LDL ======= \n");

    // Setup fields
    ldl_field *in, *tmp, *out;
    in = alloc_clover_ldl(&glattice);
    tmp = alloc_clover_ldl(&glattice);
    out = alloc_clover_ldl(&glattice);

    random_clover_ldl_cpu(in);
    lprintf("SANITY CHECK", 0, "[In field CPU copy norm unequal zero: %0.2e]\n", sqnorm_clover_ldl_cpu(in));

    // Convert twice
    to_gpu_format_clover_ldl(tmp, in);
    fill_buffers_clover_ldl(tmp);
    to_cpu_format_clover_ldl(out, tmp);

    lprintf("SANITY CHECK", 0, "[In and outfield sqnorm unequal zero and equal to each other: in %0.2e out %0.2e]\n",
            sqnorm_clover_ldl_cpu(in), sqnorm_clover_ldl_cpu(out));

    // Assert fields are equal over sqnorm
    sub_assign_clover_ldl_cpu(out, in);
    double diff_norm = sqnorm_clover_ldl_cpu(out);

    // Free and return
    free_clover_ldl(in);
    free_clover_ldl(tmp);
    free_clover_ldl(out);
    return check_diff_norm_zero(diff_norm);
}