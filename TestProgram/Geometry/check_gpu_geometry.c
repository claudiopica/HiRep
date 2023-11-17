/*******************************************************************************
*
* Check that the GPU reading and writing functions defined in suN.h 
* are bijective.
*
* NOCOMPILE = !WITH_GPU
* 
*******************************************************************************/

#include "libhr.h"

// TODO: Test for n > 1

// Double precision
int test_write_read_gauge_field();
int test_write_read_gauge_field_f();
int test_write_read_suNg_scalar_field();
int test_write_read_suNg_av_field();
int test_write_read_gtransf();
int test_write_read_clover_term();
int test_write_read_clover_force();
int test_write_read_spinor_field();
int test_write_read_scalar_field();
int test_write_read_ldl_field();
int test_write_read_staple_field();

// Single precision
int test_write_read_gauge_field_flt();
int test_write_read_gauge_field_f_flt();
int test_write_read_spinor_field_flt();

int main(int argc, char *argv[]) {
    // Init
    int return_val = 0;

    // Setup process and communication
    logger_map("DEBUG", "debug");
    setup_process(&argc, &argv);
    test_setup();

    // Double Precision Tests
    return_val += test_write_read_gauge_field();
    return_val += test_write_read_gauge_field_f();
    return_val += test_write_read_suNg_scalar_field();
    return_val += test_write_read_suNg_av_field();
    return_val += test_write_read_gtransf();
    return_val += test_write_read_ldl_field();
    return_val += test_write_read_clover_term();
    return_val += test_write_read_clover_force();
    return_val += test_write_read_staple_field();

    return_val += test_write_read_spinor_field();
    return_val += test_write_read_scalar_field();

    //Single Precision Tests
    return_val += test_write_read_gauge_field_flt();
    return_val += test_write_read_gauge_field_f_flt();
    return_val += test_write_read_spinor_field_flt();

    // Finalize and return
    finalize_process();
    return return_val;
}

int test_write_read_gauge_field() {
    lprintf("INFO", 0, " ======= TEST GAUGE FIELD ======= ");
    int return_val = 0;
    suNg_field *in, *gpu_format, *out;

    in = alloc_suNg_field(&glattice);
    out = alloc_suNg_field(&glattice);
    gpu_format = alloc_suNg_field(&glattice);

    random_u(in);
    lprintf("SANITY CHECK", 0, "[In field norm unequal zero: %0.2e]\n", sqnorm_suNg_field_cpu(in));

    suNg *in_mat, *out_mat;
    _MASTER_FOR(in->type, ix) {
        for (int comp = 0; comp < 4; comp++) {
            in_mat = _4FIELD_AT(in, ix, comp);
            out_mat = _4FIELD_AT(out, ix, comp);
            write_suNg_field_gpu(in_mat, gpu_format, ix, comp);
            read_suNg_field_gpu(out_mat, gpu_format, ix, comp);
        }
    }

    lprintf("SANITY CHECK", 0, "[Sanity check in field norm unequal zero: %0.2e]\n", sqnorm_suNg_field_cpu(in));
    lprintf("SANITY CHECK", 0, "[Sanity check out field norm unequal zero: %0.2e]\n", sqnorm_suNg_field_cpu(out));
    sub_assign_suNg_field_cpu(out, in);
    double diff_norm = sqnorm_suNg_field_cpu(out);
    return_val += check_diff_norm_zero(diff_norm);

    free_suNg_field(in);
    free_suNg_field(out);
    free_suNg_field(gpu_format);
    return return_val;
}

int test_write_read_gauge_field_flt() {
    lprintf("INFO", 0, " ======= TEST GAUGE FIELD SINGLE PRECISION ======= ");
    int return_val = 0;
    suNg_field_flt *in, *gpu_format, *out;

    in = alloc_suNg_field_flt(&glattice);
    out = alloc_suNg_field_flt(&glattice);
    gpu_format = alloc_suNg_field_flt(&glattice);

    random_suNg_field_flt_cpu(in);
    lprintf("SANITY CHECK", 0, "[In field norm unequal zero: %0.2e]\n", sqnorm_suNg_field_flt_cpu(in));

    suNg_flt *in_mat, *out_mat;
    _MASTER_FOR(in->type, ix) {
        for (int comp = 0; comp < 4; ++comp) {
            in_mat = _4FIELD_AT(in, ix, comp);
            out_mat = _4FIELD_AT(out, ix, comp);
            write_suNg_field_flt_gpu(in_mat, gpu_format, ix, comp);
            read_suNg_field_flt_gpu(out_mat, gpu_format, ix, comp);
        }
    }

    lprintf("SANITY CHECK", 0, "[Sanity check in field norm unequal zero: %0.2e]\n", sqnorm_suNg_field_flt_cpu(in));
    lprintf("SANITY CHECK", 0, "[Sanity check out field norm unequal zero: %0.2e]\n", sqnorm_suNg_field_flt_cpu(out));
    sub_assign_suNg_field_flt_cpu(out, in);
    double diff_norm = sqnorm_suNg_field_flt_cpu(out);
    check_diff_norm_zero(diff_norm);

    free_suNg_field_flt(in);
    free_suNg_field_flt(out);
    free_suNg_field_flt(gpu_format);
    return return_val;
}

int test_write_read_gauge_field_f() {
    lprintf("INFO", 0, " ======= TEST GAUGE FIELD FUNDAMENTAL REP ======= ");
    int return_val = 0;
    suNf_field *in, *gpu_format, *out;

    in = alloc_suNf_field(&glattice);
    out = alloc_suNf_field(&glattice);
    gpu_format = alloc_suNf_field(&glattice);

    random_u_f(in);
    lprintf("SANITY CHECK", 0, "[In field norm unequal zero: %0.2e]\n", sqnorm_suNf_field_cpu(in));

    suNf *in_mat, *out_mat;
    _MASTER_FOR(in->type, ix) {
        for (int comp = 0; comp < 4; ++comp) {
            in_mat = _4FIELD_AT(in, ix, comp);
            out_mat = _4FIELD_AT(out, ix, comp);
            write_suNf_field_gpu(in_mat, gpu_format, ix, comp);
            read_suNf_field_gpu(out_mat, gpu_format, ix, comp);
        }
    }

    lprintf("SANITY CHECK", 0, "[Sanity check in field norm unequal zero: %0.2e]\n", sqnorm_suNf_field_cpu(in));
    lprintf("SANITY CHECK", 0, "[Sanity check out field norm unequal zero: %0.2e]\n", sqnorm_suNf_field_cpu(out));
    sub_assign_suNf_field_cpu(out, in);
    double diff_norm = sqnorm_suNf_field_cpu(out);

    check_diff_norm_zero(diff_norm);

    free_suNf_field(in);
    free_suNf_field(out);
    free_suNf_field(gpu_format);
    return return_val;
}

int test_write_read_gauge_field_f_flt() {
    lprintf("INFO", 0, " ======= TEST GAUGE FIELD FUNDAMENTAL REP ======= ");
    int return_val = 0;
    suNf_field_flt *in, *gpu_format, *out;

    in = alloc_suNf_field_flt(&glattice);
    out = alloc_suNf_field_flt(&glattice);
    gpu_format = alloc_suNf_field_flt(&glattice);

    random_suNf_field_flt_cpu(in);
    lprintf("SANITY CHECK", 0, "[In field norm unequal zero: %0.2e]\n", sqnorm_suNf_field_flt_cpu(in));

    suNf_flt *in_mat, *out_mat;
    _MASTER_FOR(in->type, ix) {
        for (int comp = 0; comp < 4; ++comp) {
            in_mat = _4FIELD_AT(in, ix, comp);
            out_mat = _4FIELD_AT(out, ix, comp);
            write_suNf_field_flt_gpu(in_mat, gpu_format, ix, comp);
            read_suNf_field_flt_gpu(out_mat, gpu_format, ix, comp);
        }
    }

    lprintf("SANITY CHECK", 0, "[Sanity check in field norm unequal zero: %0.2e]\n", sqnorm_suNf_field_flt_cpu(in));
    lprintf("SANITY CHECK", 0, "[Sanity check out field norm unequal zero: %0.2e]\n", sqnorm_suNf_field_flt_cpu(out));
    sub_assign_suNf_field_flt_cpu(out, in);
    double diff_norm = sqnorm_suNf_field_flt_cpu(out);

    check_diff_norm_zero(diff_norm);

    free_suNf_field_flt(in);
    free_suNf_field_flt(out);
    free_suNf_field_flt(gpu_format);
    return return_val;
}

int test_write_read_spinor_field() {
    lprintf("INFO", 0, " ======= TEST SPINOR FIELD ======= ");
    int return_val = 0;
    spinor_field *in, *gpu_format, *out;

    in = alloc_spinor_field(1, &glattice);
    gpu_format = alloc_spinor_field(1, &glattice);
    out = alloc_spinor_field(1, &glattice);

    gaussian_spinor_field(in);

    lprintf("SANITY CHECK", 0, "[Sanity check in field norm unequal zero: %0.2e]\n", spinor_field_sqnorm_f_cpu(in));

    suNf_spinor *in_spinor, *out_spinor;
    _MASTER_FOR(in->type, ix) {
        in_spinor = _FIELD_AT(in, ix);
        out_spinor = _FIELD_AT(in, ix);
        write_spinor_field_gpu(in_spinor, gpu_format, ix, 0);
        read_spinor_field_gpu(out_spinor, gpu_format, ix, 0);
    }

    lprintf("SANITY CHECK", 0, "[Sanity check out field norm unequal zero: %0.2e]\n", spinor_field_sqnorm_f_cpu(out));
    spinor_field_sub_assign_f_cpu(out, in);
    double diff_norm = spinor_field_sqnorm_f_cpu(out);
    check_diff_norm_zero(diff_norm);

    free_spinor_field(in);
    free_spinor_field(gpu_format);
    free_spinor_field(out);
    return return_val;
}

int test_write_read_spinor_field_flt() {
    lprintf("INFO", 0, " ======= TEST SPINOR FIELD SINGLE PRECISION ======= \n");
    int return_val = 0;
    spinor_field_flt *in, *gpu_format, *out;

    in = alloc_spinor_field_flt(1, &glattice);
    gpu_format = alloc_spinor_field_flt(1, &glattice);
    out = alloc_spinor_field_flt(1, &glattice);

    gaussian_spinor_field_flt(in);

    lprintf("SANITY CHECK", 0, "[Sanity check in field norm unequal zero: %0.2e]\n", spinor_field_sqnorm_f_flt_cpu(in));

    suNf_spinor_flt *in_spinor, *out_spinor;
    _MASTER_FOR(in->type, ix) {
        in_spinor = _FIELD_AT(in, ix);
        out_spinor = _FIELD_AT(out, ix);
        write_spinor_field_flt_gpu(in_spinor, gpu_format, ix, 0);
        read_spinor_field_flt_gpu(out_spinor, gpu_format, ix, 0);
    }

    lprintf("SANITY CHECK", 0, "[Sanity check out field norm unequal zero: %0.2e]\n", spinor_field_sqnorm_f_flt_cpu(out));
    spinor_field_sub_assign_f_flt_cpu(out, in);
    double diff_norm = spinor_field_sqnorm_f_flt_cpu(out);
    check_diff_norm_zero(diff_norm);

    free_spinor_field_flt(in);
    free_spinor_field_flt(gpu_format);
    free_spinor_field_flt(out);
    return return_val;
}

int test_write_read_scalar_field() {
    lprintf("INFO", 0, " ======= TEST SCALAR FIELD ======= ");
    int return_val = 0;
    scalar_field *in, *gpu_format, *out;

    in = alloc_scalar_field(1, &glattice);
    gpu_format = alloc_scalar_field(1, &glattice);
    out = alloc_scalar_field(1, &glattice);

    random_scalar_field_cpu(in);
    lprintf("SANITY CHECK", 0, "[Sanity check in field norm unequal zero: %0.2e]\n", sqnorm_scalar_field_cpu(in));

    double *in_site, *out_site;
    _MASTER_FOR(in->type, ix) {
        in_site = _FIELD_AT(in, ix);
        out_site = _FIELD_AT(out, ix);
        write_scalar_field_gpu(in_site, gpu_format, ix, 0);
        read_scalar_field_gpu(out_site, gpu_format, ix, 0);
    }

    lprintf("SANITY CHECK", 0, "[Sanity check in field norm unequal zero: %0.2e]\n", sqnorm_scalar_field_cpu(in));
    lprintf("SANITY CHECK", 0, "[Sanity check out field norm unequal zero: %0.2e]\n", sqnorm_scalar_field_cpu(out));
    sub_assign_scalar_field_cpu(out, in);
    double diff_norm = sqnorm_scalar_field_cpu(out);

    check_diff_norm_zero(diff_norm);

    free_scalar_field(in);
    free_scalar_field(gpu_format);
    free_scalar_field(out);
    return return_val;
}

int test_write_read_suNg_scalar_field() {
    lprintf("INFO", 0, " ======= TEST SU(NG) SCALAR FIELD ======= ");
    int return_val = 0;
    suNg_scalar_field *in, *gpu_format, *out;

    in = alloc_suNg_scalar_field(&glattice);
    out = alloc_suNg_scalar_field(&glattice);
    gpu_format = alloc_suNg_scalar_field(&glattice);

    random_suNg_scalar_field_cpu(in);
    lprintf("SANITY CHECK", 0, "[In field norm unequal zero: %0.2e]\n", sqnorm_suNg_scalar_field_cpu(in));

    suNg_vector *in_vec, *out_vec;
    _MASTER_FOR(in->type, ix) {
        in_vec = _FIELD_AT(in, ix);
        out_vec = _FIELD_AT(out, ix);
        write_suNg_scalar_field_gpu(in_vec, gpu_format, ix, 0);
        read_suNg_scalar_field_gpu(out_vec, gpu_format, ix, 0);
    }

    lprintf("SANITY CHECK", 0, "[Sanity check in field norm unequal zero: %0.2e]\n", sqnorm_suNg_scalar_field_cpu(in));
    lprintf("SANITY CHECK", 0, "[Sanity check out field norm unequal zero: %0.2e]\n", sqnorm_suNg_scalar_field_cpu(out));
    sub_assign_suNg_scalar_field_cpu(out, in);
    double diff_norm = sqnorm_suNg_scalar_field_cpu(out);

    check_diff_norm_zero(diff_norm);

    free_suNg_scalar_field(in);
    free_suNg_scalar_field(gpu_format);
    free_suNg_scalar_field(out);
    return return_val;
}

int test_write_read_suNg_av_field() {
    lprintf("INFO", 0, " ======= TEST suNg_av_field ======= ");
    int return_val = 0;
    suNg_av_field *in, *gpu_format, *out;

    in = alloc_suNg_av_field(&glattice);
    out = alloc_suNg_av_field(&glattice);
    gpu_format = alloc_suNg_av_field(&glattice);

    random_suNg_av_field_cpu(in);
    lprintf("SANITY CHECK", 0, "[In field norm unequal zero: %0.2e]\n", sqnorm_suNg_av_field_cpu(in));

    suNg_algebra_vector *in_mat, *out_mat;
    _MASTER_FOR(in->type, ix) {
        for (int comp = 0; comp < 4; ++comp) {
            in_mat = _4FIELD_AT(in, ix, comp);
            out_mat = _4FIELD_AT(out, ix, comp);
            write_suNg_av_field_gpu(in_mat, gpu_format, ix, comp);
            read_suNg_av_field_gpu(out_mat, gpu_format, ix, comp);
        }
    }

    lprintf("SANITY CHECK", 0, "[Sanity check in field norm unequal zero: %0.2e]\n", sqnorm_suNg_av_field_cpu(in));
    lprintf("SANITY CHECK", 0, "[Sanity check out field norm unequal zero: %0.2e]\n", sqnorm_suNg_av_field_cpu(out));
    sub_assign_suNg_av_field_cpu(out, in);
    double diff_norm = sqnorm_suNg_av_field_cpu(out);

    check_diff_norm_zero(diff_norm);

    free_suNg_av_field(in);
    free_suNg_av_field(out);
    free_suNg_av_field(gpu_format);
    return return_val;
}

int test_write_read_gtransf() {
    lprintf("INFO", 0, " ======= TEST GTRANSF ======= ");
    int return_val = 0;
    gtransf *in, *gpu_format, *out;

    in = alloc_gtransf(&glattice);
    out = alloc_gtransf(&glattice);
    gpu_format = alloc_gtransf(&glattice);

    random_gtransf_cpu(in);
    lprintf("SANITY CHECK", 0, "[In field norm unequal zero: %0.2e]\n", sqnorm_gtransf_cpu(in));

    suNg *in_mat, *out_mat;
    _MASTER_FOR(in->type, ix) {
        in_mat = _FIELD_AT(in, ix);
        out_mat = _FIELD_AT(out, ix);
        write_gtransf_gpu(in_mat, gpu_format, ix, 0);
        read_gtransf_gpu(out_mat, gpu_format, ix, 0);
    }

    lprintf("SANITY CHECK", 0, "[Sanity check in field norm unequal zero: %0.2e]\n", sqnorm_gtransf_cpu(in));
    lprintf("SANITY CHECK", 0, "[Sanity check out field norm unequal zero: %0.2e]\n", sqnorm_gtransf_cpu(out));
    sub_assign_gtransf_cpu(out, in);
    double diff_norm = sqnorm_gtransf_cpu(out);

    check_diff_norm_zero(diff_norm);

    free_gtransf(in);
    free_gtransf(out);
    free_gtransf(gpu_format);
    return return_val;
}

int test_write_read_clover_term() {
    lprintf("INFO", 0, " ======= TEST CLOVER TERM ======= ");
    int return_val = 0;
    clover_term *in, *gpu_format, *out;

    in = alloc_clover_term(&glattice);
    out = alloc_clover_term(&glattice);
    gpu_format = alloc_clover_term(&glattice);

    random_clover_term_cpu(in);
    lprintf("SANITY CHECK", 0, "[In field norm unequal zero: %0.2e]\n", sqnorm_clover_term_cpu(in));

    suNfc *in_mat, *out_mat;
    _MASTER_FOR(in->type, ix) {
        for (int comp = 0; comp < 4; ++comp) {
            in_mat = _4FIELD_AT(in, ix, comp);
            out_mat = _4FIELD_AT(out, ix, comp);
            write_clover_term_gpu(in_mat, gpu_format, ix, comp);
            read_clover_term_gpu(out_mat, gpu_format, ix, comp);
        }
    }

    lprintf("SANITY CHECK", 0, "[Sanity check in field norm unequal zero: %0.2e]\n", sqnorm_clover_term_cpu(in));
    lprintf("SANITY CHECK", 0, "[Sanity check out field norm unequal zero: %0.2e]\n", sqnorm_clover_term_cpu(out));
    sub_assign_clover_term_cpu(out, in);
    double diff_norm = sqnorm_clover_term_cpu(out);

    check_diff_norm_zero(diff_norm);

    free_clover_term(in);
    free_clover_term(out);
    free_clover_term(gpu_format);
    return return_val;
}

int test_write_read_clover_force() {
    lprintf("INFO", 0, " ======= TEST CLOVER FORCE ======= ");
    int return_val = 0;
    clover_force *in, *gpu_format, *out;

    in = alloc_clover_force(&glattice);
    out = alloc_clover_force(&glattice);
    gpu_format = alloc_clover_force(&glattice);

    random_clover_force_cpu(in);
    lprintf("SANITY CHECK", 0, "[In field norm unequal zero: %0.2e]\n", sqnorm_clover_force_cpu(in));

    suNf *in_mat, *out_mat;
    _MASTER_FOR(in->type, ix) {
        for (int comp = 0; comp < 6; ++comp) {
            in_mat = _6FIELD_AT(in, ix, comp);
            out_mat = _6FIELD_AT(out, ix, comp);
            write_clover_force_gpu(in_mat, gpu_format, ix, comp);
            read_clover_force_gpu(out_mat, gpu_format, ix, comp);
        }
    }

    lprintf("SANITY CHECK", 0, "[Sanity check in field norm unequal zero: %0.2e]\n", sqnorm_clover_force_cpu(in));
    lprintf("SANITY CHECK", 0, "[Sanity check out field norm unequal zero: %0.2e]\n", sqnorm_clover_force_cpu(out));
    sub_assign_clover_force_cpu(out, in);
    double diff_norm = sqnorm_clover_force_cpu(out);

    check_diff_norm_zero(diff_norm);

    free_clover_force(in);
    free_clover_force(out);
    free_clover_force(gpu_format);
    return return_val;
}

int test_write_read_ldl_field() {
    lprintf("INFO", 0, " ======= TEST LDL FIELD ======= ");
    int return_val = 0;
    ldl_field *in, *gpu_format, *out;

    in = alloc_ldl_field(&glattice);
    out = alloc_ldl_field(&glattice);
    gpu_format = alloc_ldl_field(&glattice);

    random_ldl_field_cpu(in);
    lprintf("SANITY CHECK", 0, "[In field norm unequal zero: %0.2e]\n", sqnorm_ldl_field_cpu(in)); // sqnorm Not implemented

    ldl_t *site_in, *site_out;
    _MASTER_FOR(in->type, ix) {
        site_in = _FIELD_AT(in, ix);
        site_out = _FIELD_AT(out, ix);
        write_ldl_field_gpu(site_in, gpu_format, ix, 0);
        read_ldl_field_gpu(site_out, gpu_format, ix, 0);
    }

    lprintf("SANITY CHECK", 0, "[Sanity check in field norm unequal zero: %0.2e]\n", sqnorm_ldl_field_cpu(in));
    lprintf("SANITY CHECK", 0, "[Sanity check out field norm unequal zero: %0.2e]\n", sqnorm_ldl_field_cpu(out));

    sub_assign_ldl_field_cpu(out, in);
    double diff_norm = sqnorm_ldl_field_cpu(out);
    check_diff_norm_zero(diff_norm);

    free_ldl_field(in);
    free_ldl_field(out);
    free_ldl_field(gpu_format);
    return return_val;
}

int test_write_read_staple_field() {
    lprintf("INFO", 0, " ======= TEST STAPLE FIELD ======= ");
    int return_val = 0;
    staple_field *in, *gpu_format, *out;

    in = alloc_staple_field(&glattice);
    out = alloc_staple_field(&glattice);
    gpu_format = alloc_staple_field(&glattice);

    random_staple_field_cpu(in);
    lprintf("SANITY CHECK", 0, "[In field norm unequal zero: %0.2e]\n", sqnorm_staple_field_cpu(in));

    suNg *in_mat, *out_mat;
    _MASTER_FOR(in->type, ix) {
        for (int comp = 0; comp < 3; comp++) {
            in_mat = _3FIELD_AT(in, ix, comp);
            out_mat = _3FIELD_AT(out, ix, comp);
            write_staple_field_gpu(in_mat, gpu_format, ix, comp);
            read_staple_field_gpu(out_mat, gpu_format, ix, comp);
        }
    }

    lprintf("SANITY CHECK", 0, "[Sanity check in field norm unequal zero: %0.2e]\n", sqnorm_staple_field_cpu(in));
    lprintf("SANITY CHECK", 0, "[Sanity check out field norm unequal zero: %0.2e]\n", sqnorm_staple_field_cpu(out));
    sub_assign_staple_field_cpu(out, in);
    double diff_norm = sqnorm_staple_field_cpu(out);
    return_val += check_diff_norm_zero(diff_norm);

    free_staple_field(in);
    free_staple_field(out);
    free_staple_field(gpu_format);
    return return_val;
}
