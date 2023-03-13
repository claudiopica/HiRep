/*******************************************************************************
*
* Single-double conversion utils tests
*
*******************************************************************************/

#include "libhr.h"

int test_bijectivity_spinors();
int test_add_assign();
int test_gfield();
int test_gfield_f();

int main(int argc, char *argv[]) {
    // Setup
    int return_val = 0;
    setup_process(&argc, &argv);
    setup_gauge_fields();

    // Test block
    return_val += test_bijectivity_spinors();
    return_val += test_add_assign();
    return_val += test_gfield();

#ifdef DPHI_FLT
    return_val += test_gfield_f();
#endif

    // Finalize
    finalize_process();
    return return_val;
}

int test_bijectivity_spinors() {
    // Setup fields
    lprintf("TEST", 0, "Testing spinor assign\n");
    int return_val = 0;
    double sqnorm, sqnorm_flt;
    spinor_field *in, *out;
    spinor_field_flt *in_flt, *out_flt;

    in = alloc_spinor_field_f(1, &glattice);
    out = alloc_spinor_field_f(1, &glattice);
    in_flt = alloc_spinor_field_f_flt(1, &glattice);
    out_flt = alloc_spinor_field_f_flt(1, &glattice);

    gaussian_spinor_field(in);
    gaussian_spinor_field_flt(in_flt);

    /*
        Assign to float then back to double. Then check
        1. The result is sitewise identical to the input field
           up to double precision
        2. The single precision intermediate field has the same
           square norm as the input field to single precision
    */

    assign_sd2s(out_flt, in);
    assign_s2sd(out, out_flt);
    spinor_field_sub_assign_f(out, in);
    return_val += check_diff_norm(spinor_field_sqnorm_f(out), 1.e-11);

    sqnorm = spinor_field_sqnorm_f(in);
    sqnorm_flt = spinor_field_sqnorm_f_flt(out_flt);
    return_val += check_diff_norm(sqnorm_flt - sqnorm, 1.e-4);

    /*
        Assign to double then back to float. Check
        1. The result is sitewise identical to the input field
           to double precision
        2. The double precision intermediate field has the same
           square norm as the input field to single precision.
    */

    assign_s2sd(out, in_flt);
    assign_sd2s(out_flt, out);
    spinor_field_sub_assign_f_flt(out_flt, in_flt);
    return_val += check_diff_norm(spinor_field_sqnorm_f_flt(out_flt), 1.e-11);

    sqnorm_flt = spinor_field_sqnorm_f_flt(in_flt);
    sqnorm = spinor_field_sqnorm_f(out);
    return_val += check_diff_norm(sqnorm_flt - sqnorm, 1.e-4);

    // Free fields
    free_spinor_field_f(in);
    free_spinor_field_f(out);
    free_spinor_field_f_flt(in_flt);
    free_spinor_field_f_flt(out_flt);
    return return_val;
}

int test_add_assign() {
    // Setup fields
    lprintf("TEST", 0, "Testing spinor add assign\n");
    int return_val = 0;
    spinor_field *in, *out;
    spinor_field_flt *in_flt, *out_flt;

    in = alloc_spinor_field_f(1, &glattice);
    out = alloc_spinor_field_f(1, &glattice);
    in_flt = alloc_spinor_field_f_flt(1, &glattice);
    out_flt = alloc_spinor_field_f_flt(1, &glattice);

    /*
        First assign to single precision then add assign
        back to double precision. The result should be
        identical to a field multiplication by 2 to
        double precision.
    */

    gaussian_spinor_field(in);
    spinor_field_mul_f(out, 2.0, in);

    assign_sd2s(in_flt, in);
    add_assign_s2sd(in, in_flt);
    spinor_field_sub_assign_f(out, in);
    return_val += check_diff_norm(spinor_field_sqnorm_f(out), 1.e-11);

    /*
        Assign random single precision field to double
        precision, then add assign back to single
        precision. The result should be identical to
        a field multiplication by 2. 
    */

    gaussian_spinor_field_flt(in_flt);
    spinor_field_mul_f_flt(out_flt, 2.0, in_flt);

    assign_s2sd(in, in_flt);
    add_assign_sd2s(in_flt, in);
    spinor_field_sub_assign_f_flt(out_flt, in_flt);
    return_val += check_diff_norm_zero(spinor_field_sqnorm_f_flt(out_flt));

    // Free fields
    free_spinor_field_f(in);
    free_spinor_field_f(out);
    free_spinor_field_f_flt(in_flt);
    free_spinor_field_f_flt(out_flt);

    return return_val;
}

int test_gfield() {
    // Setup fields
    lprintf("TEST", 0, "Testing gfield assign\n");
    int return_val = 0;
    random_u(u_gauge);
    suNg_field *u_gauge_copy = alloc_gfield(&glattice);
    u_gauge_flt = alloc_gfield_flt(&glattice);

    /*
        Create copy of current gauge field, then assign 
        gauge field to single precision field. Assign
        single precision gauge field back to double precision.
        The result should be sitewise identical to the
        original gauge field to double precision.
    */

#ifdef WITH_GPU
    cudaMemcpy(u_gauge_copy->gpu_ptr, u_gauge->gpu_ptr, 4 * u_gauge->type->gsize_gauge * sizeof(suNg),
               cudaMemcpyDeviceToDevice);
#else
    copy_gfield_cpu(u_gauge_copy, u_gauge);
#endif
    assign_ud2u();
    assign_u2ud();

#ifdef WITH_GPU
    copy_from_gpu_gfield(u_gauge);
    copy_from_gpu_gfield(u_gauge_copy);
#endif

    sub_assign_gfield_cpu(u_gauge_copy, u_gauge);
    return_val += check_diff_norm(sqnorm_gfield_cpu(u_gauge_copy), 1.e-11);

    // Free fields
    free_gfield(u_gauge_copy);
    free_gfield_flt(u_gauge_flt);
    return return_val;
}

int test_gfield_f() {
    // Setup fields
    lprintf("TEST", 0, "Testing gfield_f assign\n");
    int return_val = 0;
    random_u_f(u_gauge_f);
    suNf_field *u_gauge_copy = alloc_gfield_f(&glattice);

    /*
        Create copy of current gauge field, then assign 
        gauge field to single precision field. Assign
        single precision gauge field back to double precision.
        The result should be sitewise identical to the
        original gauge field to double precision.
    */

#ifdef WITH_GPU
    cudaMemcpy(u_gauge_copy->gpu_ptr, u_gauge_f->gpu_ptr, 4 * u_gauge_f->type->gsize_gauge * sizeof(suNf),
               cudaMemcpyDeviceToDevice);
#else
    copy_gfield_f_cpu(u_gauge_copy, u_gauge_f);
#endif
    assign_ud2u_f();
    assign_u2ud_f();

#ifdef WITH_GPU
    copy_from_gpu_gfield_f(u_gauge_f);
    copy_from_gpu_gfield_f(u_gauge_copy);
#endif

    sub_assign_gfield_f_cpu(u_gauge_copy, u_gauge_f);
    return_val += check_diff_norm(sqnorm_gfield_f_cpu(u_gauge_copy), 1.e-11);

    // Free fields
    free_gfield_f(u_gauge_copy);
    return return_val;
}
