/*******************************************************************************
*
* Single-double conversion utils tests
*
*******************************************************************************/

#include "libhr.h"

int test_bijectivity_spinors();
int test_add_assign();
int test_suNg_field();
int test_suNf_field();

int main(int argc, char *argv[]) {
    // Setup
    int return_val = 0;
    setup_process(&argc, &argv);
    setup_gauge_fields();
    represent_gauge_field();
    random_u(u_gauge);
    random_u_f(u_gauge_f);

    // Test block
    return_val += test_bijectivity_spinors();
    return_val += test_add_assign();
    return_val += test_suNg_field();

#ifdef DPHI_FLT
    return_val += test_suNf_field();
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

    in = alloc_spinor_field(1, &glattice);
    out = alloc_spinor_field(1, &glattice);
    in_flt = alloc_spinor_field_flt(1, &glattice);
    out_flt = alloc_spinor_field_flt(1, &glattice);

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
    sub_assign_spinor_field(out, in);
    return_val += check_diff_norm(sqnorm_spinor_field(out), 1.e-10);

    sqnorm = sqnorm_spinor_field(in);
    sqnorm_flt = sqnorm_spinor_field_flt(out_flt);
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
    sub_assign_spinor_field_flt(out_flt, in_flt);
    return_val += check_diff_norm(sqnorm_spinor_field_flt(out_flt), 1.e-10);

    sqnorm_flt = sqnorm_spinor_field_flt(in_flt);
    sqnorm = sqnorm_spinor_field(out);
    return_val += check_diff_norm(sqnorm_flt - sqnorm, 1.e-4);

    // Free fields
    free_spinor_field(in);
    free_spinor_field(out);
    free_spinor_field_flt(in_flt);
    free_spinor_field_flt(out_flt);
    return return_val;
}

int test_add_assign() {
    // Setup fields
    lprintf("TEST", 0, "Testing spinor add assign\n");
    int return_val = 0;
    spinor_field *in, *out;
    spinor_field_flt *in_flt, *out_flt;

    in = alloc_spinor_field(1, &glattice);
    out = alloc_spinor_field(1, &glattice);
    in_flt = alloc_spinor_field_flt(1, &glattice);
    out_flt = alloc_spinor_field_flt(1, &glattice);

    /*
        First assign to single precision then add assign
        back to double precision. The result should be
        identical to a field multiplication by 2 to
        double precision.
    */

    gaussian_spinor_field(in);
    mul_spinor_field(out, 2.0, in);

    assign_sd2s(in_flt, in);
    add_assign_s2sd(in, in_flt);
    sub_assign_spinor_field(out, in);
    return_val += check_diff_norm(sqnorm_spinor_field(out), 1.e-10);

    /*
        Assign random single precision field to double
        precision, then add assign back to single
        precision. The result should be identical to
        a field multiplication by 2. 
    */

    gaussian_spinor_field_flt(in_flt);
    mul_spinor_field_flt(out_flt, 2.0, in_flt);

    assign_s2sd(in, in_flt);
    add_assign_sd2s(in_flt, in);
    sub_assign_spinor_field_flt(out_flt, in_flt);
    return_val += check_diff_norm_zero(sqnorm_spinor_field_flt(out_flt));

    // Free fields
    free_spinor_field(in);
    free_spinor_field(out);
    free_spinor_field_flt(in_flt);
    free_spinor_field_flt(out_flt);

    return return_val;
}

int test_suNg_field() {
    // Setup fields
    lprintf("TEST", 0, "Testing suNg_field assign\n");
    int return_val = 0;
    suNg_field *u_gauge_copy = alloc_suNg_field(&glattice);
    u_gauge_flt = alloc_suNg_field_flt(&glattice);

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
    copy_suNg_field_cpu(u_gauge_copy, u_gauge);
#endif
    assign_ud2u();
    assign_u2ud();

#ifdef WITH_GPU
    copy_from_gpu_suNg_field(u_gauge);
    copy_from_gpu_suNg_field(u_gauge_copy);
#endif

    sub_assign_suNg_field_cpu(u_gauge_copy, u_gauge);
    return_val += check_diff_norm(sqnorm_suNg_field_cpu(u_gauge_copy), 1.e-10);

    // Free fields
    free_suNg_field(u_gauge_copy);
    free_suNg_field_flt(u_gauge_flt);
    return return_val;
}

int test_suNf_field() {
    // Setup fields
    lprintf("TEST", 0, "Testing suNf_field assign\n");
    int return_val = 0;
    suNf_field *u_gauge_copy = alloc_suNf_field(&glattice);

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
    copy_suNf_field_cpu(u_gauge_copy, u_gauge_f);
#endif
    assign_ud2u_f();
    assign_u2ud_f();

#ifdef WITH_GPU
    copy_from_gpu_suNf_field(u_gauge_f);
    copy_from_gpu_suNf_field(u_gauge_copy);
#endif

    sub_assign_suNf_field_cpu(u_gauge_copy, u_gauge_f);
    return_val += check_diff_norm(sqnorm_suNf_field_cpu(u_gauge_copy), 1.e-10);

    // Free fields
    free_suNf_field(u_gauge_copy);
    return return_val;
}
