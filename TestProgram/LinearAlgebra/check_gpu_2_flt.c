/*******************************************************************************
*
* NOCOMPILE= !WITH_GPU
*
* Test of modules (Single precision)
*
*******************************************************************************/

#include "libhr.h"

int errors = 0; // count the number of errors during this test unit

int main(int argc, char *argv[]) {
    /* setup process id and communications */
    //logger_setlevel(0,10000);
    logger_map("DEBUG", "debug");
    setup_process(&argc, &argv);

    const int niter = 1;
    // Allocate memory for CPU and GPU spinor fields
    // add 2 for the output results used in the macro TEST
    int ninputs = 3; //max number of inputs
    spinor_field_flt *in;
    in = alloc_spinor_field_flt(ninputs + 2, &glattice);

    for (int k = 0; k < niter; k++) {
        lprintf("TEST", 0, "Loop #%d\n=====================================================\n", k);
        _TEST_GPU_OP_FLT(errors, "s1=s2", 1, in, in + 1, spinor_field_copy_f_flt(out, &in[0]);
                         spinor_field_copy_f_flt_cpu(out, &in[0]););

        _TEST_RED_OP_FLT(errors, "Re<s1,s2>", 2, in, float abs1 = spinor_field_prod_re_f_flt(&in[0], &in[1]);
                         float abs2 = spinor_field_prod_re_f_flt_cpu(&in[0], &in[1]););

        _TEST_RED_OP_FLT(errors, "Im<s1,s2>", 2, in, float abs1 = spinor_field_prod_im_f_flt(&in[0], &in[1]);
                         float abs2 = spinor_field_prod_im_f_flt_cpu(&in[0], &in[1]););

        _TEST_RED_OP_FLT(errors, "<s1,s2>", 2, in, hr_complex c1 = spinor_field_prod_f_flt(&in[0], &in[1]);
                         hr_complex c2 = spinor_field_prod_f_flt_cpu(&in[0], &in[1]); float abs1 = _complex_prod_re(c1, c1);
                         float abs2 = _complex_prod_re(c2, c2););

        _TEST_RED_OP_FLT(errors, "Re<g5*s1,s2>", 2, in, float abs1 = spinor_field_g5_prod_re_f_flt(&in[0], &in[1]);
                         float abs2 = spinor_field_g5_prod_re_f_flt_cpu(&in[0], &in[1]););

        _TEST_RED_OP_FLT(errors, "Im<g5*s1,s2>", 2, in, float abs1 = spinor_field_g5_prod_im_f_flt(&in[0], &in[1]);
                         float abs2 = spinor_field_g5_prod_im_f_flt_cpu(&in[0], &in[1]););

        _TEST_RED_OP_FLT(errors, "|s1|^2", 1, in, float abs1 = spinor_field_sqnorm_f_flt(&in[0]);
                         float abs2 = spinor_field_sqnorm_f_flt_cpu(&in[0]););

        _TEST_GPU_OP_FLT(errors, "s2+=r*s1", 2, in, in + 1, float r = 10.0; spinor_field_mul_add_assign_f_flt(out, r, &in[0]);
                         spinor_field_mul_add_assign_f_flt_cpu(out, r, &in[0]););

        _TEST_GPU_OP_FLT(errors, "s2+=c*s1", 2, in, in + 1, hr_complex_flt c = 2.0 + 1.0 * I;
                         spinor_field_mulc_add_assign_f_flt(out, c, &in[0]);
                         spinor_field_mulc_add_assign_f_flt_cpu(out, c, &in[0]););

        _TEST_GPU_OP_FLT(errors, "s2+=c*g5*s1", 2, in, in + 1, hr_complex_flt c = 2.f + 3.f * I;
                         spinor_field_g5_mulc_add_assign_f_flt(out, c, &in[0]);
                         spinor_field_g5_mulc_add_assign_f_flt_cpu(out, c, &in[0]););

        _TEST_GPU_OP_FLT(errors, "s2=r*s1", 1, in, in + 1, float r = 10.0; spinor_field_mul_f_flt(out, r, &in[0]);
                         spinor_field_mul_f_flt_cpu(out, r, &in[0]););

        _TEST_GPU_OP_FLT(errors, "s2=c*s1", 1, in, in + 1, hr_complex_flt c = 2.0 + 1.0 * I;
                         spinor_field_mulc_f_flt(out, c, &in[0]); spinor_field_mulc_f_flt_cpu(out, c, &in[0]););

        _TEST_GPU_OP_FLT(errors, "s3=s1+s2", 2, in, in + 2, spinor_field_add_f_flt(out, &in[0], &in[1]);
                         spinor_field_add_f_flt_cpu(out, &in[0], &in[1]););

        _TEST_GPU_OP_FLT(errors, "s3=s1-s2", 2, in, in + 2, spinor_field_sub_f_flt(out, &in[0], &in[1]);
                         spinor_field_sub_f_flt_cpu(out, &in[0], &in[1]););

        _TEST_GPU_OP_FLT(errors, "s2+=s1", 2, in, in + 1, spinor_field_add_assign_f_flt(out, &in[0]);
                         spinor_field_add_assign_f_flt_cpu(out, &in[0]););

        _TEST_GPU_OP_FLT(errors, "s2-=s1", 2, in, in + 1, spinor_field_sub_assign_f_flt(out, &in[0]);
                         spinor_field_sub_assign_f_flt_cpu(out, &in[0]););

        _TEST_GPU_OP_FLT(errors, "s1=0", 0, in, in, spinor_field_zero_f_flt(out); spinor_field_zero_f_flt_cpu(out););

        _TEST_GPU_OP_FLT(errors, "s2=-s1", 1, in, in + 1, spinor_field_minus_f_flt(out, &in[0]);
                         spinor_field_minus_f_flt_cpu(out, &in[0]););

        _TEST_GPU_OP_FLT(errors, "s3=r*s1+s*s2", 2, in, in + 2, float r = 2.f; float s = 3.f;
                         spinor_field_lc_f_flt(out, r, &in[0], s, &in[1]);
                         spinor_field_lc_f_flt_cpu(out, r, &in[0], s, &in[1]););

        _TEST_GPU_OP_FLT(errors, "s3+=r*s1+s*s2", 3, in, in + 2, float r = 2.f; float s = 3.f;
                         spinor_field_lc_add_assign_f_flt(out, r, &in[0], s, &in[1]);
                         spinor_field_lc_add_assign_f_flt_cpu(out, r, &in[0], s, &in[1]););

        _TEST_GPU_OP_FLT(errors, "s3=c*s1+d*s2", 2, in, in + 2, hr_complex_flt c = 2.f + 3.f * I;
                         hr_complex_flt d = 3.f + 4.f * I; spinor_field_clc_f_flt(out, c, &in[0], d, &in[1]);
                         spinor_field_clc_f_flt_cpu(out, c, &in[0], d, &in[1]););

        _TEST_GPU_OP_FLT(errors, "s3+=c*s1+d*s2", 3, in, in + 2, hr_complex_flt c = 2.f + 3.f * I;
                         hr_complex_flt d = 3.f + 4.f * I; spinor_field_clc_add_assign_f_flt(out, c, &in[0], d, &in[1]);
                         spinor_field_clc_add_assign_f_flt_cpu(out, c, &in[0], d, &in[1]););

        _TEST_GPU_OP_FLT(errors, "s2=g5*s1", 1, in, in + 1, spinor_field_g5_f_flt(out, &in[0]);
                         spinor_field_g5_f_flt_cpu(out, &in[0]););

        _TEST_GPU_OP_FLT(errors, "s1=g5*s1", 1, in, in, spinor_field_g5_assign_f_flt(&in[0]);
                         spinor_field_g5_assign_f_flt_cpu(&in[0]););

        _TEST_GPU_OP_FLT(errors, "s1=g0*s2", 1, in, in + 1, spinor_field_g0_f_flt(out, &in[0]);
                         spinor_field_g0_f_flt_cpu(out, &in[0]););

        _TEST_GPU_OP_FLT(errors, "s1=g1*s2", 1, in, in + 1, spinor_field_g1_f_flt(out, &in[0]);
                         spinor_field_g1_f_flt_cpu(out, &in[0]););

        _TEST_GPU_OP_FLT(errors, "s1=g2*s2", 1, in, in + 1, spinor_field_g2_f_flt(out, &in[0]);
                         spinor_field_g2_f_flt_cpu(out, &in[0]););

        _TEST_GPU_OP_FLT(errors, "s1=g3*s2", 1, in, in + 1, spinor_field_g3_f_flt(out, &in[0]);
                         spinor_field_g3_f_flt_cpu(out, &in[0]););

        _TEST_GPU_OP_FLT(errors, "lc1", 1, in, in + 1, double c1 = 2.1; spinor_field_lc1_f_flt(c1, out, &in[0]);
                         spinor_field_lc1_f_flt_cpu(c1, out, &in[0]););

        _TEST_GPU_OP_FLT(errors, "lc2", 1, in, in + 1, double c1 = 2.4; double c2 = 4.3;
                         spinor_field_lc2_f_flt(c1, c2, out, &in[0]); spinor_field_lc2_f_flt_cpu(c1, c2, out, &in[0]););

        _TEST_GPU_OP_FLT(errors, "lc3", 3, in, in + 2, double c1 = 2.4; double c2 = 4.3;
                         spinor_field_lc3_f_flt(c1, c2, &in[0], &in[1], &in[2]);
                         spinor_field_lc3_f_flt_cpu(c1, c2, &in[0], &in[1], &in[2]););
    }

    free_spinor_field_flt(in);
    finalize_process();

    return errors;
}
