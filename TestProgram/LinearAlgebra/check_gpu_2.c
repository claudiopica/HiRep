/*******************************************************************************
*
* NOCOMPILE= !WITH_GPU
*
* Test of modules
*
*******************************************************************************/

#include "libhr.h"

int errors = 0; // count the number of errors during this test uni

int main(int argc, char *argv[]) {
    /* setup process id and communications */
    //logger_setlevel(0,10000);
    logger_map("DEBUG", "debug");
    setup_process(&argc, &argv);

    const int niter = 1;
    // Allocate memory for CPU and GPU spinor fields
    // add 2 for the output results used in the macro TEST
    int ninputs = 3; //max number of inputs
    spinor_field *in;
    in = alloc_spinor_field_f(ninputs + 2, &glattice);

    for (int k = 0; k < niter; k++) {
        lprintf("TEST", 0, "Loop #%d\n=====================================================\n", k);
        _TEST_GPU_OP(errors, "s1=s2", 1, in, in + 1, spinor_field_copy_f(out, &in[0]); spinor_field_copy_f_cpu(out, &in[0]););

        _TEST_RED_OP(errors, "Re<s1,s2>", 2, in, double abs1 = spinor_field_prod_re_f(&in[0], &in[1]);
                     double abs2 = spinor_field_prod_re_f_cpu(&in[0], &in[1]););

        _TEST_RED_OP(errors, "Im<s1,s2>", 2, in, double abs1 = spinor_field_prod_im_f(&in[0], &in[1]);
                     double abs2 = spinor_field_prod_im_f_cpu(&in[0], &in[1]););

        _TEST_RED_OP(errors, "<s1,s2>", 2, in, hr_complex c1 = spinor_field_prod_f(&in[0], &in[1]);
                     hr_complex c2 = spinor_field_prod_f_cpu(&in[0], &in[1]); double abs1 = _complex_prod_re(c1, c1);
                     double abs2 = _complex_prod_re(c2, c2););

        _TEST_RED_OP(errors, "Re<g5*s1,s2>", 2, in, double abs1 = spinor_field_g5_prod_re_f(&in[0], &in[1]);
                     double abs2 = spinor_field_g5_prod_re_f_cpu(&in[0], &in[1]););

        _TEST_RED_OP(errors, "Im<g5*s1,s2>", 2, in, double abs1 = spinor_field_g5_prod_im_f(&in[0], &in[1]);
                     double abs2 = spinor_field_g5_prod_im_f_cpu(&in[0], &in[1]););

        _TEST_RED_OP(errors, "|s1|^2", 1, in, double abs1 = spinor_field_sqnorm_f(&in[0]);
                     double abs2 = spinor_field_sqnorm_f_cpu(&in[0]););

        _TEST_GPU_OP(errors, "s2+=r*s1", 2, in, in + 1, double r = 10.0; spinor_field_mul_add_assign_f(out, r, &in[0]);
                     spinor_field_mul_add_assign_f_cpu(out, r, &in[0]););

        _TEST_GPU_OP(errors, "s2+=c*s1", 2, in, in + 1, hr_complex c = 2.0 + 1.0 * I;
                     spinor_field_mulc_add_assign_f(out, c, &in[0]); spinor_field_mulc_add_assign_f_cpu(out, c, &in[0]););

        _TEST_GPU_OP(errors, "s2+=c*g5*s1", 2, in, in + 1, hr_complex c = 2.0 + 3.0 * I;
                     spinor_field_g5_mulc_add_assign_f(out, c, &in[0]); spinor_field_g5_mulc_add_assign_f_cpu(out, c, &in[0]););

        _TEST_GPU_OP(errors, "s2=r*s1", 1, in, in + 1, double r = 10.0; spinor_field_mul_f(out, r, &in[0]);
                     spinor_field_mul_f_cpu(out, r, &in[0]););

        _TEST_GPU_OP(errors, "s2=c*s1", 1, in, in + 1, hr_complex c = 2.0 + 1.0 * I; spinor_field_mulc_f(out, c, &in[0]);
                     spinor_field_mulc_f_cpu(out, c, &in[0]););

        _TEST_GPU_OP(errors, "s3=s1+s2", 2, in, in + 2, spinor_field_add_f(out, &in[0], &in[1]);
                     spinor_field_add_f_cpu(out, &in[0], &in[1]););

        _TEST_GPU_OP(errors, "s3=s1-s2", 2, in, in + 2, spinor_field_sub_f(out, &in[0], &in[1]);
                     spinor_field_sub_f_cpu(out, &in[0], &in[1]););

        _TEST_GPU_OP(errors, "s2+=s1", 2, in, in + 1, spinor_field_add_assign_f(out, &in[0]);
                     spinor_field_add_assign_f_cpu(out, &in[0]););

        _TEST_GPU_OP(errors, "s2-=s1", 2, in, in + 1, spinor_field_sub_assign_f(out, &in[0]);
                     spinor_field_sub_assign_f_cpu(out, &in[0]););

        _TEST_GPU_OP(errors, "s1=0", 0, in, in, spinor_field_zero_f(out); spinor_field_zero_f_cpu(out););

        _TEST_GPU_OP(errors, "s2=-s1", 1, in, in + 1, spinor_field_minus_f(out, &in[0]);
                     spinor_field_minus_f_cpu(out, &in[0]););

        _TEST_GPU_OP(errors, "s3=r*s1+s*s2", 2, in, in + 2, double r = 2.0; double s = 3.0;
                     spinor_field_lc_f(out, r, &in[0], s, &in[1]); spinor_field_lc_f_cpu(out, r, &in[0], s, &in[1]););

        _TEST_GPU_OP(errors, "s3+=r*s1+s*s2", 3, in, in + 2, double r = 2.0; double s = 3.0;
                     spinor_field_lc_add_assign_f(out, r, &in[0], s, &in[1]);
                     spinor_field_lc_add_assign_f_cpu(out, r, &in[0], s, &in[1]););

        _TEST_GPU_OP(errors, "s3=c*s1+d*s2", 2, in, in + 2, hr_complex c = 2.0 + 3.0 * I; hr_complex d = 3.0 + 4.0 * I;
                     spinor_field_clc_f(out, c, &in[0], d, &in[1]); spinor_field_clc_f_cpu(out, c, &in[0], d, &in[1]););

        _TEST_GPU_OP(errors, "s3+=c*s1+d*s2", 3, in, in + 2, hr_complex c = 2.0 + 3.0 * I; hr_complex d = 3.0 + 4.0 * I;
                     spinor_field_clc_add_assign_f(out, c, &in[0], d, &in[1]);
                     spinor_field_clc_add_assign_f_cpu(out, c, &in[0], d, &in[1]););

        _TEST_GPU_OP(errors, "s2=g5*s1", 1, in, in + 1, spinor_field_g5_f(out, &in[0]); spinor_field_g5_f_cpu(out, &in[0]););

        _TEST_GPU_OP(errors, "s1=g5*s1", 1, in, in, spinor_field_g5_assign_f(&in[0]); spinor_field_g5_assign_f_cpu(&in[0]););

        _TEST_GPU_OP(errors, "s1=g0*s2", 1, in, in + 1, spinor_field_g0_f(out, &in[0]); spinor_field_g0_f_cpu(out, &in[0]););

        _TEST_GPU_OP(errors, "s1=g1*s2", 1, in, in + 1, spinor_field_g1_f(out, &in[0]); spinor_field_g1_f_cpu(out, &in[0]););

        _TEST_GPU_OP(errors, "s1=g2*s2", 1, in, in + 1, spinor_field_g2_f(out, &in[0]); spinor_field_g2_f_cpu(out, &in[0]););

        _TEST_GPU_OP(errors, "s1=g3*s2", 1, in, in + 1, spinor_field_g3_f(out, &in[0]); spinor_field_g3_f_cpu(out, &in[0]););

        _TEST_GPU_OP(errors, "lc1", 1, in, in + 1, double c1 = 2.1; spinor_field_lc1_f(c1, out, &in[0]);
                     spinor_field_lc1_f_cpu(c1, out, &in[0]););

        _TEST_GPU_OP(errors, "lc2", 1, in, in + 1, double c1 = 2.4; double c2 = 4.3; spinor_field_lc2_f(c1, c2, out, &in[0]);
                     spinor_field_lc2_f_cpu(c1, c2, out, &in[0]););

        _TEST_GPU_OP(errors, "lc3", 3, in, in + 2, double c1 = 2.4; double c2 = 4.3;
                     spinor_field_lc3_f(c1, c2, &in[0], &in[1], &in[2]);
                     spinor_field_lc3_f_cpu(c1, c2, &in[0], &in[1], &in[2]););
    }

    free_spinor_field_f(in);
    finalize_process();

    return errors;
}
