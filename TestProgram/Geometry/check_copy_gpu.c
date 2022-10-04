/*******************************************************************************
*
* NOCOMPILE= !WITH_GPU
*
* Check that the conversion functions ..._togpuformat and ..._tocpuformat
* are bijective
*
*******************************************************************************/

#define MAIN_PROGRAM

#include "suN.h"
#include "suN_types.h"
#include "setup.h"
#include "global.h"
#include "linear_algebra.h"
#include "logger.h"
#include "random.h"
#include "memory.h"
#include "update.h"
#include "geometry.h"

int test_convert_back_forth_spinor_field();
int test_convert_back_forth_spinor_field_flt();
int test_convert_back_forth_gfield_f();
int test_convert_back_forth_gfield();

int main(int argc, char *argv[])
{
    // Init
    int return_val = 0;

    // Setup process and communication
    logger_map("DEBUG", "debug");
    setup_process(&argc, &argv);

    // Run tests
    return_val += test_convert_back_forth_spinor_field();
    return_val += test_convert_back_forth_spinor_field_flt();
    return_val += test_convert_back_forth_gfield_f();
    return_val += test_convert_back_forth_gfield();

    // Finalize and return
    finalize_process();
    return return_val;
}

int test_convert_back_forth_spinor_field()
{
    lprintf("INFO", 0, " ======= TEST SPINOR FIELD ======= \n");
    int return_val = 0;
    spinor_field *in, *tmp, *out;
    in = alloc_spinor_field_f(1, &glattice);
    tmp = alloc_spinor_field_f(1, &glattice);
    out = alloc_spinor_field_f(1, &glattice);
    gaussian_spinor_field(in);

    // Save transformed field in CPU copy of tmp field
    spinor_field_togpuformat(tmp, in);

    // Sanity checks that the CPU copy of in field
    // and CPU copy of the tmp field have non-zero square norms
    lprintf("SANITY CHECK", 0, "[In field CPU copy norm unequal zero: %0.15lf]\n", spinor_field_sqnorm_f_cpu(in));
    lprintf("SANITY CHECK", 0, "[Tmp field CPU copy norm unequal zero: %0.15lf]\n", spinor_field_sqnorm_f_cpu(tmp));

    // Transform back to out field
    spinor_field_tocpuformat(out, tmp);

    spinor_field_sub_assign_f_cpu(out, in);
    double diff_norm = spinor_field_sqnorm_f_cpu(out);

    if (diff_norm != 0)
    {
        lprintf("RESULT", 0, "FAILED \n");
        return_val = 1;
    }
    else
    {
        lprintf("RESULT", 0, "OK \n");
        return_val = 0;
    }
    lprintf("RESULT", 0, "[Diff norm %0.2e]\n", diff_norm);
}

int test_convert_back_forth_spinor_field_flt()
{
    lprintf("INFO", 0, " ======= TEST SPINOR FIELD SINGLE PRECISION ======= \n");
    int return_val = 0;
    spinor_field_flt *in, *tmp, *out;
    in = alloc_spinor_field_f_flt(1, &glattice);
    tmp = alloc_spinor_field_f_flt(1, &glattice);
    out = alloc_spinor_field_f_flt(1, &glattice);
    gaussian_spinor_field_flt(in);

    // Save transformed field in CPU copy of tmp field
    spinor_field_togpuformat_flt(tmp, in);

    // Sanity checks that the CPU copy of in field
    // and CPU copy of the tmp field have non-zero square norms
    lprintf("SANITY CHECK", 0, "[In field CPU copy norm unequal zero: %0.15lf]\n", spinor_field_sqnorm_f_flt_cpu(in));
    lprintf("SANITY CHECK", 0, "[Tmp field CPU copy norm unequal zero: %0.15lf]\n", spinor_field_sqnorm_f_flt_cpu(tmp));

    // Transform back to out field
    spinor_field_tocpuformat_flt(out, tmp);

    spinor_field_sub_assign_f_flt_cpu(out, in);
    double diff_norm = spinor_field_sqnorm_f_flt_cpu(out);

    if (diff_norm != 0)
    {
        lprintf("RESULT", 0, "FAILED \n");
        return_val = 1;
    }
    else
    {
        lprintf("RESULT", 0, "OK \n");
        return_val = 0;
    }
    lprintf("RESULT", 0, "[Diff norm %0.2e]\n", diff_norm);
}

int test_convert_back_forth_gfield_f()
{
    lprintf("INFO", 0, " ======= TEST GAUGE FIELD IN FUNDAMENTAL REP ======= \n");
    int return_val = 0;
    suNf_field *in, *tmp, *out;
    in = alloc_gfield_f(&glattice);
    tmp = alloc_gfield_f(&glattice);
    out = alloc_gfield_f(&glattice);
    random_u_f(in);

    // Save transformed field in CPU copy of tmp field
    gfield_togpuformat_f(tmp, in);

    // Transform back to out field
    gfield_tocpuformat_f(out, tmp);

    suNf in_mat, tmp_mat, out_mat;
    double sqnorm = 0.0;
    double sqnorm_in_check = 0.0;
    double sqnorm_tmp_check = 0.0;
    double sqnorm_out_check = 0.0;
    double diff_norm = 0.0;

    _MASTER_FOR(in->type, ix)
    {
        in_mat = *(in->ptr+ix);
        out_mat = *(out->ptr+ix);
        tmp_mat = *(tmp->ptr+ix);

        _suNf_sqnorm(sqnorm, in_mat);
        sqnorm_in_check += sqnorm;

        _suNf_sqnorm(sqnorm, tmp_mat);
        sqnorm_tmp_check += sqnorm;

        _suNf_sqnorm(sqnorm, out_mat);
        sqnorm_out_check += sqnorm;

        _suNg_sub_assign(out_mat, in_mat);
        _suNg_sqnorm(sqnorm, out_mat);
        diff_norm += sqnorm;
    }

    lprintf("SANITY CHECK", 0, "[Tmp sqnorm unequal zero: %0.2e]\n", sqnorm_tmp_check);
    lprintf("SANITY CHECK", 0, "[In and outfield sqnorms unequal zero and equal to each other: in %0.2e out %0.2e]\n",
                    sqnorm_in_check, sqnorm_out_check);

    if (diff_norm != 0)
    {
        lprintf("RESULT", 0, "FAILED \n");
        return_val = 1;
    }
    else
    {
        lprintf("RESULT", 0, "OK \n");
        return_val = 0;
    }
    lprintf("RESULT", 0, "[Diff norm %0.2e]\n", diff_norm);
}

int test_convert_back_forth_gfield()
{
    lprintf("INFO", 0, " ======= TEST GAUGE FIELD ======= \n");
    int return_val = 0;
    suNg_field *in, *tmp, *out;
    in = alloc_gfield(&glattice);
    tmp = alloc_gfield(&glattice);
    out = alloc_gfield(&glattice);
    random_u(in);

    // Save transformed field in CPU copy of tmp field
    gfield_togpuformat(tmp, in);

    // Transform back to out field
    gfield_tocpuformat(out, tmp);

    suNg in_mat, tmp_mat, out_mat;
    double sqnorm = 0.0;
    double sqnorm_in_check = 0.0;
    double sqnorm_tmp_check = 0.0;
    double sqnorm_out_check = 0.0;
    double diff_norm = 0.0;

    _MASTER_FOR(in->type, ix)
    {
        in_mat = *(in->ptr+ix);
        out_mat = *(out->ptr+ix);
        tmp_mat = *(tmp->ptr+ix);

        _suNf_sqnorm(sqnorm, in_mat);
        sqnorm_in_check += sqnorm;

        _suNf_sqnorm(sqnorm, tmp_mat);
        sqnorm_tmp_check += sqnorm;

        _suNf_sqnorm(sqnorm, out_mat);
        sqnorm_out_check += sqnorm;

        _suNg_sub_assign(out_mat, in_mat);
        _suNg_sqnorm(sqnorm, out_mat);
        diff_norm += sqnorm;
    }

    lprintf("SANITY CHECK", 0, "[Tmp sqnorm unequal zero: %0.2e]\n", sqnorm_tmp_check);
    lprintf("SANITY CHECK", 0, "[In and outfield sqnorms unequal zero and equal to each other: in %0.2e out %0.2e]\n",
                    sqnorm_in_check, sqnorm_out_check);

    if (diff_norm != 0)
    {
        lprintf("RESULT", 0, "FAILED \n");
        return_val = 1;
    }
    else
    {
        lprintf("RESULT", 0, "OK \n");
        return_val = 0;
    }
    lprintf("RESULT", 0, "[Diff norm %0.2e]\n", diff_norm);
}
