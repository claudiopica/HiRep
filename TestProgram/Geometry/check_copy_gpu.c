/*******************************************************************************
*
* NOCOMPILE= !WITH_GPU
*
* Check that the conversion functions from GPU to CPU format and back 
* are bijective
*
*******************************************************************************/

#define MAIN_PROGRAM

#include "suN.h"
#include "suN_types.h"
#include "setup.h"
#include "global.h"
#include "linear_algebra.h"
#include "basis_linear_algebra.h"
#include "logger.h"
#include "random.h"
#include "memory.h"
#include "update.h"
#include "geometry.h"
#include "gpu_geometry.h"
#include "hr_complex.h"

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
    lprintf("SANITY CHECK", 0, "[In field CPU copy norm unequal zero: %0.15lf]\n", spinor_field_sqnorm_f_cpu(in));

    to_gpu_format_spinor_field_f(tmp, in);
    to_cpu_format_spinor_field_f(out, tmp);
    
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
    return return_val;
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
    lprintf("SANITY CHECK", 0, "[In field CPU copy norm unequal zero: %0.15lf]\n", spinor_field_sqnorm_f_flt_cpu(in));

    // Save transformed field in CPU copy of tmp field
    to_gpu_format_spinor_field_f_flt(tmp, in);
    to_cpu_format_spinor_field_f_flt(out, tmp);

    lprintf("SANITY CHECK", 0, "[In and outfield sqnorms unequal zero and equal to each other: in %0.2e out %0.2e]\n", 
                    spinor_field_sqnorm_f_flt_cpu(in), spinor_field_sqnorm_f_flt_cpu(out));
    
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
    lprintf("INFO", 0, " ======= TEST GAUGE FIELD REPRESENTED ======= \n");
    int return_val = 0;
    int vol4h = T*X*Y*Z/2;
    suNf_field *in, *tmp, *out;
    in = alloc_gfield_f(&glattice);
    tmp = alloc_gfield_f(&glattice);
    out = alloc_gfield_f(&glattice);
    random_u_f(in);
    lprintf("SANITY CHECK", 0, "[In field CPU copy norm unequal zero: %0.15lf]\n", sqnorm_gfield_f_cpu(in));

    // Save transformed field in CPU copy of tmp field
    to_gpu_format_gfield_f(tmp, in);
    to_cpu_format_gfield_f(out, tmp);

    lprintf("SANITY CHECK", 0, "[In and outfield sqnorms unequal zero and equal to each other: in %0.2e out %0.2e]\n", 
                    sqnorm_gfield_f_cpu(in), sqnorm_gfield_f_cpu(out));

    sub_assign_gfield_f_cpu(out, in);
    double diff_norm = sqnorm_gfield_f_cpu(out);

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
    return return_val;
}

int test_convert_back_forth_gfield() 
{
    lprintf("INFO", 0, " ======= TEST GAUGE FIELD ======= \n");
    int return_val = 0;
    int vol4h = T*X*Y*Z/2;
    suNg_field *in, *tmp, *out;
    in = alloc_gfield(&glattice);
    tmp = alloc_gfield(&glattice);
    out = alloc_gfield(&glattice);
    random_u(in);
    lprintf("SANITY CHECK", 0, "[In field CPU copy norm unequal zero: %0.15lf]\n", sqnorm_gfield_cpu(in));

    to_gpu_format_gfield(tmp, in);
    to_cpu_format_gfield(out, tmp);

    lprintf("SANITY CHECK", 0, "[In and outfield sqnorms unequal zero and equal to each other: in %0.2e out %0.2e]\n", 
                    sqnorm_gfield_cpu(in), sqnorm_gfield_cpu(out));

    sub_assign_gfield_cpu(out,in);
    double diff_norm = sqnorm_gfield_cpu(out);

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
    return return_val;
}
