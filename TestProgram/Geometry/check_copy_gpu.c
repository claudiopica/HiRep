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

int main(int argc, char *argv[]) 
{
    // Init
    int return_val = 0;

    // Setup process and communication
    logger_map("DEBUG", "debug");
    setup_process(&argc, &argv);

    // Run tests
    return_val += test_convert_back_forth_spinor_field();

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
    togpuformat_spinor_field_f(tmp, in);

    // Sanity checks that the CPU copy of in field 
    // and CPU copy of the tmp field have non-zero square norms
    lprintf("SANITY CHECK", 0, "[In field CPU copy norm unequal zero: %0.15lf]\n", spinor_field_sqnorm_f_cpu(in));
    lprintf("SANITY CHECK", 0, "[Tmp field CPU copy norm unequal zero: %0.15lf]\n", spinor_field_sqnorm_f_cpu(tmp));

    // Transform back to out field
    //tocpuformat_spinor_field_f(out, tmp);
    
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