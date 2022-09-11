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

int main(int argc, char *argv[]) 
{
    // Init
    int return_val = 0;

    // Setup process and communication
    logger_map("DEBUG", "debug");
    setup_process(&argc, &argv);

    //return_val += test();

    // Finalize and return
    finalize_process();
    return return_val;
}

int test_convert_back_forth_spinor_field() 
{
    lprintf("INFO", 0, " ======= TEST SPINOR FIELD ======= \n");
    int return_val = 0;
    spinor_field *in, *out;
    in = alloc_spinor_field_f(1, &glattice);
    out = alloc_spinor_field_f(1, &glattice);
    gaussian_spinor_field(in);

    lprintf("SANITY CHECK", 0, "[In field norm unequal zero: %0.15lf]\n", spinor_field_sqnorm_f_cpu(in));

    // Copy in to out
    spinor_field_mul_f(out, 1.0, in);
    lprintf("SANITY CHECK", 0, "[In field and out field identical norms: in: %0.15lf out: %0.15lf]\n", 
                spinor_field_sqnorm_f_cpu(in), spinor_field_sqnorm_f_cpu(out));

    // Copy out field to gpu
    togpuformat_spinor_field_f(out, in);
    lprintf("SANITY CHECK", 0, "[Out field GPU copy not zero: %0.15lf]\n", spinor_field_sqnorm_f(out));

    // Overwrite cpu copy of out field
    spinor_field_zero_f_cpu(out);
    lprintf("SANITY CHECK", 0, "[Out field CPU copy overwritten to zero: %0.15lf]\n", spinor_field_sqnorm_f_cpu(out));

    // Copy back and check that its the same as the infield
    //fromgpuformat_spinor_field_f(in);
    
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