/*******************************************************************************
*
* NOCOMPILE= !WITH_GPU
* NOCOMPILE= !WITH_MPI
*
* Check that copy sync and buffer comms execute without errors
*
*******************************************************************************/

#include "libhr.h"

/* Double precision tests */
int test_sync_spinor_field();
int test_sync_gfield_f();
int test_sync_gfield();

int test_comms_spinor_field();
int test_comms_gfield_f();
int test_comms_gfield();

/* Single precision tests */
int test_sync_spinor_field_flt();
int test_sync_gfield_f_flt();
int test_sync_gfield_flt();

int test_comms_spinor_field_flt();
int test_comms_gfield_f_flt();
int test_comms_gfield_flt();

int main(int argc, char *argv[]) 
{
    // Init
    int return_val = 0;

    // Setup process and communication
    logger_map("DEBUG", "debug");
    setup_process(&argc, &argv);
    test_setup();

    // Run tests
      /* Double precision */
    return_val += test_sync_spinor_field();
    /*return_val += test_sync_gfield_f();
    return_val += test_sync_gfield();
    return_val += test_comms_spinor_field();
    return_val += test_comms_gfield_f();
    return_val += test_comms_gfield();*/

      /* Single precision */
    /*return_val += test_sync_spinor_field_flt();
    return_val += test_sync_gfield_f_flt();
    return_val += test_sync_gfield_flt();
    return_val += test_comms_spinor_field_flt();
    return_val += test_comms_gfield_f_flt();
    return_val += test_comms_gfield_flt();*/

    // Finalize and return
    finalize_process();
    return return_val;
}


int test_sync_spinor_field() 
{
    lprintf("INFO", 0, " ======= TEST SPINOR FIELD ======= \n");

    int return_val = 1;
    spinor_field *f = alloc_spinor_field_f(1, &glattice);
    gaussian_spinor_field(f);
    lprintf("SANITY CHECK", 0, "[In field CPU copy norm unequal zero: %0.2e]\n", spinor_field_sqnorm_f_cpu(f));

    sync_gpu_spinor_field_f(f);

    double sqnorm = spinor_field_sqnorm_f_cpu(f);
    lprintf("SANITY CHECK", 0, "[In field CPU copy norm unequal zero: %0.2e]\n", sqnorm);

    if (isfinite(sqnorm) && sqnorm > 1e-17) {
        return_val = 0;
    } else {
        return_val = 1;
    }

    return return_val;
}

int test_sync_gfield_f() 
{
    return 1;
}

int test_sync_gfield() 
{
    return 1;
}

int test_comms_spinor_field() 
{
    return 1;
}

int test_comms_gfield_f() 
{
    return 1;
}

int test_comms_gfield() 
{
    return 1;
}

int test_sync_spinor_field_flt() 
{
    return 1;
}

int test_sync_gfield_f_flt() 
{
    return 1;
}

int test_sync_gfield_flt() 
{
    return 1;
}

int test_comms_spinor_field_flt() 
{
    return 1;
}

int test_comms_gfield_f_flt() 
{
    return 1;
}

int test_comms_gfield_flt() 
{
    return 1;
}



