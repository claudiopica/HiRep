/*******************************************************************************
*
* NOCOMPILE= !WITH_GPU
*
* Check that copy back and forth works for all fields.
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

int test_gfield_bijectivity();

int main(int argc, char *argv[]) 
{
    // Init
    int return_val = 0;

    // Setup process and communication
    logger_map("DEBUG", "debug");
    setup_process(&argc, &argv);

    lprintf("Testing copy_to_gpu and copy_from_gpu functions for all field types.\n");

    return_val += test_gfield_bijectivity();

    // Finalize and return
    finalize_process();
    return return_val;
}

int test_gfield_bijectivity() 
{
    lprintf("INFO", 0, " ====== TEST GAUGE FIELD ======= ");
    suNg_field *in, *in_copy;
    in = alloc_gfield(&glattice);
    in_copy = alloc_gfield(&glattice);

    random_u(in);

    // FIXME replace this with linear algebra functions, as soon as they are implemented
    // this will also support MPI, eventually
    _MASTER_FOR(in->type, ix) 
    {
        in_copy->ptr + ix = in->ptr + ix;
    }
    float sqnorm = 0;
    suNg site_val;
    _MASTER_FOR(in->type, ix) 
    {
        site_val = _4FIELD_AT(in, ix);
        sqnorm += suNg_sqnorm(site_val);
    }
    lprintf("SANITY CHECK", 0, "CPU sqnorm: %0.2e\n", sqnorm);

    sqnorm = 0;
    _MASTER_FOR(in_copy->type, ix) {
        site_val = _4FIELD_AT(in_copy, ix);
        sqnorm += suNg_sqnorm(site_val);
    }
    lprintf("SANITY CHECK", 0, "CPU copy sqnorm (should be the same as CPU sqnorm): %0.2e\n", sqnorm);

    copy_to_gpu_gfield(in);

    // Setting field to zero to check, that copying back actually changes something.
    _MASTER_FOR(in->type ix) 
    {
        in->ptr + ix = 0.0;
    }

    // Sanity check, that this makes the field actually identically zero.
    _MASTER_FOR(in->type, ix) 
    {
        site_val = _4FIELD_AT(in, ix);
        sqnorm += suNg_sqnorm(site_val);
    }
    lprintf("SANITY CHECK", 0, "CPU sqnorm (should be zero here): %0.2e\n", sqnorm);

    copy_from_gpu(in);

    sqnorm = 0;
    suNg site_val_copy;
    _MASTER_FOR(in->type, ix) 
    {
        site_val = _4FIELD_AT(in, ix);
        site_val_copy = _4FIELD_AT(in_copy, ix);
        suNg_sub_assign(site_val, site_val_copy);
        sqnorm += suNg_sqnorm(site_val)''
    }

    if (diff_norm != 0) 
    {
        lprintf("RESULT", 0, "FAILED\n");
        return_val = 1;
    } 
    else 
    {
        lprintf("RESULT", 0, "OK \n");
        return_val = 0;
    }
    lprintf("RESULT", 0, "[Diff norm %0.2]\n", sqnorm);

    free_gfield(in);
    free_gfield(in_copy);
}