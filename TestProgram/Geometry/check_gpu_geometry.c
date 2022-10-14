/*******************************************************************************
*
* NOCOMPILE= !WITH_GPU
*
* Check that the GPU reading and writing functions defined in suN.h 
* are bijective.
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
#include "gpu_geometry.h"
#include "basis_linear_algebra.h"

int test_write_read_gauge_field_f();
int test_write_read_gauge_field();
int test_write_read_spinor_field_f();
int test_write_read_spinor_field_f_flt();

int main(int argc, char *argv[]) 
{
    // Init
    int return_val = 0;

    // Setup process and communication
    logger_map("DEBUG", "debug");
    setup_process(&argc, &argv);

    return_val += test_write_read_spinor_field_f();
    return_val += test_write_read_spinor_field_f_flt();
    return_val += test_write_read_gauge_field_f();
    return_val += test_write_read_gauge_field();

    // Finalize and return
    finalize_process();
    return return_val;
}

int test_write_read_gauge_field()
{
    lprintf("INFO", 0, " ======= TEST GAUGE FIELD ======= ");
    int vol4h = T*X*Y*Z/2;
    int return_val = 0;
    suNg_field *in, *gpu_format, *out;

    in = alloc_gfield(&glattice);
    out = alloc_gfield(&glattice);
    gpu_format = alloc_gfield(&glattice);
    random_u(in);
    random_u(out);

    suNg *in_mat, *block_start, *out_mat;    
    _PIECE_FOR(in->type, ixp) 
    {
        block_start = _4FIELD_BLK(gpu_format, ixp);
        _SITE_FOR(in->type, ixp, ix) 
        {
            int ix_loc = ix % vol4h;
            for (int comp = 0; comp < 4; comp++) 
            {
                in_mat = _4FIELD_AT(in, ix, comp);
                out_mat = _4FIELD_AT(out, ix, comp);
                write_gpu_suNg(vol4h, (*in_mat), block_start, ix_loc, comp);
                read_gpu_suNg(vol4h, (*out_mat), block_start, ix_loc, comp);
            }
        }
    }

    lprintf("SANITY CHECK", 0, "[Sanity check in field norm unequal zero: %0.15lf]\n", sqnorm_gfield_cpu(in));
    lprintf("SANITY CHECK", 0, "[Sanity check out field norm unequal zero: %0.15lf]\n", sqnorm_gfield_cpu(out));
    sub_assign_gfield_cpu(out, in);
    double diff_norm = sqnorm_gfield_cpu(out);

    // Since this is just a copy they have to be identical
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

    free_gfield(in);
    free_gfield(out);
    free_gfield(gpu_format);
    return return_val;
}

int test_write_read_gauge_field_f()
{
    lprintf("INFO", 0, " ======= TEST GAUGE FIELD FUNDAMENTAL REP ======= ");
    int vol4h = T*X*Y*Z/2;
    int return_val = 0;
    suNf_field *in, *gpu_format, *out;

    in = alloc_gfield_f(&glattice);
    out = alloc_gfield_f(&glattice);
    gpu_format = alloc_gfield_f(&glattice);
    random_u_f(in);

    suNf *in_mat, *block_start, *out_mat;
    _PIECE_FOR(in->type, ixp) 
    {
        block_start = _4FIELD_BLK(gpu_format, ixp);
        _SITE_FOR(in->type, ixp, ix) 
        {
            int ix_loc = ix % vol4h;
            for (int comp = 0; comp < 4; comp++) 
            {
                in_mat = _4FIELD_AT(in, ix, comp);
                out_mat = _4FIELD_AT(out, ix, comp);
                write_gpu_suNf(vol4h, (*in_mat), block_start, ix_loc, comp);
                read_gpu_suNf(vol4h, (*out_mat), block_start, ix_loc, comp);
            }
        }
    }

    lprintf("SANITY CHECK", 0, "[Sanity check in field norm unequal zero: %0.15lf]\n", sqnorm_gfield_f_cpu(in));
    lprintf("SANITY CHECK", 0, "[Sanity check out field norm unequal zero: %0.15lf]\n", sqnorm_gfield_f_cpu(out));
    sub_assign_gfield_f_cpu(out, in);
    double diff_norm = sqnorm_gfield_f_cpu(out);

    // Since this is just a copy they have to be identical
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

    free_gfield_f(in);
    free_gfield_f(out);
    free_gfield_f(gpu_format);
    return return_val;
}

int test_write_read_spinor_field_f() 
{
    lprintf("INFO", 0, " ======= TEST SPINOR FIELD ======= ");
    int vol4h = T*X*Y*Z/2;
    int return_val = 0;
    spinor_field *in, *gpu_format, *out;

    in = alloc_spinor_field_f(1, &glattice);
    gpu_format = alloc_spinor_field_f(1, &glattice);
    out = alloc_spinor_field_f(1, &glattice);
    gaussian_spinor_field(in);

    lprintf("SANITY CHECK", 0, "[Sanity check in field norm unequal zero: %0.15lf]\n", spinor_field_sqnorm_f_cpu(in));

    suNf_spinor *in_spinor, *block_start, *out_spinor;
    _PIECE_FOR(in->type, ixp) 
    {
        block_start = _FIELD_BLK(gpu_format, ixp);
        _SITE_FOR(in->type, ixp, ix) 
        {
            in_spinor = _FIELD_AT(in, ix);
            out_spinor = _FIELD_AT(out, ix);
            int ix_loc = ix % vol4h;
            for (int comp = 0; comp < 4; ++comp) 
            {
                write_gpu_suNf_spinor(vol4h, (*in_spinor).c[comp], block_start, ix_loc, comp);
                read_gpu_suNf_spinor(vol4h, (*out_spinor).c[comp], block_start, ix_loc, comp);
            }
        } 
    }

    lprintf("SANITY CHECK", 0, "[Sanity check out field norm unequal zero: %0.15lf]\n", spinor_field_sqnorm_f_cpu(out));
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

    free_spinor_field_f(in);
    free_spinor_field_f(gpu_format);
    free_spinor_field_f(out);
    return return_val;
}

int test_write_read_spinor_field_f_flt() 
{
    lprintf("INFO", 0, " ======= TEST SPINOR FIELD SINGLE PRECISION ======= ");
    int vol4h = T*X*Y*Z/2;
    int return_val = 0;
    spinor_field_flt *in, *gpu_format, *out;

    in = alloc_spinor_field_f_flt(1, &glattice);
    gpu_format = alloc_spinor_field_f_flt(1, &glattice);
    out = alloc_spinor_field_f_flt(1, &glattice);
    gaussian_spinor_field_flt(in);
    lprintf("SANITY CHECK", 0, "[Sanity check in field norm unequal zero: %0.2e]\n", spinor_field_sqnorm_f_flt_cpu(in));

    suNf_spinor_flt *in_spinor, *block_start, *out_spinor;

    _PIECE_FOR(in->type, ixp) 
    {
        block_start = _FIELD_BLK(gpu_format, ixp);
        _SITE_FOR(in->type, ixp, ix) 
        {
            in_spinor = _FIELD_AT(in, ix);
            out_spinor = _FIELD_AT(out, ix);
            int ix_loc = ix % vol4h;
            for (int comp = 0; comp < 4; ++comp) 
            {
                write_gpu_suNf_spinor_flt(vol4h, (*in_spinor).c[comp], block_start, ix_loc, comp);
                read_gpu_suNf_spinor_flt(vol4h, (*out_spinor).c[comp], block_start, ix_loc, comp);
            }
        }
    }

    lprintf("SANITY CHECK", 0, "[Sanity check out field norm unequal zero: %0.15lf]\n", spinor_field_sqnorm_f_flt_cpu(out));
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

    free_spinor_field_f_flt(in);
    free_spinor_field_f_flt(gpu_format);
    free_spinor_field_f_flt(out);
    return return_val;
}