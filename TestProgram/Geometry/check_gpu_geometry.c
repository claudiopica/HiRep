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
#include <stdio.h>
#include "mpi.h"
#include "hr_complex.h"

// TODO: get gaussian spinor fields to work with MPI-CUDA
// TODO: &glat_even, &glat_odd
// TODO: Test for n > 1
// TODO: Get scalar fields to pass

// Double precision
int test_write_read_gauge_field_f();
int test_write_read_gauge_field();
int test_write_read_spinor_field_f();
int test_write_read_spinor_field_f_vector_wise();
int test_write_read_sfield();
int test_write_read_scalar_field();

// Single precision
int test_write_read_gauge_field_flt();
int test_write_read_spinor_field_f_flt();
int test_write_read_spinor_field_f_flt_vector_wise();

int main(int argc, char *argv[]) 
{
    // Init
    int return_val = 0;

    // Setup process and communication
    logger_map("DEBUG", "debug");
    setup_process(&argc, &argv);
    test_setup();

    // Double Precision Tests
    return_val += test_write_read_gauge_field_f();
    return_val += test_write_read_gauge_field();
    return_val += test_write_read_spinor_field_f();
    return_val += test_write_read_spinor_field_f_vector_wise();
    //return_val += test_write_read_sfield();
    //return_val += test_write_read_scalar_field();

    //Single Precision Tests
    return_val += test_write_read_gauge_field_flt();
    //return_val += test_write_read_spinor_field_f_flt();
    return_val += test_write_read_spinor_field_f_flt_vector_wise();


    // Finalize and return
    finalize_process();
    return return_val;
}

int test_write_read_gauge_field()
{
    lprintf("INFO", 0, " ======= TEST GAUGE FIELD ======= ");
    int return_val = 0;
    suNg_field *in, *gpu_format, *out;

    in = alloc_gfield(&glattice);
    out = alloc_gfield(&glattice);
    gpu_format = alloc_gfield(&glattice);

    random_u(in); 
    lprintf("SANITY CHECK", 0, "[In field norm unequal zero: %0.2e]\n", sqnorm_gfield_cpu(in));

    suNg *in_mat, *block_start, *out_mat;   
    int stride = 0; 
    _PIECE_FOR(in->type, ixp) 
    {
        block_start = _4FIELD_BLK(gpu_format, ixp);
        stride = in->type->master_end[ixp] - in->type->master_end[ixp] +1;
        _SITE_FOR(in->type, ixp, ix) 
        {
            int ix_loc = _GPU_IDX_TO_LOCAL(in, ix, ixp);
            for (int comp = 0; comp < 4; comp++) 
            {
                in_mat = _4FIELD_AT(in, ix, comp);
                out_mat = _4FIELD_AT(out, ix, comp);
                write_gpu_suNg(stride, (*in_mat), block_start, ix_loc, comp);
                read_gpu_suNg(stride, (*out_mat), block_start, ix_loc, comp);
            }
        }
    }

    lprintf("SANITY CHECK", 0, "[Sanity check in field norm unequal zero: %0.2e]\n", sqnorm_gfield_cpu(in));
    lprintf("SANITY CHECK", 0, "[Sanity check out field norm unequal zero: %0.2e]\n", sqnorm_gfield_cpu(out));
    sub_assign_gfield_cpu(out, in);
    double diff_norm = sqnorm_gfield_cpu(out);

    check_diff_norm_zero(diff_norm);

    free_gfield(in);
    free_gfield(out);
    free_gfield(gpu_format);
    return return_val;
}

int test_write_read_gauge_field_flt()
{
    lprintf("INFO", 0, " ======= TEST GAUGE FIELD SINGLE PRECISION ======= ");
    int return_val = 0;
    suNg_field_flt *in, *gpu_format, *out;

    in = alloc_gfield_flt(&glattice);
    out = alloc_gfield_flt(&glattice);
    gpu_format = alloc_gfield_flt(&glattice);

    random_gfield_flt_cpu(in);
    lprintf("SANITY CHECK", 0, "[In field norm unequal zero: %0.2e]\n", sqnorm_gfield_flt_cpu(in));

    suNg_flt *in_mat, *block_start, *out_mat;  
    int stride = 0;  
    _PIECE_FOR(in->type, ixp) 
    {
        block_start = _4FIELD_BLK(gpu_format, ixp);
        stride = in->type->master_end[ixp] - in->type->master_start[ixp] +1;
        _SITE_FOR(in->type, ixp, ix) 
        {
            int ix_loc = _GPU_IDX_TO_LOCAL(in, ix, ixp);
            for (int comp = 0; comp < 4; comp++) 
            {
                in_mat = _4FIELD_AT(in, ix, comp);
                out_mat = _4FIELD_AT(out, ix, comp);
                write_gpu_suNg_flt(stride, (*in_mat), block_start, ix_loc, comp);
                read_gpu_suNg_flt(stride, (*out_mat), block_start, ix_loc, comp);
            }
        }
    }

    lprintf("SANITY CHECK", 0, "[Sanity check in field norm unequal zero: %0.2e]\n", sqnorm_gfield_flt_cpu(in));
    lprintf("SANITY CHECK", 0, "[Sanity check out field norm unequal zero: %0.2e]\n", sqnorm_gfield_flt_cpu(out));
    sub_assign_gfield_flt_cpu(out, in);
    double diff_norm = sqnorm_gfield_flt_cpu(out);
    check_diff_norm_zero(diff_norm);

    free_gfield_flt(in);
    free_gfield_flt(out);
    free_gfield_flt(gpu_format);
    return return_val;
}

int test_write_read_gauge_field_f()
{
    lprintf("INFO", 0, " ======= TEST GAUGE FIELD FUNDAMENTAL REP ======= ");
    int return_val = 0;
    suNf_field *in, *gpu_format, *out;

    in = alloc_gfield_f(&glattice);
    out = alloc_gfield_f(&glattice);
    gpu_format = alloc_gfield_f(&glattice);

    random_u_f(in);
    lprintf("SANITY CHECK", 0, "[In field norm unequal zero: %0.2e]\n", sqnorm_gfield_f_cpu(in));

    suNf *in_mat, *block_start, *out_mat;
    int stride = 0;
    _PIECE_FOR(in->type, ixp) 
    {
        block_start = _4FIELD_BLK(gpu_format, ixp);
        stride = in->type->master_end[ixp] - in->type->master_start[ixp] + 1;
        _SITE_FOR(in->type, ixp, ix) 
        {
            int ix_loc = _GPU_IDX_TO_LOCAL(in, ix, ixp);
            for (int comp = 0; comp < 4; comp++) 
            {
                in_mat = _4FIELD_AT(in, ix, comp);
                out_mat = _4FIELD_AT(out, ix, comp);
                write_gpu_suNf(stride, (*in_mat), block_start, ix_loc, comp);
                read_gpu_suNf(stride, (*out_mat), block_start, ix_loc, comp);
            }
        }
    }

    lprintf("SANITY CHECK", 0, "[Sanity check in field norm unequal zero: %0.2e]\n", sqnorm_gfield_f_cpu(in));
    lprintf("SANITY CHECK", 0, "[Sanity check out field norm unequal zero: %0.2e]\n", sqnorm_gfield_f_cpu(out));
    sub_assign_gfield_f_cpu(out, in);
    double diff_norm = sqnorm_gfield_f_cpu(out);

    check_diff_norm_zero(diff_norm);

    free_gfield_f(in);
    free_gfield_f(out);
    free_gfield_f(gpu_format);
    return return_val;
}

int test_write_read_spinor_field_f() 
{
    lprintf("INFO", 0, " ======= TEST SPINOR FIELD ======= ");
    int return_val = 0;
    spinor_field *in, *gpu_format, *out;

    in = alloc_spinor_field_f(1, &glattice);
    gpu_format = alloc_spinor_field_f(1, &glattice);
    out = alloc_spinor_field_f(1, &glattice);

    //gaussian_spinor_field(in);
    random_spinor_field_f_cpu(in);

    lprintf("SANITY CHECK", 0, "[Sanity check in field norm unequal zero: %0.2e]\n", spinor_field_sqnorm_f_cpu(in));

    suNf_spinor *in_spinor, *block_start, *out_spinor;
    int stride = 0;
    _PIECE_FOR(in->type, ixp)
    {
        block_start = _FIELD_BLK(gpu_format, ixp);
        stride = in->type->master_end[ixp] - in->type->master_start[ixp] + 1;
        _SITE_FOR(in->type, ixp, ix) 
        {
            in_spinor = _FIELD_AT(in, ix);
            out_spinor = _FIELD_AT(out, ix);
            int ix_loc = _GPU_IDX_TO_LOCAL(in, ix, ixp);
            write_gpu_suNf_spinor(stride, (*in_spinor), block_start, ix_loc, 0);
            read_gpu_suNf_spinor(stride, (*out_spinor), block_start, ix_loc, 0);
        } 
    }

    lprintf("SANITY CHECK", 0, "[Sanity check out field norm unequal zero: %0.2e]\n", spinor_field_sqnorm_f_cpu(out));
    spinor_field_sub_assign_f_cpu(out, in);
    double diff_norm = spinor_field_sqnorm_f_cpu(out);
    check_diff_norm_zero(diff_norm);

    free_spinor_field_f(in);
    free_spinor_field_f(gpu_format);
    free_spinor_field_f(out);
    return return_val;
}

int test_write_read_spinor_field_f_vector_wise() 
{
    lprintf("INFO", 0, " ======= TEST SPINOR FIELD II ======= ");
    int return_val = 0;
    spinor_field *in, *gpu_format, *out;

    in = alloc_spinor_field_f(1, &glattice);
    gpu_format = alloc_spinor_field_f(1, &glattice);
    out = alloc_spinor_field_f(1, &glattice);

    //gaussian_spinor_field(in);
    random_spinor_field_f_cpu(in);

    lprintf("SANITY CHECK", 0, "[Sanity check in field norm unequal zero: %0.2e]\n", spinor_field_sqnorm_f_cpu(in));

    suNf_spinor *in_spinor, *block_start, *out_spinor;
    int stride = 0;
    _PIECE_FOR(in->type, ixp)
    {
        block_start = _FIELD_BLK(gpu_format, ixp);
        stride = in->type->master_end[ixp] - in->type->master_start[ixp] + 1;
        _SITE_FOR(in->type, ixp, ix) 
        {
            in_spinor = _FIELD_AT(in, ix);
            out_spinor = _FIELD_AT(out, ix);
            int ix_loc = _GPU_IDX_TO_LOCAL(in, ix, ixp);
            for (int comp = 0; comp < 4; ++comp) 
            {
                write_gpu_suNf_vector(stride, (*in_spinor).c[comp], block_start, ix_loc, comp);
                read_gpu_suNf_vector(stride, (*out_spinor).c[comp], block_start, ix_loc, comp);
            }
        } 
    }

    lprintf("SANITY CHECK", 0, "[Sanity check out field norm unequal zero: %0.2e]\n", spinor_field_sqnorm_f_cpu(out));
    spinor_field_sub_assign_f_cpu(out, in);
    double diff_norm = spinor_field_sqnorm_f_cpu(out);
    check_diff_norm_zero(diff_norm);

    free_spinor_field_f(in);
    free_spinor_field_f(gpu_format);
    free_spinor_field_f(out);
    return return_val;
}

int test_write_read_spinor_field_f_flt() 
{
    lprintf("INFO", 0, " ======= TEST SPINOR FIELD SINGLE PRECISION ======= ");
    int return_val = 0;
    spinor_field_flt *in, *gpu_format, *out;

    in = alloc_spinor_field_f_flt(1, &glattice);
    gpu_format = alloc_spinor_field_f_flt(1, &glattice);
    out = alloc_spinor_field_f_flt(1, &glattice);

    //gaussian_spinor_field_flt(in);
    random_spinor_field_f_flt_cpu(in);

    lprintf("SANITY CHECK", 0, "[Sanity check in field norm unequal zero: %0.2e]\n", spinor_field_sqnorm_f_flt_cpu(in));

    suNf_spinor_flt *in_spinor, *block_start, *out_spinor;
    int stride = 0;
    _PIECE_FOR(in->type, ixp) 
    {
        block_start = _FIELD_BLK(gpu_format, ixp);
        stride = in->type->master_end[ixp] - in->type->master_start[ixp] + 1;
        _SITE_FOR(in->type, ixp, ix) 
        {
            in_spinor = _FIELD_AT(in, ix);
            out_spinor = _FIELD_AT(out, ix);
            int ix_loc = _GPU_IDX_TO_LOCAL(in, ix, ixp);
            write_gpu_suNf_spinor_flt(stride, (*in_spinor), block_start, ix_loc, 0);
            read_gpu_suNf_spinor_flt(stride, (*out_spinor), block_start, ix_loc, 0);
        }
    }

    lprintf("SANITY CHECK", 0, "[Sanity check out field norm unequal zero: %0.2e]\n", spinor_field_sqnorm_f_flt_cpu(out));
    spinor_field_sub_assign_f_flt_cpu(out, in);
    double diff_norm = spinor_field_sqnorm_f_flt_cpu(out);
    check_diff_norm_zero(diff_norm);

    free_spinor_field_f_flt(in);
    free_spinor_field_f_flt(gpu_format);
    free_spinor_field_f_flt(out);
    return return_val;
}

int test_write_read_spinor_field_f_flt_vector_wise() 
{
    lprintf("INFO", 0, " ======= TEST SPINOR FIELD SINGLE PRECISION II ======= ");
    int return_val = 0;
    spinor_field_flt *in, *gpu_format, *out;

    in = alloc_spinor_field_f_flt(1, &glattice);
    gpu_format = alloc_spinor_field_f_flt(1, &glattice);
    out = alloc_spinor_field_f_flt(1, &glattice);

    gaussian_spinor_field_flt(in);
    //random_spinor_field_f_flt_cpu(in);

    lprintf("SANITY CHECK", 0, "[Sanity check in field norm unequal zero: %0.2e]\n", spinor_field_sqnorm_f_flt_cpu(in));

    suNf_spinor_flt *in_spinor, *block_start, *out_spinor;
    int stride = 0;
    _PIECE_FOR(in->type, ixp) 
    {
        block_start = _FIELD_BLK(gpu_format, ixp);
        stride = in->type->master_end[ixp] - in->type->master_start[ixp] + 1;
        _SITE_FOR(in->type, ixp, ix) 
        {
            in_spinor = _FIELD_AT(in, ix);
            out_spinor = _FIELD_AT(out, ix);
            int ix_loc = _GPU_IDX_TO_LOCAL(in, ix, ixp);
            for (int comp = 0; comp < 4; ++comp)
            {
                write_gpu_suNf_vector(stride, (*in_spinor).c[comp], block_start, ix_loc, comp);
                read_gpu_suNf_vector(stride, (*out_spinor).c[comp], block_start, ix_loc, comp);
            }
        }
    }

    lprintf("SANITY CHECK", 0, "[Sanity check out field norm unequal zero: %0.2e]\n", spinor_field_sqnorm_f_flt_cpu(out));
    /*spinor_field_sub_assign_f_flt_cpu(out, in);
    double diff_norm = spinor_field_sqnorm_f_flt_cpu(out);
    check_diff_norm_zero(diff_norm);*/

    free_spinor_field_f_flt(in);
    free_spinor_field_f_flt(gpu_format);
    free_spinor_field_f_flt(out);
    return return_val;
}

int test_write_read_sfield() 
{
    lprintf("INFO", 0, " ======= TEST SCALAR FIELD ======= ");
    int return_val = 0;
    scalar_field *in, *gpu_format, *out;

    in = alloc_sfield(1, &glattice);
    gpu_format = alloc_sfield(1, &glattice);
    out = alloc_sfield(1, &glattice);

    random_sfield_cpu(in);
    lprintf("SANITY CHECK", 0, "[Sanity check in field norm unequal zero: %0.2e]\n", sqnorm_sfield_cpu(in));

    double *in_site, *block_start, *out_site;
    int stride = 0;
    _PIECE_FOR(in->type, ixp) 
    {
        block_start = _FIELD_BLK(gpu_format, ixp);
        stride = in->type->master_end[ixp] - in->type->master_start[ixp] + 1;
        _SITE_FOR(in->type, ixp, ix) 
        {
            int ix_loc = _GPU_IDX_TO_LOCAL(in, ix, ixp);
            in_site = _FIELD_AT(in, ix);
            out_site = _FIELD_AT(out, ix);
            write_gpu_double(stride, (*in_site), block_start, ix_loc, 0);
            read_gpu_double(stride, (*out_site), block_start, ix_loc, 0);
        }
    }

    lprintf("SANITY CHECK", 0, "[Sanity check in field norm unequal zero: %0.2e]\n", sqnorm_sfield_cpu(in));
    lprintf("SANITY CHECK", 0, "[Sanity check out field norm unequal zero: %0.2e]\n", sqnorm_sfield_cpu(out));
    sub_assign_sfield_cpu(out, in);
    double diff_norm = sqnorm_sfield_cpu(out);

    check_diff_norm_zero(diff_norm);

    free_sfield(in);
    free_sfield(gpu_format);
    free_sfield(out);
    return return_val;
}

int test_write_read_scalar_field() 
{
    lprintf("INFO", 0, " ======= TEST SU(NG) SCALAR FIELD ======= ");
    int return_val = 0;
    suNg_scalar_field *in, *gpu_format, *out;

    in = alloc_scalar_field(&glattice);
    out = alloc_scalar_field(&glattice);
    gpu_format = alloc_scalar_field(&glattice);

    random_scalar_field_cpu(in);
    lprintf("SANITY CHECK", 0, "[In field norm unequal zero: %0.2e]\n", sqnorm_scalar_field_cpu(in));

    suNg_vector *in_vec, *block_start, *out_vec;
    int stride = 0;
    _PIECE_FOR(in->type, ixp) 
    {
        block_start = _FIELD_BLK(gpu_format, ixp);
        stride = in->type->master_end[ixp] - in->type->master_end[ixp] + 1;
        _SITE_FOR(in->type, ixp, ix) 
        {
            int ix_loc = _GPU_IDX_TO_LOCAL(in, ix, ixp);
            in_vec = _FIELD_AT(in, ix);
            out_vec = _FIELD_AT(out, ix);
            write_gpu_suNg_vector(stride, (*in_vec), block_start, ix_loc, 0);
            read_gpu_suNg_vector(stride, (*out_vec), block_start, ix_loc, 0);
        }
    }

    lprintf("SANITY CHECK", 0, "[Sanity check in field norm unequal zero: %0.2e]\n", sqnorm_scalar_field_cpu(in));
    lprintf("SANITY CHECK", 0, "[Sanity check out field norm unequal zero: %0.2e]\n", sqnorm_scalar_field_cpu(out));
    sub_assign_scalar_field_cpu(out, in);
    double diff_norm = sqnorm_scalar_field_cpu(out);

    check_diff_norm_zero(diff_norm);

    free_scalar_field(in);
    free_scalar_field(gpu_format);
    free_scalar_field(out);
    return return_val;
}
