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

//int test_write_read_spinor_field_f_vector_wise();
int test_write_read_gauge_field_f();
int test_write_read_gauge_field();
//int test_write_read_spinor_field_f_flt_vector_wise();
int test_write_read_spinor_field_f_spinor_wise();

int main(int argc, char *argv[]) 
{
    // Init
    int return_val = 0;

    // Setup process and communication
    logger_map("DEBUG", "debug");
    setup_process(&argc, &argv);

    //return_val += test_write_read_spinor_field_f_vector_wise();
    return_val += test_write_read_spinor_field_f_spinor_wise();
    return_val += test_write_read_gauge_field_f();
    return_val += test_write_read_gauge_field();

    // Single precision
    //return_val += test_write_read_spinor_field_f_flt_vector_wise();

    // Finalize and return
    finalize_process();
    return return_val;
}

int test_write_read_gauge_field()
{
    lprintf("INFO", 0, " ======= TEST GAUGE FIELD ======= ");
    int vol4h = T*X*Y*Z/2;
    int return_val = 0;
    double diff_norm = 0.0;
    double sqnorm = 0.0;
    double sqnorm_in_check = 0.0; // For sanity checks
    double sqnorm_out_check = 0.0; // For sanity checks
    suNg_field *in, *gpu_format, *out;
    in = alloc_gfield(&glattice);
    out = alloc_gfield(&glattice);
    gpu_format = alloc_gfield(&glattice);
    random_u(in);

    suNg *in_mat, *out_mat;
    int dim = sizeof(in->ptr)/sizeof(double);
    
    _MASTER_FOR(in->type, ix) 
    {
        for (int comp = 0; comp < dim; comp++) 
        {
            in_mat = _4FIELD_AT(in, ix, comp);
            out_mat = _4FIELD_AT(out, ix, comp);
            write_gpu_suNg(vol4h, (*in_mat), gpu_format->ptr, ix, comp);
            read_gpu_suNg(vol4h, (*out_mat), gpu_format->ptr, ix, comp);

            _suNg_sqnorm(sqnorm, (*in_mat));
            sqnorm_in_check += sqnorm;

            _suNg_sqnorm(sqnorm, (*out_mat));
            sqnorm_out_check += sqnorm;

            _suNg_sub_assign((*out_mat), (*in_mat));
            _suNg_sqnorm(sqnorm, (*out_mat));
            diff_norm += sqnorm;
        }
    }

    lprintf("SANITY CHECK", 0, "[Sanity check in field norm unequal zero: %0.15lf]\n", sqnorm_in_check);
    lprintf("SANITY CHECK", 0, "[Sanity check out field norm unequal zero: %0.15lf]\n", sqnorm_out_check);

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
    free_gfield(gpu_format);
    free_gfield(out);
    return return_val;
}

int test_write_read_gauge_field_f()
{
    lprintf("INFO", 0, " ======= TEST GAUGE FIELD FUNDAMENTAL REP ======= ");
    int vol4h = T*X*Y*Z/2;
    int return_val = 0;
    double diff_norm = 0.0;
    double sqnorm = 0.0;
    double sqnorm_in_check = 0.0; // For sanity checks
    double sqnorm_out_check = 0.0; // For sanity checks
    suNf_field *in, *gpu_format, *out;
    in = alloc_gfield_f(&glattice);
    out = alloc_gfield_f(&glattice);
    gpu_format = alloc_gfield_f(&glattice);
    random_u_f(in);
    random_u_f(out);

    suNf *in_mat, *out_mat, *tmp_mat;
    int dim = sizeof(in->ptr)/sizeof(double);

    _MASTER_FOR(in->type, ix) 
    {
        for (int comp = 0; comp < dim; comp++) 
        {
            in_mat = _4FIELD_AT(in, ix, comp);
            out_mat = _4FIELD_AT(in, ix, comp);
            write_gpu_suNf(vol4h, (*in_mat), gpu_format->ptr, ix, comp);
            read_gpu_suNf(vol4h, (*out_mat), gpu_format->ptr, ix, comp);

            _suNf_sqnorm(sqnorm, (*in_mat));
            sqnorm_in_check += sqnorm;

            _suNf_sqnorm(sqnorm, (*out_mat));
            sqnorm_out_check += sqnorm;

            _suNf_sub_assign((*out_mat), (*in_mat));
            _suNf_sqnorm(sqnorm, (*out_mat));
            diff_norm += sqnorm;
        }
    }

    lprintf("SANITY CHECK", 0, "[Sanity check in field norm unequal zero: %0.15lf]\n", sqnorm_in_check);
    lprintf("SANITY CHECK", 0, "[Sanity check out field norm unequal zero: %0.15lf]\n", sqnorm_out_check);

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
    free_gfield_f(gpu_format);
    free_gfield_f(out);
    return return_val;
}

int test_write_read_spinor_field_f_spinor_wise() 
{
    int vol4h = T*X*Y*Z/2;
    int return_val = 0;
    spinor_field *in, *gpu_format, *out;
    in = alloc_spinor_field_f(1, &glattice);
    gpu_format = alloc_spinor_field_f(1, &glattice);
    out = alloc_spinor_field_f(1, &glattice);
    gaussian_spinor_field(in);
    lprintf("SANITY CHECK", 0, "[Sanity check]");

    suNf_spinor *in_spinor, *out_spinor;

    _MASTER_FOR(in->type, ix) 
    {
        in_spinor = _FIELD_AT(in, ix);
        out_spinor = _FIELD_AT(out, ix);
        for (int comp = 0; comp < 4; ++comp) 
        {
            write_gpu_suNf_spinor(vol4h, (*in_spinor), gpu_format->ptr, ix, comp);
            read_gpu_suNf_spinor(vol4h, (*out_spinor), gpu_format->ptr, ix, comp);
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

    free_spinor_field_f(in);
    free_spinor_field_f(gpu_format);
    free_spinor_field_f(out);
    return return_val;
}

/*int test_write_read_spinor_field_f_vector_wise()
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

    suNf_spinor *site;

    _PIECE_FOR(in->type, ixp)
    {
        int eo_stride = in->type->master_start[1];
        _SITE_FOR(in->type, ixp, ix)  
        {
            site = _FIELD_AT(in, ix);
            for (int comp=0; comp < 4; comp++) 
            {
                int local_ix = ix % vol4h;
                write_gpu_suNf_vector(vol4h, (*site).c[comp], gpu_format->ptr+ixp*eo_stride, local_ix, comp);
                read_gpu_suNf_vector(vol4h, (*site).c[comp], gpu_format->ptr+ixp*eo_stride, local_ix, comp);
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
} */

/*int test_write_read_spinor_field_f_flt_vector_wise()
{
    lprintf("INFO", 0, " ======= TEST SPINOR FIELD SINGLE PRECISION ======= ");
    int vol4h = T*X*Y*Z/2;
    int return_val = 0;
    spinor_field_flt *in, *gpu_format, *out;
    in = alloc_spinor_field_f_flt(1, &glattice);
    gpu_format = alloc_spinor_field_f_flt(1, &glattice);
    out = alloc_spinor_field_f_flt(1, &glattice);
    gaussian_spinor_field_flt(in);
    lprintf("SANITY CHECK", 0, "[Sanity check in field norm unequal zero: %0.15lf]\n", spinor_field_sqnorm_f_flt_cpu(in));

    suNf_vector_flt in_vec, out_vec;

    _MASTER_FOR(in->type, ix) 
    {
        for (int comp=0; comp < 4; comp++) 
        {
            write_gpu_suNf_vector_flt(vol4h, (*(in->ptr+ix)).c[comp], gpu_format->ptr, ix, comp);
            read_gpu_suNf_vector_flt(vol4h, (*(out->ptr+ix)).c[comp], gpu_format->ptr, ix, comp);
        }
    }

    lprintf("SANITY CHECK", 0, "[Sanity check out field norm unequal zero: %0.15lf]\n", spinor_field_sqnorm_f_flt_cpu(out));
    spinor_field_sub_assign_f_flt_cpu(out, in);
    double diff_norm = spinor_field_sqnorm_f_flt_cpu(out);

    if (diff_norm != 0) {
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
} */