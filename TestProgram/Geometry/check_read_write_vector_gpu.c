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

int test_write_read_spinor_field_f();
int test_write_read_gauge_field_f();
int test_write_read_gauge_field();
int test_write_read_spinor_field_f_flt();

int main(int argc, char *argv[]) 
{
    // Init
    int return_val = 0;

    // Setup process and communication
    logger_map("DEBUG", "debug");
    setup_process(&argc, &argv);

    return_val += test_write_read_spinor_field_f();
    return_val += test_write_read_gauge_field_f();
    return_val += test_write_read_gauge_field();

    // Single precision
    return_val += test_write_read_spinor_field_f_flt();

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

    suNg in_mat, out_mat;
    int dim = sizeof(in->ptr)/sizeof(double);
    
    _MASTER_FOR(in->type, ix) 
    {
        for (int comp = 0; comp < dim; comp++) 
        {
            in_mat = *(in->ptr+ix);
            out_mat = *(out->ptr+ix);
            write_gpu_suNg(vol4h, in_mat, gpu_format->ptr, ix, comp);
            read_gpu_suNg(vol4h, out_mat, gpu_format->ptr, ix, comp);

            _suNg_sqnorm(sqnorm, in_mat);
            sqnorm_in_check += sqnorm;

            _suNg_sqnorm(sqnorm, out_mat);
            sqnorm_out_check += sqnorm;

            _suNg_sub_assign(out_mat, in_mat);
            _suNg_sqnorm(sqnorm, out_mat);
            diff_norm += sqnorm;
        }
    }

    lprintf("INFO", 0, "[Sanity check in field norm unequal zero: %0.15lf]\n", sqnorm_in_check);
    lprintf("INFO", 0, "[Sanity check out field norm unequal zero: %0.15lf]\n", sqnorm_out_check);

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

    suNf in_mat, out_mat;
    int dim = sizeof(in->ptr)/sizeof(double);
    
    _MASTER_FOR(in->type, ix) 
    {
        for (int comp = 0; comp < dim; comp++) 
        {
            in_mat = *(in->ptr+ix);
            out_mat = *(out->ptr+ix);
            write_gpu_suNf(vol4h, in_mat, gpu_format->ptr, ix, comp);
            read_gpu_suNf(vol4h, out_mat, gpu_format->ptr, ix, comp);


            _suNf_sqnorm(sqnorm, in_mat);
            sqnorm_in_check += sqnorm;

            _suNf_sqnorm(sqnorm, out_mat);
            sqnorm_out_check += sqnorm;

            _suNf_sub_assign(out_mat, in_mat);
            _suNf_sqnorm(sqnorm, out_mat);
            diff_norm += sqnorm;
        }
    }

    lprintf("INFO", 0, "[Sanity check in field norm unequal zero: %0.15lf]\n", sqnorm_in_check);
    lprintf("INFO", 0, "[Sanity check out field norm unequal zero: %0.15lf]\n", sqnorm_out_check);

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
    lprintf("INFO", 0, "[Sanity check in field norm unequal zero: %0.15lf]\n", spinor_field_sqnorm_f_cpu(in));

    suNf_vector in_vec, out_vec;

    _MASTER_FOR(in->type, ix) 
    {
        for (int comp=0; comp < 4; comp++) 
        {
            write_gpu_suNf_vector(vol4h, (*(in->ptr+ix)).c[comp], gpu_format->ptr, ix, comp);
            read_gpu_suNf_vector(vol4h, (*(out->ptr+ix)).c[comp], gpu_format->ptr, ix, comp);
        }
    }

    lprintf("INFO", 0, "[Sanity check out field norm unequal zero: %0.15lf]\n", spinor_field_sqnorm_f_cpu(out));
    spinor_field_sub_assign_f_cpu(out, in);
    double diff_norm = spinor_field_sqnorm_f_cpu(out);

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
    lprintf("INFO", 0, "[Sanity check in field norm unequal zero: %0.15lf]\n", spinor_field_sqnorm_f_flt_cpu(in));

    suNf_vector_flt in_vec, out_vec;

    _MASTER_FOR(in->type, ix) 
    {
        for (int comp=0; comp < 4; comp++) 
        {
            write_gpu_suNf_vector_flt(vol4h, (*(in->ptr+ix)).c[comp], gpu_format->ptr, ix, comp);
            read_gpu_suNf_vector_flt(vol4h, (*(out->ptr+ix)).c[comp], gpu_format->ptr, ix, comp);
        }
    }

    lprintf("INFO", 0, "[Sanity check out field norm unequal zero: %0.15lf]\n", spinor_field_sqnorm_f_flt_cpu(out));
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
} 