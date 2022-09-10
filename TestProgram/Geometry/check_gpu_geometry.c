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
//int test_write_read_gauge_field();

int main(int argc, char *argv[]) 
{
    // Init
    int return_val = 0;

    // Setup process and communication
    logger_map("DEBUG", "debug");
    setup_process(&argc, &argv);

    return_val += test_write_read_spinor_field_f();
    return_val += test_write_read_gauge_field_f();
    //return_val += test_write_read_gauge_field();


    // Finalize and return
    finalize_process();
    return return_val;
}

/*int test_write_read_gauge_suNg_av() 
{
    int vol4h = T*X*Y*Z/2;
    int return_val = 0;
    suNg_av_field *in, *gpu_format, *out;
    in = alloc_avfield(&glattice);
    out = alloc_avfield(&glattice);
    gpu_format = alloc_avfield(&glattice);

    suNg_algebra_vector in_vec, out_vec;
    int dim = 4*sizeof(in->ptr)/sizeof(double);

    _PIECE_FOR(in->type, ixp) 
    {
        _SITE_FOR(in->type, ixp, ix) 
        {
            for (int comp = 0; comp < dim; comp++) 
            {
                _suNg_av_read_gpu(vol4h, (*(in->ptr+ix)).c[comp], gpu_format->ptr, ix, comp);
                _suNg_av_read_gpu(vol4h, (*(out->ptr+ix)).c[comp], gpu_format->ptr, ix, comp);
            }
        }
        
        

}*/

/*int test_write_read_gauge_field()
{
    int vol4h = T*X*Y*Z/2;
    int return_val = 0;
    suNg_field *in, *gpu_format, *out;
    in = alloc_gfield(&glattice);
    out = alloc_gfield(&glattice);
    gpu_format = alloc_gfield(&glattice);

    suNg in_mat, out_mat;
    int dim = sizeof(in->ptr)/sizeof(double);
    
    _PIECE_FOR(in->type, ixp) 
    {
        _SITE_FOR(in->type, ixp, ix) 
        {
            for (int comp = 0; comp < dim; comp++) 
            {
                _suNg_write_gpu(vol4h, (*(in->ptr+ix)).c[comp], gpu_format->ptr, ix, comp);
                _suNg_read_gpu(vol4h, (*(out->ptr+ix)).c[comp], gpu_format->ptr, ix, comp);
            }
        }
    }

    suNg_field_sub_assign_cpu(out, in);
    double diff_norm = spinor_field_sqnorm_cpu(out);

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
}*/

int test_write_read_gauge_field_f()
{
    int vol4h = T*X*Y*Z/2;
    int return_val = 0;
    double diff_norm = 0.0;
    double sqnorm = 0.0;
    suNf_field *in, *gpu_format, *out;
    in = alloc_gfield_f(&glattice);
    out = alloc_gfield_f(&glattice);
    gpu_format = alloc_gfield_f(&glattice);

    suNf in_mat, out_mat;
    int dim = sizeof(in->ptr)/sizeof(double);
    
    _PIECE_FOR(in->type, ixp) 
    {
        _SITE_FOR(in->type, ixp, ix) 
        {
            for (int comp = 0; comp < dim; comp++) 
            {
                in_mat = *(in->ptr+ix);
                out_mat = *(out->ptr+ix);
                _suNf_write_gpu(vol4h, in_mat, gpu_format->ptr, ix, comp);
                _suNf_read_gpu(vol4h, out_mat, gpu_format->ptr, ix, comp);


                _suNf_sub_assign(out_mat, in_mat);
                _suNf_sqnorm(sqnorm, out_mat);
                diff_norm += sqnorm;
            }
        }
    }

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
    int vol4h = T*X*Y*Z/2;
    int return_val = 0;
    spinor_field *in, *gpu_format, *out;
    in = alloc_spinor_field_f(1, &glattice);
    gpu_format = alloc_spinor_field_f(1, &glattice);
    out = alloc_spinor_field_f(1, &glattice);
    gaussian_spinor_field(in);

    suNf_vector in_vec, out_vec;

    _PIECE_FOR(in->type, ixp) 
    {
        _SITE_FOR(in->type, ixp, ix) 
        {
            for (int comp=0; comp < 4; comp++) 
            {
                _suNf_write_spinor_gpu(vol4h, (*(in->ptr+ix)).c[comp], gpu_format->ptr, ix, comp);
                _suNf_read_spinor_gpu(vol4h, (*(out->ptr+ix)).c[comp], gpu_format->ptr, ix, comp);
            }
        }
    }

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