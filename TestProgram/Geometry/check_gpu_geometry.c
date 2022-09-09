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

int main(int argc, char *argv[]) 
{
    // Init
    int pass = 1;

    // Setup process and communication
    logger_map("DEBUG", "debug");
    setup_process(&argc, &argv);

    pass &= test_write_read_spinor_field_f();
    
    // Finalize and return
    finalize_process();
    return pass;
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
                _suNf_write_spinor_gpu(vol4h/2, (*(in->ptr+ix)).c[comp], gpu_format->ptr, ix, comp);
                _suNf_read_spinor_gpu(vol4h/2, (*(out->ptr+ix)).c[comp], gpu_format->ptr, ix, comp);
            }
        }
    }

    spinor_field_sub_assign_f_cpu(out, in);
    double diff_norm = spinor_field_sqnorm_f_cpu(out);

    if (diff_norm > 1e-14) {
        lprintf("RESULT", 0, "FAILED \n");
        return_val = 1;
    } 
    else 
    {
        lprintf("RESULT", 0, "OK \n");
        return_val = 0;
    }
    lprintf("RESULT", 0, "[Diff norm gpu-cpu %0.2e]\n", diff_norm);

    free_spinor_field_f(in);
    free_spinor_field_f(out);
    return return_val;
} 