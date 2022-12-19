#define MAIN_PROGRAM

// Copy sendbuf_gpu_ptr to sendbuf_ptr and check that it is the same
// as syncing on the CPU

#include <stdbool.h>
#include "geometry.h"
#include "memory.h"
#include "update.h"
#include "global.h"
#include "suN_types.h"
#include "dirac.h"
#include "linear_algebra.h"
#include "inverters.h"
#include "logger.h"
#include "setup.h"
#include "hr_complex.h"
#include "random.h"
#include "representation.h"
#include "communications.h"

int run_test();

int main(int argc, char *argv[]) 
{
    int return_val = 0;
    setup_process(&argc, &argv);
    setup_gauge_fields();
    random_u(u_gauge);
    represent_gauge_field();
    copy_to_gpu_gfield_f(u_gauge_f);

    return_val += run_test();

    finalize_process();
    return return_val;
}

int run_test() 
{
    lprintf("INFO", 0, "[Testing GPU Communications]\n");

    spinor_field *s, *q;
    int return_val = 0;

    s = alloc_spinor_field_f(1, &glattice);
    q = alloc_spinor_field_f(1, &glattice);

    gaussian_spinor_field(s);
    spinor_field_mul_f(q, 1.0, s); // Create copy
    copy_to_gpu_spinor_field_f(s);

    start_sendrecv_gpu_spinor_field_f(s);
    complete_sendrecv_gpu_spinor_field_f(s);

    start_sf_sendrecv(q);
    complete_sf_sendrecv(q);

    int total_buffer_length = 0;
    for (int i = 0; i < s->type->nbuffers_spinor; ++i) {
        total_buffer_length += s->type->sbuf_len[i];
    }

    /*for (int i = 0; i < total_buffer_length; ++i) {
        suNf_spinor spinor_s = *(q->sendbuf_ptr + i);
        printf("spinor comp: %0.2 + %0.2e", creal(spinor_s.c[0].c[0]), cimag(spinor_s.c[0].c[0]));
    }*/

    //printf("Total buffer length: %d\n", total_buffer_length);
    sync_gpu_spinor_field_f(s);

    suNf_spinor* host_copy = (suNf_spinor*)malloc(total_buffer_length*sizeof(suNf_spinor));
    cudaMemcpy(host_copy, s->sendbuf_gpu_ptr, total_buffer_length*sizeof(suNf_spinor), cudaMemcpyDeviceToHost);

    for (int i = 0; i < s->type->nbuffers_spinor; ++i) 
    {
        int start = s->type->sbuf_start[i];
        int stride = s->type->sbuf_len[i];
        for (int ix = 0; ix < stride; ++ix) 
        {
            suNf_spinor *out, *in, *cpu_cmp;
            out = (suNf_spinor*)malloc(sizeof(suNf_spinor));
            in = host_copy + start;

            cpu_cmp = q->sendbuf_ptr + ix + start;
            read_gpu_suNf_spinor(stride, 
                                (*out), 
                                in, 
                                ix, 0);
            printf("base index: %d, ix: %d, ix_loc: %d, spinor comp: %0.2e + i%0.2e\n", start, ix+start, ix, creal((*out).c[0].c[0]), cimag((*out).c[0].c[0]));
            printf("CPU cmp: %0.2e + i%0.2e\n", creal((*cpu_cmp).c[0].c[0]), cimag((*cpu_cmp).c[0].c[0]));
        }
    }

    /*int total_buffers_gauge;
    for (int i = 0; i < s->type->nbuffers_gauge; ++i) {
        total_buffers_gauge += s->type->sbuf_len[i];
    }

    setup_gauge_fields();
    random_u(u_gauge);
    represent_gauge_field();
    copy_to_gpu_gfield_f(u_gauge_f);
    sync_gpu_gfield_f(u_gauge_f);
    cudaMemcpy(u_gauge_f->sendbuf_ptr, u_gauge_f->sendbuf_gpu_ptr, total_buffers_gauge*sizeof(suNf_spinor), cudaMemcpyDeviceToHost);

    for (int i = 0; i < s->type->nbuffers_spinor; ++i) 
    {
        int start = s->type->sbuf_start[i];
        int stride = s->type->sbuf_len[i];
        for (int ix = 0; ix < stride; ++ix) 
        {
            suNf *out, *in;
            out = (suNf*)malloc(sizeof(suNf));
            in  = u_gauge_f->sendbuf_ptr + start;
            read_gpu_suNf(stride, (*out), in, ix, 0);
            printf("base index: %d, ix: %d, gauge comp: %0.2e + i%0.2e\n", start, ix, creal((*out).c[0]), cimag((*out).c[0]));
        }
    }*/

    return 0;
}