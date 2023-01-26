/*******************************************************************************
*
* NOCOMPILE= !WITH_GPU
*
* Speed test of Dirac Operator for GPU
*
*******************************************************************************/

#include "libhr.h"

int main(int argc, char* argv[]) {
    setup_process(&argc, &argv);
    setup_gauge_fields();
    int n_warmup = 20;
    int number_of_operations = 200;
    u_gauge_f_flt = alloc_gfield_f_flt(&glattice);
    struct timeval start, end, etime;

    spinor_field* in = alloc_spinor_field_f(1, &glattice);
    spinor_field* out = alloc_spinor_field_f(1, &glattice);

    for (int i = 0; i < n_warmup; ++i) {
        to_gpu_format_spinor_field_f(out, in);
        to_cpu_format_spinor_field_f(in, out);
    }

    gettimeofday(&start, 0);
    for (int i = 0; i < number_of_operations; ++i) {
        to_gpu_format_spinor_field_f(out, in);
        to_cpu_format_spinor_field_f(in, out);
    }
    gettimeofday(&end, 0);
    timeval_subtract(&etime, &end, &start);
    lprintf("RESULT",0,"Time: [%ld sec %ld usec]\n",etime.tv_sec,etime.tv_usec);

    finalize_process();
    return 0;
}