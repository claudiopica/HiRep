/*******************************************************************************
 *
 * Speed test of reduction operations (incl. plaquette)
 *
 * NOCOMPILE= !WITH_GPU
 *
 *******************************************************************************/

#include "libhr.h"

int main(int argc, char *argv[]) {
    int n_warmup = 100;
    double time_target = 5000.;
    int ninputs = 1;
    int n_reps = 0;
    Timer clock;

    setup_process(&argc, &argv);
    setup_random_gauge_fields();

    spinor_field *in = alloc_spinor_field(1, &glattice);
    gaussian_spinor_field(in);

    lprintf("REDUCTION BENCHMARKS", 1, "Testing naive CUDA reduction of doubles\n");
    int flopsite = flops_per_site(CUDA_REDUCTION) * sizeof(*u_gauge->ptr) / sizeof(double);
    int bytesite = bytes_per_site(CUDA_REDUCTION) * sizeof(*u_gauge->ptr) / sizeof(double);

    // The global sum is the most efficient reduction we can do
    // TODO; only do reduction, bypass the host to device copy
    _WARMUP_SPEEDTEST(clock, n_warmup, time_target, n_reps, global_sum_gpu_double((double *)u_gauge->ptr, 1));
    _RUN_SPEEDTEST(clock, n_warmup, time_target, n_reps, flopsite, bytesite, global_sum_gpu_double((double *)u_gauge->ptr, 1));

    // How about the spinor field sqnorm?
    lprintf("REDUCTION BENCHMARKS", 1, "Testing spinor field sqnorm\n");
    flopsite = flops_per_site(SF_SQNORM);
    bytesite = bytes_per_site(SF_SQNORM);
    _WARMUP_SPEEDTEST(clock, n_warmup, time_target, n_reps, sqnorm(in));
    _RUN_SPEEDTEST(clock, n_warmup, time_target, n_reps, flopsite, bytesite, sqnorm(in));

    flopsite = flops_per_site(PLAQUETTE);
    bytesite = bytes_per_site(PLAQUETTE);
    lprintf("REDUCTION BENCHMARKS", 1, "Testing average plaquette\n");
    _WARMUP_SPEEDTEST(clock, n_warmup, time_target, n_reps, avr_plaquette());
    _RUN_SPEEDTEST(clock, n_warmup, time_target, n_reps, flopsite, bytesite, avr_plaquette());

    return 0;
}
