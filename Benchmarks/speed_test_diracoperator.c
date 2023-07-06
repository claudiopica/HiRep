/*******************************************************************************
 *
 * Speed test of Dirac Operator for CPU double and single precision
 *
 *******************************************************************************/

#include "libhr.h"

int main(int argc, char *argv[]) {
    int n_warmup = 100;
    double time_target = 5000.; //number of milliseconds the benchmark will run
    int ninputs = 1;
    int noutputs = 1;
    int n_reps = 0;
    Timer clock;

    setup_process(&argc, &argv);
    setup_random_gauge_fields();

    spinor_field *in = alloc_spinor_field_f(ninputs + noutputs, &glattice);
    setup_random_fields(ninputs + noutputs, in);

    // Check speed diracoperator
    int flopsite = flops_per_site(DPHI_CORE);
    int bytesite = bytes_per_site(DPHI_CORE);

    // speed test Dirac operator
    _WARMUP_SPEEDTEST(clock, n_warmup, time_target, n_reps, Dphi_(&in[0], &in[1]));
    _RUN_SPEEDTEST(clock, n_warmup, time_target, n_reps, flopsite, bytesite, Dphi_(&in[0], &in[1]));

    //measure timer resolution
    evaluate_timer_resolution(clock);

    free_spinor_field_f(in);
    finalize_process();
    return 0;
}
