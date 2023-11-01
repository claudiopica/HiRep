/*******************************************************************************
 *
 * NOCOMPILE= !WITH_CLOVER && !WITH_EXPCLOVER
 * 
 * Speed test of Dirac Operator for CPU double and single precision
 *
 *******************************************************************************/

#include "libhr.h"

static double hmass = 0.1;

#define __cphi Cphi(-hmass, &in[0], &in[1])
#define __cphi_inv_ Cphi_diag_inv(-hmass, &in[0], &in[1])
#define __cphi_ Cphi_diag(-hmass, &in[0], &in[1])

int main(int argc, char *argv[]) {
    int n_warmup = 100;
    double time_target = 5000.; //number of milliseconds the benchmark will run
    int ninputs = 1;
    int noutputs = 1;
    int n_reps = 0;
    Timer clock;

    setup_process(&argc, &argv);
    setup_random_gauge_fields();
    setup_clover();

    spinor_field *in = alloc_spinor_field(ninputs + noutputs, &glattice);
    setup_random_fields(ninputs + noutputs, in);

    lprintf("LA TEST", 0, "Speedtesting application of clover-improved dirac operator\n");
    int flopsite = flops_per_site(CPHI);
    int bytesite = bytes_per_site(CPHI);

    _WARMUP_SPEEDTEST(clock, n_warmup, time_target, n_reps, __cphi);
    _RUN_SPEEDTEST(clock, n_warmup, time_target, n_reps, flopsite, bytesite, __cphi);

    lprintf("LA TEST", 0, "Speedtesting application of clover term\n");
    flopsite = flops_per_site(CPHI_CORE);
    bytesite = bytes_per_site(CPHI_CORE);

    _WARMUP_SPEEDTEST(clock, n_warmup, time_target, n_reps, __cphi_);
    _RUN_SPEEDTEST(clock, n_warmup, time_target, n_reps, flopsite, bytesite, __cphi_);

    lprintf("LA TEST", 0, "Speedtesting application of clover term inverse\n");
    flopsite = flops_per_site(CPHI_INV);
    bytesite = bytes_per_site(CPHI_INV);

    _WARMUP_SPEEDTEST(clock, n_warmup, time_target, n_reps, __cphi_inv_);
    _RUN_SPEEDTEST(clock, n_warmup, time_target, n_reps, flopsite, bytesite, __cphi_inv_);

    //measure timer resolution
    evaluate_timer_resolution(clock);

    free_spinor_field(in);
    finalize_process();
    return 0;
}
