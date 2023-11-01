/*******************************************************************************
 *
 * NOCOMPILE= !WITH_CLOVER && !WITH_EXPCLOVER
 * NOCOMPILE= !DPHI_FLT
 * 
 * Speed test of Dirac Operator for CPU double and single precision
 *
 *******************************************************************************/

#include "libhr.h"

static double hmass = 0.1;

#define __cphi Cphi_flt(-hmass, &in[0], &in[1])
#define __cphi_inv_ Cphi_diag_inv_flt(-hmass, &in[0], &in[1])
#define __cphi_ Cphi_diag_flt(-hmass, &in[0], &in[1])

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

    spinor_field_flt *in = alloc_spinor_field_flt(ninputs + noutputs, &glattice);
    setup_random_fields_flt(ninputs + noutputs, in);

    u_gauge_f_flt = alloc_suNf_field_flt(&glattice);

    lprintf("LA TEST", 0, "Speedtesting application of single precision clover-improved dirac operator\n");
    int flopsite = flops_per_site(CPHI_FLT);
    int bytesite = bytes_per_site(CPHI_FLT);

    _WARMUP_SPEEDTEST(clock, n_warmup, time_target, n_reps, __cphi);
    _RUN_SPEEDTEST(clock, n_warmup, time_target, n_reps, flopsite, bytesite, __cphi);

    lprintf("LA TEST", 0, "Speedtesting application of single precision clover term\n");
    flopsite = flops_per_site(CPHI_FLT_CORE);
    bytesite = bytes_per_site(CPHI_FLT_CORE);

    _WARMUP_SPEEDTEST(clock, n_warmup, time_target, n_reps, __cphi_);
    _RUN_SPEEDTEST(clock, n_warmup, time_target, n_reps, flopsite, bytesite, __cphi_);

    lprintf("LA TEST", 0, "Speedtesting application of single precision clover term inverse\n");
    flopsite = flops_per_site(CPHI_INV_FLT);
    bytesite = bytes_per_site(CPHI_INV_FLT);

    _WARMUP_SPEEDTEST(clock, n_warmup, time_target, n_reps, __cphi_inv_);
    _RUN_SPEEDTEST(clock, n_warmup, time_target, n_reps, flopsite, bytesite, __cphi_inv_);

    free_spinor_field_flt(in);
    free_suNf_field_flt(u_gauge_f_flt);
    free_suNg_field_flt(u_gauge_flt);
    finalize_process();
    return 0;
}
