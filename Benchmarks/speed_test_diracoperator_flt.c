/*******************************************************************************
 *
 * NOCOMPILE= !DPHI_FLT
 *
 * Speed test of dirac operator for single precision
 *
 *******************************************************************************/

#include "libhr.h"

int main(int argc, char *argv[]) {
    int n_warmup = 1000;
    double time_target = 5000.;
    int ninputs = 1;
    int noutputs = 1;
    int n_reps = 0;
    Timer clock;

    setup_process(&argc, &argv);
    setup_gauge_fields();

    spinor_field_flt *in = alloc_spinor_field_flt(ninputs + noutputs, &glattice);
    setup_random_fields_flt(ninputs + noutputs, in);

    u_gauge_f_flt = alloc_suNf_field_flt(&glattice);

    lprintf("LA TEST", 0, "Speedtesting application of single precision massless dirac operator (hopping term)\n");
    int flopsite = flops_per_site(DPHI_CORE_FLT);
    int bytesite = bytes_per_site(DPHI_CORE_FLT);

    _WARMUP_SPEEDTEST(clock, n_warmup, time_target, n_reps, Dphi_flt_(&in[0], &in[1]));
    _RUN_SPEEDTEST(clock, n_warmup, time_target, n_reps, flopsite, bytesite, Dphi_flt_(&in[0], &in[1]));

    free_spinor_field_flt(in);
    free_suNf_field_flt(u_gauge_f_flt);
    free_suNg_field_flt(u_gauge_flt);
    finalize_process();
    exit(0);
}
