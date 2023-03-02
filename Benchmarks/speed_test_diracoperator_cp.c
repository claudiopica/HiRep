/*******************************************************************************
 *
 * NOCOMPILE= WITH_GPU
 *
 * Speed test of Dirac Operator for CPU double and single precision
 *
 *******************************************************************************/

#include "libhr.h"

int main(int argc, char *argv[]) {
    int n_warmup = 100;
    double time_target = 5000.; //number of milliseconds the benchmark will run

    setup_process(&argc, &argv);

    setup_gauge_fields();

    lprintf("MAIN", 0, "n_warmup = %d\n", n_warmup);

    /* allocate memory */
    lprintf("MAIN", 0, "Allocating spinor field\n");
    spinor_field *s0 = alloc_spinor_field_f(2, &glattice);
    spinor_field *s1 = alloc_spinor_field_f(2, &glattice);

    lprintf("MAIN", 0, "Randomizing spinor field...\n");
    gaussian_spinor_field(s0);
    gaussian_spinor_field(s1);

    double res1 = spinor_field_sqnorm_f(s0);
    double res2 = spinor_field_sqnorm_f(s1);

    lprintf("LA_TEST", 0, "Square norm of the random spinor fields %lf and %lf\nThey must be different from zero\n", res1,
            res2);

    lprintf("MAIN", 0, "Generating a random gauge field... ");

    random_u(u_gauge);
    represent_gauge_field();

    lprintf("MAIN", 0, "done.\n");

    // Check speed diracoperator

#if defined(REPR_ADJOINT)
    int flopsite = 8 * NF * (7 + 8 * NF);
#else
    int flopsite = 8 * NF * (7 + 16 * NF);
#endif
    int bytesite = 40 * sizeof(suNf_vector) + 8 * sizeof(suNf); // 8 spinors read + 1 spinor read+write + 8 gauge matrices read

    lprintf("LA TEST", 0, "Flop per site = %d\n", flopsite);
    lprintf("LA TEST", 0, "Byte per site = %d\n", bytesite);
    lprintf("LA TEST", 0, "Dirac data movement = %e MB\n", (double)bytesite / (1024. * 1024.) * GLB_VOLUME);

    // speed test Dirac operator
    lprintf("LA TEST", 0, "Warmup application of the Diracoperator %d times.\n", n_warmup);
    Timer clock;
    timer_set(&clock);
    for (int i = 0; i < n_warmup; ++i) {
        Dphi_(s1, s0);
    }
    double elapsed = timer_lap(&clock) * 1.e-3; //time in milliseconds

    int n_reps = (int)(n_warmup * 1.01 * (time_target / elapsed));
    bcast_int(&n_reps, 1);

    lprintf("LA TEST", 0, "reps: %d , total time: %lf msec, time single: %lf usec\n", n_warmup, elapsed,
            elapsed / n_warmup * 1000.);

    lprintf("LA TEST", 0, "\nEvaluating the massless Diracoperator.\n");
    do {
        lprintf("LA TEST", 0, "Trying reps: %d\n", n_reps);

        elapsed = timer_lap(&clock) * 1.e-3; //time in milliseconds
        for (int i = 0; i < n_reps; ++i) {
            Dphi_(s1, s0);
        }
        elapsed = timer_lap(&clock) * 1.e-3; //time in milliseconds
        n_reps = (int)((double)(n_reps * 1.01 * time_target) / elapsed);
        bcast_int(&n_reps, 1);
    } while (elapsed < time_target * .95);

    lprintf("LA TEST", 0,
            "Massless Diracoperator reps: %d , total time: %lf msec, time single: %lf usec, GFLOPS: %1.6g , BAND: %1.6g GB/s\n",
            n_reps, elapsed, elapsed / n_reps * 1000., (((double)n_reps * GLB_VOLUME) * flopsite) / elapsed / 1.e6,
            (((double)n_reps * GLB_VOLUME) * bytesite) / elapsed / 1.e6);

    //measure timer resolution
    for (int i = 0; i < 1; ++i) {
        elapsed = timer_lap(&clock); //time in microseconds
        elapsed = timer_lap(&clock); //time in microseconds
        lprintf("LA_TEST", 0, "Timer resolution = %lf usec\n", elapsed);
    }
    lprintf("LA_TEST", 0, "Nominal timer resolution = %lf usec\n", timer_res());

    free_spinor_field_f(s0);
    free_spinor_field_f(s1);

    finalize_process();
    return 0;
}
