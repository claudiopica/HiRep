/*******************************************************************************
 *
 * Test random number generation
 *
 *******************************************************************************/

#include "libhr.h"
#include <math.h>

static int n_warmup = 50;
static double time_target = 500.;
int max_iterations = 10;

#define min(a, b) (a > b) ? b : a

void setup_random_fields(int n, spinor_field s[n]) {
    for (int i = 0; i < n; i++) {
        gaussian_spinor_field(&s[i]);
    }
}

#ifdef WITH_GPU
#define synchronize cudaDeviceSynchronize();
#else
#define synchronize
#endif

#define _SPEED_TEST_RNG(_random, _random_init, _type)                              \
    do {                                                                           \
        Timer clock;                                                               \
        timer_set(&clock);                                                         \
        int n_batch = 1e5;                                                         \
        _type r[n_batch];                                                          \
        _random_init(1, 42);                                                       \
        timer_lap(&clock);                                                         \
        for (int i = 0; i < n_warmup; ++i) {                                       \
            _random(r, n_batch);                                                   \
            synchronize;                                                           \
        }                                                                          \
        double elapsed = timer_lap(&clock) * 1.e-3;                                \
                                                                                   \
        int n_reps = (int)(n_warmup * 1.01 * (time_target / elapsed));             \
        global_sum_int(&n_reps, 1);                                                \
                                                                                   \
        lprintf("WARMUP", 0,                                                       \
                "Needed %lf msec to generate %d random numbers. "                  \
                "Time single: %lf usec\n",                                         \
                n_warmup *n_batch, elapsed, elapsed / n_warmup / n_batch * 1000.); \
                                                                                   \
        int i = 0;                                                                 \
        do {                                                                       \
            lprintf("LA TEST", 0, "Trying reps: %d\n", n_reps);                    \
            timer_lap(&clock);                                                     \
            for (int i = 0; i < n_reps; ++i) {                                     \
                _random(r, n_batch);                                               \
                synchronize;                                                       \
            }                                                                      \
            elapsed = timer_lap(&clock) * 1.e-3;                                   \
            n_reps = (int)((double)(n_reps * 1.01 * time_target) / elapsed);       \
            global_sum_int(&n_reps, 1);                                            \
            i++;                                                                   \
        } while (elapsed < time_target * .95 && i < max_iterations);               \
        lprintf("RESULT", 0,                                                       \
                "RNG speed: %d samples generated in %lf, "                         \
                "corresponding to a rate of %e GSamples/s\n",                      \
                n_reps *n_batch, elapsed, n_reps *n_batch / elapsed * 1e-6);       \
    } while (0);

#define _SPEED_TEST_RANDOM_FIELD(_name, _in, _test)                                                          \
    do {                                                                                                     \
        int bytesite = 2 * sizeof(*(_in->ptr));                                                              \
        lprintf("LA TEST", 0, "Byte per site = %d\n", bytesite);                                             \
        lprintf("LA TEST", 0, "Generated data =  %e MB\n", (double)bytesite / (1024. * 1024.) * GLB_VOLUME); \
        Timer clock;                                                                                         \
        timer_set(&clock);                                                                                   \
        lprintf("WARMUP", 0, "Warmup application of %s %d times\n", _name, n_warmup);                        \
        double elapsed = timer_lap(&clock) * 1.e-3;                                                          \
        for (int i = 0; i < n_warmup; ++i) {                                                                 \
            _test                                                                                            \
        };                                                                                                   \
        elapsed = timer_lap(&clock) * 1.e-3;                                                                 \
        int n_reps = (int)(n_warmup * 1.01 * time_target / elapsed);                                         \
        bcast_int(&n_reps, 1);                                                                               \
        lprintf("WARMUP", 0,                                                                                 \
                "reps: %d, "                                                                                 \
                "total time: %lf msec, "                                                                     \
                "time single: %lf usec\n",                                                                   \
                n_warmup, elapsed, elapsed / n_warmup * 1000.);                                              \
        lprintf("LA TEST", 0, "Generating %s\n", _name);                                                     \
        do {                                                                                                 \
            lprintf("LA TEST", 0, "Trying reps: %d\n", n_reps);                                              \
            elapsed = timer_lap(&clock) * 1.e-3;                                                             \
            for (int i = 0; i < n_reps; ++i) {                                                               \
                _test                                                                                        \
            };                                                                                               \
            elapsed = timer_lap(&clock) * 1.e-3;                                                             \
            n_reps = (int)((double)(n_reps * 1.01 * time_target) / elapsed);                                 \
            bcast_int(&n_reps, 1);                                                                           \
        } while (elapsed < time_target * .95);                                                               \
                                                                                                             \
        double vol_reps = (double)n_reps * GLB_VOLUME;                                                       \
        lprintf("RESULT", 0,                                                                                 \
                "%s reps: %d, "                                                                              \
                "total time: %lf msec, "                                                                     \
                "time single: %lf usec, "                                                                    \
                "BANDWIDTH: %1.6g GB/s\n",                                                                   \
                _name, n_reps, elapsed, elapsed / n_reps * 1000., (vol_reps * bytesite) / elapsed / 1.e6);   \
    } while (0);

int main(int argc, char *argv[]) {
    setup_process(&argc, &argv);
    lprintf("MAIN", 0, "n_warmup = %d\n", n_warmup);

    _SPEED_TEST_RNG(ranlxd, rlxd_init, double);
    _SPEED_TEST_RNG(ranlxs, rlxs_init, float);

    spinor_field *in = alloc_spinor_field(1, &glattice);
    _SPEED_TEST_RANDOM_FIELD("Gaussian spinor field", in, gaussian_spinor_field(in););

    spinor_field_flt *in_flt = alloc_spinor_field_flt(1, &glattice);
    _SPEED_TEST_RANDOM_FIELD("Gaussian spinor field single precision", in_flt, gaussian_spinor_field_flt(in_flt););

    suNg_field *u = alloc_suNg_field(&glattice);
    _SPEED_TEST_RANDOM_FIELD("Random gauge field", in, random_u(u););

    suNg_av_field *avf = alloc_suNg_av_field(&glattice);
    _SPEED_TEST_RANDOM_FIELD("Gaussian momenta", in, gaussian_momenta(avf););

    suNg_scalar_field *sf = alloc_suNg_scalar_field(&glattice);
    _SPEED_TEST_RANDOM_FIELD("Gaussian scalar momenta", in, gaussian_scalar_momenta(sf););

    finalize_process();
    return 0;
}