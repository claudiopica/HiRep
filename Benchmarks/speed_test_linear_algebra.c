/*******************************************************************************
*
* Speed test of double precision linear algebra functions
*
*******************************************************************************/

#include "libhr.h"

static int n_warmup = 50;
#ifdef WITH_GPU
static double time_target = 50.;
#else
static double time_target = 500.;
#endif

#ifdef WITH_GPU
#define synchronize cudaDeviceSynchronize();
#else
#define synchronize
#endif

void setup_random_fields(int n, spinor_field s[n]) {
    for (int i = 0; i < n; i++) {
        gaussian_spinor_field(&s[i]);
    }
}

int bytes_per_site(int ninputs, int noutputs, int sitesize) {
    return (ninputs + 2 * noutputs) * sitesize;
}

#define _PRINT_SETUP(_ninputs, _noutputs, _flopsite, _in)                                                  \
    do {                                                                                                   \
        int bytesite = bytes_per_site(_ninputs, _noutputs, sizeof((&_in[0])->ptr));                        \
        lprintf("LA TEST", 0, "Flop per size = %d\n", _flopsite);                                          \
        lprintf("LA TEST", 0, "Byte per site = %d\n", bytesite);                                           \
        lprintf("LA TEST", 0, "Data movement = %e MB\n", (double)bytesite / (1024. * 1024.) * GLB_VOLUME); \
    } while (0);

#define _WARMUP(_name, _elapsed, _n_reps, _clock, _in, _test, _ninputs)            \
    lprintf("WARMUP", 0, "Warmup application of %s %d times.\n", _name, n_warmup); \
    _elapsed = 0;                                                                  \
    for (int i = 0; i < n_warmup; ++i) {                                           \
        setup_random_fields(_ninputs, _in);                                        \
        timer_lap(&_clock);                                                        \
        _test;                                                                     \
        synchronize;                                                               \
        _elapsed += timer_lap(&_clock) * 1.e-3;                                    \
    }                                                                              \
    int _n_reps = (int)(n_warmup * 1.01 * (time_target / _elapsed));               \
    bcast_int(&n_reps, 1);                                                         \
    lprintf("WARMUP", 0,                                                           \
            "reps: %d, "                                                           \
            "total time: %lf msec, "                                               \
            "time single: %lf usec\n",                                             \
            n_warmup, _elapsed, _elapsed / n_warmup * 1000.);

#define _SPEEDTEST(_name, _elapsed, _n_reps, _clock, _in, _test, _ninputs)  \
    lprintf("LA TEST", 0, "Evaluating %s\n", _name);                        \
    do {                                                                    \
        lprintf("LA TEST", 0, "Trying reps: %d\n", _n_reps);                \
        _elapsed = 0;                                                       \
        for (int i = 0; i < _n_reps; ++i) {                                 \
            setup_random_fields(_ninputs, _in);                             \
            timer_lap(&_clock);                                             \
            _test;                                                          \
            synchronize;                                                    \
            _elapsed += timer_lap(&_clock) * 1.e-3;                         \
        }                                                                   \
        _n_reps = (int)((double)(_n_reps * 1.01 * time_target) / _elapsed); \
        bcast_int(&n_reps, 1);                                              \
    } while (_elapsed < time_target * .95);

#ifdef WITH_GPU
#define _PERFORMANCE_REACHED(_bandwidth)                                                             \
    do {                                                                                             \
        int n_devices;                                                                               \
        cudaGetDeviceCount(&n_devices);                                                              \
        double peak_memory_bandwidth = 0;                                                            \
        for (int i = 0; i < n_devices; i++) {                                                        \
            cudaDeviceProp prop;                                                                     \
            cudaGetDeviceProperties(&prop, i);                                                       \
            peak_memory_bandwidth += 2.0 * prop.memoryClockRate * (prop.memoryBusWidth / 8) / 1.0e6; \
        }                                                                                            \
        lprintf("RESULT", 0,                                                                         \
                "Peak Memory Bandwidth: %1.6g GB/s, "                                                \
                "Performance reached: %0.2f \%\n",                                                   \
                peak_memory_bandwidth, _bandwidth / peak_memory_bandwidth);                          \
    } while (0);

#else
#define _PERFORMANCE_REACHED(_bandwidth)
#endif

#define _PRINT_RESULT(_name, _ninputs, _noutputs, _elapsed, _n_reps, _flopsite, _in)          \
    do {                                                                                      \
        int bytesite = bytes_per_site(_ninputs, _noutputs, sizeof((*(&_in[0])->ptr)));        \
        double vol_reps = (double)_n_reps * GLB_VOLUME;                                       \
        double bandwidth = (vol_reps * bytesite) / _elapsed / 1.e6;                           \
        lprintf("RESULT", 0,                                                                  \
                "%s reps: %d, "                                                               \
                "total time: %lf msec, "                                                      \
                "time single: %lf usec, "                                                     \
                "GFLOPS: %1.6g\n"                                                             \
                "BANDWIDTH: %1.6g GB/s, ",                                                    \
                _name, _n_reps, _elapsed, _elapsed / _n_reps * 1000., _flopsite * GLB_VOLUME, \
                (vol_reps * _flopsite) / _elapsed / 1.e6, bandwidth);                         \
        _PERFORMANCE_REACHED(bandwidth);                                                      \
    } while (0);

#define _SPEED_TEST_LIN_ALG(_name, _ninputs, _noutputs, _flopsite, _in, _test)      \
    do {                                                                            \
        setup_random_fields(_ninputs, _in);                                         \
        _PRINT_SETUP(_ninputs, _noutputs, _flopsite, _in);                          \
        Timer clock;                                                                \
        double elapsed;                                                             \
        timer_set(&clock);                                                          \
        _WARMUP(_name, elapsed, n_reps, clock, _in, _test, _ninputs);               \
        _SPEEDTEST(_name, elapsed, n_reps, clock, _in, _test, _ninputs);            \
        _PRINT_RESULT(_name, _ninputs, _noutputs, elapsed, n_reps, _flopsite, _in); \
    } while (0);

int main(int argc, char *argv[]) {
    setup_process(&argc, &argv);

    int ninputs = 3;
    spinor_field *in;
    in = alloc_spinor_field(ninputs, &glattice);

    hr_complex c;
    _complex_i_add(c, 1.5, 2.5);

    Timer clock;
    timer_set(&clock);

    _SPEED_TEST_LIN_ALG("Copy spinor field", 1, 1, 0, in, spinor_field_copy_f(&in[0], &in[1]););

    _SPEED_TEST_LIN_ALG("Spinor field product", 1, 1, NF * 4 * 2, in, spinor_field_prod_f(&in[0], &in[1]););

    _SPEED_TEST_LIN_ALG("Square norm", 1, 0, NF * 4 * 2, in, spinor_field_sqnorm_f(&in[0]););

    _SPEED_TEST_LIN_ALG("g5 application", 1, 1, NF * 2 * 2, in, spinor_field_g5_f(&in[0], &in[1]););

    _SPEED_TEST_LIN_ALG("g5 mulc add assign", 1, 1, NF * 6 * 2, in, spinor_field_g5_mulc_add_assign_f(&in[0], c, &in[1]););

    // Timer resolution
    for (int i = 0; i < 1; ++i) {
        timer_lap(&clock);
        timer_lap(&clock);
        lprintf("LA TEST", 0, "Nominal timer resolution = %lf usec\n", timer_res());
    }

    free_spinor_field(in);
    finalize_process();
    return 0;
}