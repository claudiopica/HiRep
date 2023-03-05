/*******************************************************************************
*
* Speed test of strided reading
*
*******************************************************************************/

#include "libhr.h"
#include "./speed_test_geometry_gpu.hpp"

static int n_warmup = 50;
static double time_target = 50.;

static suNf_spinor spinor;
static suNf_spinor *target;

void simple_read_cpu(spinor_field *in) {
    _MASTER_FOR(in->type, ix) {
#ifdef FIXED_STRIDE
        read_gpu_suNf_spinor(0, spinor, in->ptr, ix, 0, 1);
#else
        spinor = *_FIELD_AT(in, ix);
#endif
    }
}

void simple_write_cpu(spinor_field *in) {
    _MASTER_FOR(in->type, ix) {
#ifdef FIXED_STRIDE
        write_gpu_suNf_spinor(0, spinor, in->ptr, ix, 0, 1);
#else
        target = _FIELD_AT(in, ix);
        *target = spinor;
#endif
    }
}

#ifndef WITH_GPU
#define simple_read simple_read_cpu
#define simple_write simple_write_cpu
#define synchronize
#define _PERFORMANCE_REACHED(_bandwidth)
#else
#define simple_read simple_read_gpu
#define simple_write simple_write_gpu
#define synchronize cudaDeviceSynchronize()
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
#endif

int main(int argc, char *argv[]) {
    int bytesite = sizeof(suNf_vector);
    setup_process(&argc, &argv);
    spinor_field *in = alloc_spinor_field_f(1, &glattice);
    gaussian_spinor_field(in);

    Timer clock;
    double elapsed;

    // Read

    /* WARMUP */
    lprintf("WARMUP", 0, "Warmup read %d times\n", n_warmup);
    timer_set(&clock);
    elapsed = 0;
    for (int i = 0; i < n_warmup; ++i) {
        gaussian_spinor_field(in);
        timer_lap(&clock);
        simple_read(in);
        synchronize;
        elapsed += timer_lap(&clock) * 1.e-3;
    }

    int n_reps = (int)(n_warmup * 1.01 * (time_target / elapsed));
    bcast_int(&n_reps, 1);

    lprintf("WARMUP", 0,
            "reps: %d, total time %lf msec, "
            "time single %lf usec\n",
            n_warmup, elapsed, elapsed / n_warmup * 1000.);

    /* LATENCY TEST */
    lprintf("LA TEST", 0, "Reading spinor field.\n");
    do {
        lprintf("LA TEST", 0, "Trying reps: %d\n", n_reps);

        elapsed = 0;
        for (int i = 0; i < n_reps; ++i) {
            gaussian_spinor_field(in);
            timer_lap(&clock);
            simple_read(in);
            synchronize;
            elapsed += timer_lap(&clock) * 1.e-3;
        }

        n_reps = (int)((double)(n_reps * 1.01 * time_target) / elapsed);
        bcast_int(&n_reps, 1);
    } while (elapsed < time_target * .95);

    lprintf("LA TEST", 0,
            "Reading operation: %d, total time: %lf msec, "
            "time single: %lf usec, BANDWIDTH: %1.6g GB/s\n",
            n_reps, elapsed, elapsed / n_reps * 1000., (((double)n_reps * GLB_VOLUME) * bytesite) / elapsed / 1.e6);
    _PERFORMANCE_REACHED(((((double)n_reps * GLB_VOLUME) * bytesite) / elapsed / 1.e6))

    // Write
    /* WARMUP */
    lprintf("WARMUP", 0, "Warmup write %d times\n", n_warmup);
    timer_set(&clock);
    elapsed = 0;
    for (int i = 0; i < n_warmup; ++i) {
        gaussian_spinor_field(in);
        timer_lap(&clock);
        simple_write(in);
        synchronize;
        elapsed += timer_lap(&clock) * 1.e-3;
    }

    n_reps = (int)(n_warmup * 1.01 * (time_target / elapsed));
    bcast_int(&n_reps, 1);

    lprintf("WARMUP", 0,
            "reps: %d, total time %lf msec, "
            "time single %lf usec\n",
            n_warmup, elapsed, elapsed / n_warmup * 1000.);

    /* LATENCY TEST */
    lprintf("LA TEST", 0, "Writing spinor field.\n");
    do {
        lprintf("LA TEST", 0, "Trying reps: %d\n", n_reps);

        elapsed = 0;
        for (int i = 0; i < n_reps; ++i) {
            gaussian_spinor_field(in);
            timer_lap(&clock);
            simple_write(in);
            synchronize;
            elapsed += timer_lap(&clock) * 1.e-3;
        }

        n_reps = (int)((double)(n_reps * 1.01 * time_target) / elapsed);
        bcast_int(&n_reps, 1);
    } while (elapsed < time_target * .95);

    lprintf("LA TEST", 0,
            "Writing operation: %d, total time: %lf msec, "
            "time single: %lf usec, BAND: %1.6g GB/s\n",
            n_reps, elapsed, elapsed / n_reps * 1000., (((double)n_reps * GLB_VOLUME) * 2 * bytesite) / elapsed / 1.e6);
    _PERFORMANCE_REACHED(((((double)n_reps * GLB_VOLUME) * bytesite) / elapsed / 1.e6))

    return 0;
}