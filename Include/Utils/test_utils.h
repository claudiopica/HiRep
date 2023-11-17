/***************************************************************************\
* Copyright (c) 2022, Sofie Martins, Claudio Pica                           *   
* All rights reserved.                                                      * 
\***************************************************************************/

/**
 * @file test_utils.h
 * @brief utilities to simplify testing
 */

#ifndef TEST_UTILS_H
#define TEST_UTILS_H

#include "libhr_core.h"
#include "timing.h"

#ifdef __cplusplus
extern "C" {
#endif

#ifdef WITH_GPU
#define synchronize cudaDeviceSynchronize();
#else
#define synchronize
#endif

#define _TEST_GPU_OP(_errors, _name, _ninputs, _in, _out, _test) \
    do {                                                         \
        setup_random_fields(_ninputs, _in);                      \
        spinor_field *out = _out;                                \
        spinor_field *diff = out + 1;                            \
        _test lprintf("GPU TEST", 2, "%15s: ", _name);           \
        compare_cpu_gpu(errors, out, diff);                      \
    } while (0)

#define _TEST_RED_OP(_errors, _name, _ninputs, _in, _test) \
    do {                                                   \
        setup_random_fields(_ninputs, _in);                \
        _test lprintf("GPU TEST", 2, "%15s: ", _name);     \
        compare_diff(_errors, abs1, abs2);                 \
    } while (0)

#define _TEST_GPU_OP_FLT(_errors, _name, _ninputs, _in, _out, _test) \
    do {                                                             \
        setup_random_fields_flt(_ninputs, _in);                      \
        spinor_field_flt *out = _out;                                \
        spinor_field_flt *diff = out + 1;                            \
        _test lprintf("GPU TEST", 2, "%15s: ", _name);               \
        compare_cpu_gpu_flt(_errors, out, diff);                     \
    } while (0)

#define _TEST_RED_OP_FLT(_errors, _name, _ninputs, _in, _test) \
    do {                                                       \
        setup_random_fields_flt(_ninputs, _in);                \
        _test lprintf("GPU TEST", 2, "%15s: ", _name);         \
        compare_diff_flt(_errors, abs1, abs2);                 \
    } while (0)

#define _WARMUP_SPEEDTEST(_clock, _n_warmup, _time_target, _n_reps, _operator)        \
    lprintf("LA TEST", 0, "Warmup application %d times.\n", _n_warmup);               \
    do {                                                                              \
        timer_set(&clock);                                                            \
        for (int i = 0; i < _n_warmup; ++i) {                                         \
            _operator;                                                                \
        }                                                                             \
        synchronize;                                                                  \
        double elapsed = timer_lap(&clock) * 1.e-3;                                   \
        lprintf("LA TEST", 0, "total time: %lf msec\n", elapsed);                     \
        lprintf("LA TEST", 0, "time single: %lf usec\n", elapsed / n_warmup * 1000.); \
        n_reps = (int)(n_warmup * 1.01 * (_time_target / elapsed));                   \
        bcast_int(&n_reps, 1);                                                        \
    } while (0)

#define _RUN_SPEEDTEST(_clock, _n_warmup, _time_target, _n_reps, _flopsite, _bytesite, _operator)                            \
    do {                                                                                                                     \
        double __elapsed = 0;                                                                                                \
        do {                                                                                                                 \
            timer_lap(&_clock);                                                                                              \
            for (int i = 0; i < _n_reps; ++i) {                                                                              \
                _operator;                                                                                                   \
            }                                                                                                                \
            synchronize;                                                                                                     \
            __elapsed = timer_lap(&clock) * 1.e-3;                                                                           \
            _n_reps = (int)((double)(n_reps * 1.01 * time_target) / __elapsed);                                              \
            bcast_int(&n_reps, 1);                                                                                           \
        } while (__elapsed < _time_target * .95);                                                                            \
                                                                                                                             \
        lprintf("LA TEST", 0, "Number of repetitions: %d\n", _n_reps);                                                       \
        lprintf("LA TEST", 0, "Total time: %lf msec\n", __elapsed);                                                          \
        lprintf("LA TEST", 0, "Time single: %lf usec\n", __elapsed / _n_reps * 1000.);                                       \
        lprintf("LA TEST", 0, "GFLOPS: %1.6g\n", (((double)_n_reps * GLB_VOLUME) * _flopsite) / __elapsed / 1.e6);           \
        lprintf("LA TEST", 0, "BANDWIDTH: %1.6g GB/s\n\n", (((double)_n_reps * GLB_VOLUME) * _bytesite) / __elapsed / 1.e6); \
    } while (0)

double spinor_max(suNf_spinor *s);
float spinor_max_flt(suNf_spinor_flt *s);
double spinor_field_findmax_f(spinor_field *in);
float spinor_field_findmax_f_flt(spinor_field_flt *in);
void compare_diff(int errors, double abs1, double abs2);
void compare_diff_flt(int errors, float abs1, float abs2);
void compare_cpu_gpu(int errors, spinor_field *out, spinor_field *diff);
void compare_cpu_gpu_flt(int errors, spinor_field_flt *out, spinor_field_flt *diff);

void evaluate_timer_resolution(Timer clock);
void setup_random_fields(int n, spinor_field s[]);
void setup_random_fields_flt(int n, spinor_field_flt s[]);
void spinor_field_sanity_check(int ninputs, spinor_field *in);
void spinor_field_sanity_check_flt(int ninputs, spinor_field_flt *in);

void setup_random_gauge_fields();
void setup_clover();

// Testing utils
void test_setup();
int check_diff_norm(double, double);
int check_diff_norm_zero(double);
int check_finiteness(double);

// COPY
void copy_suNg_field_cpu(suNg_field *, suNg_field *);
void copy_suNg_scalar_field_cpu(suNg_scalar_field *, suNg_scalar_field *);
void copy_suNg_field_flt_cpu(suNg_field_flt *, suNg_field_flt *);
void copy_suNf_field_cpu(suNf_field *, suNf_field *);
void copy_suNf_field_flt_cpu(suNf_field_flt *, suNf_field_flt *);
void copy_suNg_av_field_cpu(suNg_av_field *, suNg_av_field *);
void copy_scalar_field_cpu(scalar_field *, scalar_field *);
void copy_ldl_field_cpu(ldl_field *, ldl_field *);
void copy_gtransf_cpu(gtransf *, gtransf *);
void copy_clover_term_cpu(clover_term *, clover_term *);
void copy_clover_force_cpu(clover_force *, clover_force *);
void copy_staple_field_cpu(staple_field *, staple_field *);

// RANDOM
void random_spinor_field_cpu(spinor_field *);
void random_spinor_field_flt_cpu(spinor_field_flt *);
void random_suNg_field_cpu(suNg_field *);
void random_suNg_scalar_field_cpu(suNg_scalar_field *);
void random_suNg_field_flt_cpu(suNg_field_flt *);
void random_suNf_field_cpu(suNf_field *);
void random_suNf_field_flt_cpu(suNf_field_flt *);
void random_suNg_av_field_cpu(suNg_av_field *);
void random_scalar_field_cpu(scalar_field *);
void random_gtransf_cpu(gtransf *);
void random_ldl_field_cpu(ldl_field *);
void random_clover_term_cpu(clover_term *);
void random_clover_force_cpu(clover_force *);
void random_staple_field_cpu(staple_field *);

// SUB ASSIGN
void sub_assign_suNg_field_cpu(suNg_field *, suNg_field *);
void sub_assign_suNg_field_flt_cpu(suNg_field_flt *, suNg_field_flt *);
void sub_assign_suNf_field_cpu(suNf_field *, suNf_field *);
void sub_assign_suNf_field_flt_cpu(suNf_field_flt *, suNf_field_flt *);
void sub_assign_suNg_scalar_field_cpu(suNg_scalar_field *, suNg_scalar_field *);
void sub_assign_suNg_av_field_cpu(suNg_av_field *, suNg_av_field *);
void sub_assign_scalar_field_cpu(scalar_field *, scalar_field *);
void sub_assign_gtransf_cpu(gtransf *, gtransf *);
void sub_assign_ldl_field_cpu(ldl_field *, ldl_field *);
void sub_assign_clover_term_cpu(clover_term *, clover_term *);
void sub_assign_clover_force_cpu(clover_force *, clover_force *);
void sub_assign_staple_field_cpu(staple_field *, staple_field *);

// SQNORM
double sqnorm_spinor_field_cpu(spinor_field *);
float sqnorm_spinor_field_flt_cpu(spinor_field_flt *);
double sqnorm_suNg_field_cpu(suNg_field *);
double sqnorm_suNf_field_cpu(suNf_field *);
float sqnorm_suNg_field_flt_cpu(suNg_field_flt *);
float sqnorm_suNf_field_flt_cpu(suNf_field_flt *);
double sqnorm_suNg_scalar_field_cpu(suNg_scalar_field *);
double sqnorm_suNg_av_field_cpu(suNg_av_field *);
double sqnorm_gtransf_cpu(gtransf *);
double sqnorm_scalar_field_cpu(scalar_field *);
double sqnorm_ldl_field_cpu(ldl_field *);
double sqnorm_clover_term_cpu(clover_term *);
double sqnorm_clover_force_cpu(clover_force *);
double sqnorm_staple_field_cpu(staple_field *);

// SET ZERO
void zero_suNg_field_cpu(suNg_field *);
void zero_scalar_field_cpu(scalar_field *);
void zero_suNf_field_cpu(suNf_field *);
void zero_suNg_field_flt_cpu(suNg_field_flt *);
void zero_suNf_field_flt_cpu(suNf_field_flt *);
void zero_suNg_scalar_field_cpu(suNg_scalar_field *);
void zero_suNg_av_field_cpu(suNg_av_field *);
void zero_gtransf_cpu(gtransf *);
void zero_ldl_field_cpu(ldl_field *);
void zero_clover_term_cpu(clover_term *);
void zero_clover_force_cpu(clover_force *);
void zero_staple_field_cpu(staple_field *);

#ifdef __cplusplus
}
#endif
#endif
