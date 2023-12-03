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
//#include "random.h"
//#include "memory.h"

#ifdef __cplusplus
extern "C" {
#endif

static double const EPSILON_TEST = 1.e-14;
static double const EPSILON_FLT_TEST = 1.e-4;

#ifdef WITH_GPU
#define synchronize cudaDeviceSynchronize();
#else
#define synchronize
#endif

#ifdef WITH_GPU
#define gpu_copy(_s) copy_to_gpu(_s)
#define cpu_copy(_s) copy_from_gpu(_s)
#define maxnorm(_s) max_cpu(_s)
#define twonorm(_s) sqnorm_cpu(_s)
#else
#define gpu_copy(_s)
#define cpu_copy(_s)
#define maxnorm(_s) max(_s)
#define twonorm(_s) sqnorm(_s)
#endif

#define _FIELD_DESC(_s1)                                     \
    _Generic((_s1),                                          \
        spinor_field *: "spinor field",                      \
        spinor_field_flt *: "single precision spinor field", \
        scalar_field *: "scalar field",                      \
        suNg_field *: "suNg field",                          \
        suNf_field *: "suNf field",                          \
        suNfc_field *: "suNfc_field",                        \
        suNg_field_flt *: "single precision suNg_field",     \
        suNf_field_flt *: "single precision suNf_field",     \
        suNg_scalar_field *: "suNg_scalar_field",            \
        suNg_av_field *: "avfield",                          \
        gtransf *: "gtransf",                                \
        clover_term *: "clover term",                        \
        clover_force *: "clover force",                      \
        staple_field *: "staple field")

#define sanity_check(_n, _s)                                                                                   \
    for (int _k = 0; _k < _n; ++_k) {                                                                          \
        lprintf("SANITY CHECK", 10, "L2 norm of field nr %d: %lf (should be nonzero)\n", _k, sqnorm(&_s[_k])); \
    }

#define compare_cpu_cpu(_errors, _out, _diff, _tag, _precision)                       \
    do {                                                                              \
        sub_assign(_diff, _out);                                                      \
        cpu_copy(_diff);                                                              \
        double res = maxnorm(_diff);                                                  \
        double norm2 = twonorm(_diff);                                                \
        const char *msg = (res > _precision) ? ++_errors, "[FAIL]" : "[ OK ]";        \
        lprintf(_tag, 2, "%s MAX norm=%.10e L2 norm=%.10e\n", msg, res, sqrt(norm2)); \
    } while (0)

#ifdef WITH_GPU
/// @brief Compare two spinor fields in the cpu and gpu parts of out by comparing the MAX and L2 norm
/// @param out Input spinor_field. Th function compare its cpu and gpu parts
/// @param diff Additional spinor_field used for scratch work space
#define compare_cpu_gpu(_errors, _out, _diff, _tag, _precision)           \
    copy(_diff, _out);                                                    \
    copy_from_gpu(_diff);                                                 \
    sub_assign_cpu(_diff, _out);                                          \
    double res = maxnorm(_diff);                                          \
    double norm2 = sqnorm_cpu(_diff);                                     \
    const char *msg = (res > _precision) ? ++errors, "[FAIL]" : "[ OK ]"; \
    lprintf(_tag, 2, "%s MAX norm=%.10e L2 norm=%.10e\n", msg, res, sqrt(norm2));
#endif

#define setup_random_fields_m(_n, _s)                    \
    lprintf("MAIN", 10, "Setup random spinor fields\n"); \
    for (int _i = 0; _i < _n; _i++) {                    \
        random_field(&_s[_i]);                           \
        gpu_copy(&_s[_i]);                               \
        sanity_check(_n, _s);                            \
    }

#define _TEST_CPU_INV_OP(_errors, _name, _ninputs, _in, _out, _test, _tag, _precision) \
    do {                                                                               \
        setup_random_fields_m(_ninputs, _in);                                          \
        _test lprintf(_tag, 2, "%35s: ", _name);                                       \
        compare_cpu_cpu(_errors, _in, (_out), _tag, _precision);                       \
    } while (0)

#define _TEST_RED_INV_OP(_errors, _name, _ninputs, _in, _out, _test, _tag, _precision) \
    do {                                                                               \
        setup_random_fields_m(_ninputs, _in);                                          \
        _test lprintf(_tag, 2, "%35s: ", _name);                                       \
        compare_diff(_errors, abs1, abs2, _tag, _precision);                           \
    } while (0)

#define _TEST_GPU_OP(_errors, _name, _ninputs, _in, _out, _test, _tag, _prec) \
    do {                                                                      \
        setup_random_fields_m(_ninputs, _in);                                 \
        _test lprintf(_tag, 2, "%15s: ", _name);                              \
        compare_cpu_gpu(_errors, _out, (_out + 1), _tag, _prec);              \
    } while (0)

#define _TEST_RED_OP(_errors, _name, _ninputs, _in, _test, _tag, _prec) \
    do {                                                                \
        setup_random_fields_m(_ninputs, _in);                           \
        _test lprintf(_tag, 2, "%15s: ", _name);                        \
        compare_diff(_errors, abs1, abs2, _tag, _prec);                 \
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
void compare_diff(int errors, double abs1, double abs2, char tag[], double prec);
void compare_diff_flt(int errors, float abs1, float abs2);

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

// RANDOM
void random_spinor_field_cpu(spinor_field *);
void random_spinor_field_flt_cpu(spinor_field_flt *);
void random_suNg_field_cpu(suNg_field *);
void random_suNf_field_cpu(suNf_field *);
void random_suNfc_field_cpu(suNfc_field *);
void random_suNg_field_flt_cpu(suNg_field_flt *);
void random_suNf_field_flt_cpu(suNf_field_flt *);
void random_suNg_scalar_field_cpu(suNg_scalar_field *);
void random_suNg_av_field_cpu(suNg_av_field *);
void random_scalar_field_cpu(scalar_field *);
void random_gtransf_cpu(gtransf *);
void random_ldl_field_cpu(ldl_field *);
void random_clover_term_cpu(clover_term *);
void random_clover_force_cpu(clover_force *);
void random_staple_field_cpu(staple_field *);

#ifdef __cplusplus
}
#endif
#endif
