/***************************************************************************\
* Copyright (c) 2022, Sofie Martins, Claudio Pica                           *   
* All rights reserved.                                                      * 
\***************************************************************************/

#include "libhr.h"
#include <string.h>

//static double EPSILON = 1.e-14;
//static float EPSILON_FLT = 1.e-4;

// logger level 10 for a more verbose setting

double spinor_max(suNf_spinor *s) {
    double *a = (double *)s;
    double max = 0.;
    for (int i = 0; i < sizeof(suNf_spinor) / sizeof(*a); i++) {
        double v = fabs(a[i]);
        if (max < v) { max = v; }
    }
    return max;
}

float spinor_max_flt(suNf_spinor_flt *s) {
    float *a = (float *)s;
    float max = 0.;
    for (int i = 0; i < sizeof(suNf_spinor_flt) / sizeof(*a); i++) {
        float v = fabs(a[i]);
        if (max < v) { max = v; }
    }
    return max;
}

double spinor_field_findmax_f(spinor_field *in) {
    double max = 0.;
    _ONE_SPINOR_FOR(in) {
        suNf_spinor c = *_SPINOR_PTR(in);
        double v = spinor_max(&c);
        if (max < v) { max = v; }
    }
#ifdef WITH_MPI
    global_max(&max, 1);
#endif
    return max;
}

float spinor_field_findmax_f_flt(spinor_field_flt *in) {
    float max = 0.;
    _ONE_SPINOR_FOR(in) {
        suNf_spinor_flt c = *_SPINOR_PTR(in);
        float v = spinor_max_flt(&c);
        if (max < v) { max = v; }
    }
#ifdef WITH_MPI
    global_max_flt(&max, 1);
#endif
    return max;
}

/// @brief  Check if the two inputs are the same within a given relative precision of EPSILON
/// @param abs1
/// @param abs2
void compare_diff(int errors, double abs1, double abs2, char tag[], double precision) {
    double rel = fabs(abs1 - abs2) / fabs(abs1);
    const char *msg = (rel > precision || !isfinite(rel) || !isfinite(abs1) || !isfinite(abs2)) ? ++errors, "[FAIL]" : "[ OK ]";
    lprintf(tag, 2, "%s rel=%.10e abs=%.10e diff=%.10e\n", msg, rel, fabs(abs1), fabs(abs1 - abs2));
}

void evaluate_timer_resolution(Timer clock) {
    timer_lap(&clock); //time in microseconds
    double elapsed = timer_lap(&clock); //time in microseconds
    lprintf("LA TEST", 0, "Timer resolution = %lf usec\n", elapsed);
    lprintf("LA TEST", 0, "Nominal timer resolution = %lf usec\n", timer_res());
}

void setup_random_gauge_fields() {
    lprintf("MAIN", 10, "Setup random gauge fields\n");
    setup_gauge_fields();
    random_u(u_gauge);

#ifdef DPHI_FLT
    u_gauge_f_flt = alloc_suNf_field_flt(&glattice);
    u_gauge_flt = alloc_suNg_field_flt(&glattice);

    assign_ud2u();
#ifdef WITH_GPU
    assign_ud2u_cpu();
#endif

#endif

    represent_gauge_field();
#ifdef DPHI_FLT
    assign_ud2u_f();
#ifdef WITH_GPU
    assign_ud2u_f_cpu();
#endif
#endif

#if defined(WITH_GPU) && defined(ALLOCATE_REPR_GAUGE_FIELD)
    copy_from_gpu_suNf_field(u_gauge_f);
#endif

#ifdef WITH_MPI
    start_sendrecv_suNf_field(u_gauge_f);
    complete_sendrecv_suNf_field(u_gauge_f);

#endif

#if defined(WITH_GPU) && defined(ALLOCATE_REPR_GAUGE_FIELD) && defined(DPHI_FLT)
    copy_from_gpu_suNf_field_flt(u_gauge_f_flt);
    start_sendrecv_suNf_field_flt(u_gauge_f_flt);
    complete_sendrecv_suNf_field_flt(u_gauge_f_flt);
#endif

    lprintf("MAIN", 10, "done.\n");
}

void setup_clover() {
#if defined(WITH_CLOVER) || defined(WITH_EXPCLOVER)
    lprintf("MAIN", 10, "Setup clover\n");
#ifdef WITH_GPU
    copy_from_gpu_clover_term(cl_term);
#endif
    start_sendrecv_clover_term(cl_term);
    complete_sendrecv_clover_term(cl_term);
#ifdef WITH_GPU
    copy_from_gpu_clover_force(cl_force);
#endif
    lprintf("MAIN", 10, "done.\n");
#endif
}

// different function to setup gauge fields and clover
void setup_random_fields(int n, spinor_field s[]) {
    lprintf("MAIN", 10, "Setup random spinor fields\n");
    for (int i = 0; i < n; i++) {
        gaussian_spinor_field(&s[i]);
#ifdef WITH_GPU
        copy_to_gpu_spinor_field(&s[i]);
#endif
    }
    spinor_field_sanity_check(n, s);
    lprintf("MAIN", 10, "done.\n");
}

/// Generates an array of gaussian spinor fields
/// and copy the results in the gpu memory
void setup_random_fields_flt(int n, spinor_field_flt s[]) {
    lprintf("MAIN", 10, "Setup single precision random spinor fields\n");
    for (int i = 0; i < n; i++) {
        gaussian_spinor_field_flt(&s[i]);
#ifdef WITH_GPU
        copy_to_gpu_spinor_field_flt(&s[i]);
#endif
    }
    spinor_field_sanity_check_flt(n, s);
    lprintf("MAIN", 10, "done.\n");
}

void spinor_field_sanity_check(int ninputs, spinor_field *in) {
    for (int i = 0; i < ninputs; i++) {
        spinor_field *s = in + i;
        double res = sqnorm_spinor_field(s);
        lprintf("SANITY CHECK", 10, "Input spinor field nr %d square norm: %lf (should be nonzero)\n", i, res);
    }
}

void spinor_field_sanity_check_flt(int ninputs, spinor_field_flt *in) {
    for (int i = 0; i < ninputs; i++) {
        spinor_field_flt *s = in + i;
        double res = sqnorm_spinor_field_flt(s);
        lprintf("SANITY CHECK", 10, "Input spinor field nr %d square norm: %lf (should be nonzero)\n", i, res);
    }
}

// The following functions are primarily for testing purposes
// This is all for CPU

void test_setup() {
    rlxd_init(1, 205);
    rlxs_init(2, 208);
}

int check_diff_norm(double diff_norm, double tol) {
    int return_val = 0;
    if (fabs(diff_norm) > tol) {
        lprintf("RESULT", 0, "FAILED \n");
        return_val = 1;
    } else {
        lprintf("RESULT", 0, "OK \n");
        return_val = 0;
    }
    lprintf("RESULT", 0, "[Diff norm %0.2e]\n", diff_norm);
    return return_val;
}

int check_diff_norm_zero(double diff_norm) {
    int return_val = 0;
    if (fabs(diff_norm) != 0) {
        lprintf("RESULT", 0, "FAILED \n");
        return_val = 1;
    } else {
        lprintf("RESULT", 0, "OK \n");
        return_val = 0;
    }
    lprintf("RESULT", 0, "[Diff norm is %0.2e]\n", diff_norm);
    return return_val;
}

int check_not_zero(double value, double tolerance) {
    int return_val = 0;
    if (fabs(value) < tolerance) {
        lprintf("RESULT", 0, "FAILED \n");
        return_val = 1;
    } else {
        lprintf("RESULT", 0, "OK \n");
        return_val = 0;
    }
    lprintf("RESULT", 0, "[Value %0.2e is not zero up to tolerance.]\n", value);
    return return_val;
}

int check_finiteness(double value) {
    int return_val = 0;
    if (!isfinite(value)) {
        lprintf("RESULT", 0, "FAILED \n");
        return_val = 1;
    } else {
        lprintf("RESULT", 0, "OK \n");
        return_val = 0;
    }
    lprintf("RESULT", 0, "[Value is %0.2e.]\n", value);
    return return_val;
}

void rand_field_dbl(double *d, int n) {
    for (int i = 0; i < n; ++i) {
        d[i] = (double)rand() / (double)RAND_MAX;
    }
}

void rand_field_flt(float *f, int n) {
    for (int i = 0; i < n; i++) {
        f[i] = (float)rand() / (float)RAND_MAX;
    }
}

// ** RANDOM FIELDS FOR TESTING **
void random_spinor_field_cpu(spinor_field *f) {
    int n = f->type->gsize_spinor * sizeof(suNf_spinor) / sizeof(double);
    ranlxd((double *)(f->ptr), n);
}

void random_spinor_field_flt_cpu(spinor_field_flt *f) {
    int n = f->type->gsize_spinor * sizeof(suNf_spinor_flt) / sizeof(float);
    ranlxs((float *)(f->ptr), n);
}

void random_suNg_field_cpu(suNg_field *f) {
    int n = 4 * f->type->gsize_gauge * sizeof(suNg) / sizeof(double);
    ranlxd((double *)(f->ptr), n);
}

void random_suNg_scalar_field_cpu(suNg_scalar_field *f) {
    int n = f->type->gsize_gauge * sizeof(suNg_vector) / sizeof(double);
    ranlxd((double *)(f->ptr), n);
}

void random_suNg_field_flt_cpu(suNg_field_flt *f) {
    int n = 4 * f->type->gsize_gauge * sizeof(suNg_flt) / sizeof(float);
    ranlxs((float *)(f->ptr), n);
}

void random_suNf_field_flt_cpu(suNf_field_flt *f) {
    int n = 4 * f->type->gsize_gauge * sizeof(suNf_flt) / sizeof(float);
    ranlxs((float *)(f->ptr), n);
}

void random_suNf_field_cpu(suNf_field *f) {
    int n = 4 * f->type->gsize_gauge * sizeof(suNf) / sizeof(double);
    ranlxd((double *)(f->ptr), n);
}

void random_suNfc_field_cpu(suNfc_field *f) {
    int n = 4 * f->type->gsize_gauge * sizeof(suNfc) / sizeof(double);
    ranlxd((double *)(f->ptr), n);
}

void random_gfield_f_flt_cpu(suNf_field_flt *f) {
    int n = 4 * f->type->gsize_gauge * sizeof(suNf_flt) / sizeof(float);
    ranlxs((float *)(f->ptr), n);
}

void random_suNg_av_field_cpu(suNg_av_field *f) {
    int n = 4 * f->type->gsize_gauge * sizeof(suNg_algebra_vector) / sizeof(double);
    ranlxd((double *)(f->ptr), n);
}

void random_scalar_field_cpu(scalar_field *f) {
    int n = f->type->gsize_spinor;
    ranlxd((double *)(f->ptr), n);
}

void random_gtransf_cpu(gtransf *f) {
    int n = f->type->gsize_gauge * sizeof(suNg) / sizeof(double);
    ranlxd((double *)(f->ptr), n);
}

void random_ldl_field_cpu(ldl_field *f) {
    int n = f->type->gsize_gauge * sizeof(ldl_t) / sizeof(double);
    ranlxd((double *)(f->ptr), n);
}

void random_clover_term_cpu(clover_term *f) {
    int n = 4 * f->type->gsize_gauge * sizeof(suNfc) / sizeof(double);
    ranlxd((double *)(f->ptr), n);
}

void random_clover_force_cpu(clover_force *f) {
    int n = 6 * f->type->gsize_gauge * sizeof(suNf) / sizeof(double);
    ranlxd((double *)(f->ptr), n);
}

void random_staple_field_cpu(staple_field *f) {
    int n = 3 * f->type->gsize_gauge * sizeof(suNg) / sizeof(double);
    ranlxd((double *)(f->ptr), n);
}

#ifdef WITH_GPU

#endif
