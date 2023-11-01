/***************************************************************************\
* Copyright (c) 2022, Sofie Martins, Claudio Pica                           *   
* All rights reserved.                                                      * 
\***************************************************************************/

#include "libhr.h"
#include <string.h>

static double EPSILON = 1.e-14;
static float EPSILON_FLT = 1.e-4;

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
void compare_diff(int errors, double abs1, double abs2) {
    double rel = fabs(abs1 - abs2) / fabs(abs1);
    const char *msg = (rel > EPSILON) ? ++errors, "[FAIL]" : "[ OK ]";
    lprintf("GPU TEST", 2, "%s rel=%.10e abs=%.10e diff=%.10e\n", msg, rel, fabs(abs1), fabs(abs1 - abs2));
}

/// @brief  Check if the two inputs are the same within a given relative precision of EPSILON
/// @param abs1
/// @param abs2
void compare_diff_flt(int errors, float abs1, float abs2) {
    float rel = fabs(abs1 - abs2) / fabs(abs1);
    const char *msg = (rel > EPSILON_FLT) ? ++errors, "[FAIL]" : "[ OK ]";
    lprintf("GPU TEST", 2, "%s rel=%.10e abs=%.10e diff=%.10e\n", msg, rel, fabs(abs1), fabs(abs1 - abs2));
}

#ifdef WITH_GPU
/// @brief Compare two spinor fields in the cpu and gpu parts of out by comparing the MAX and L2 norm
/// @param out Input spinor_field. Th function compare its cpu and gpu parts
/// @param diff Additional spinor_field used for scratch work space
void compare_cpu_gpu(int errors, spinor_field *out, spinor_field *diff) {
    spinor_field_copy_f_gpu(diff, out);
    copy_from_gpu_spinor_field(diff);
    spinor_field_sub_assign_f_cpu(diff, out);
    double res = spinor_field_findmax_f(diff);
    double norm2 = spinor_field_sqnorm_f_cpu(diff);
    const char *msg = (res > EPSILON) ? ++errors, "[FAIL]" : "[ OK ]";
    lprintf("GPU TEST", 2, "%s MAX norm=%.10e L2 norm=%.10e\n", msg, res, sqrt(norm2));
}

/// @brief Compare two spinor fields in the cpu and gpu parts of out by comparing the MAX and L2 norm
/// @param out Input spinor_field. Th function compare its cpu and gpu parts
/// @param diff Additional spinor_field used for scratch work space
void compare_cpu_gpu_flt(int errors, spinor_field_flt *out, spinor_field_flt *diff) {
    spinor_field_copy_f_flt_gpu(diff, out);
    copy_from_gpu_spinor_field_flt(diff);
    spinor_field_sub_assign_f_flt_cpu(diff, out);
    float res = spinor_field_findmax_f_flt(diff);
    float norm2 = spinor_field_sqnorm_f_flt_cpu(diff);
    const char *msg = (res > EPSILON_FLT) ? ++errors, "[FAIL]" : "[ OK ]";
    lprintf("GPU TEST", 2, "%s MAX norm=%.10e L2 norm=%.10e\n", msg, res, sqrt(norm2));
}

#endif

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
        double res = spinor_field_sqnorm_f(s);
        lprintf("SANITY CHECK", 10, "Input spinor field nr %d square norm: %lf (should be nonzero)\n", i, res);
    }
}

void spinor_field_sanity_check_flt(int ninputs, spinor_field_flt *in) {
    for (int i = 0; i < ninputs; i++) {
        spinor_field_flt *s = in + i;
        double res = spinor_field_sqnorm_f_flt(s);
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

void copy_suNg_field_cpu(suNg_field *out, suNg_field *in) {
    memcpy(out->ptr, in->ptr, 4 * out->type->gsize_gauge * sizeof(suNg));
}

void copy_suNg_scalar_field_cpu(suNg_scalar_field *out, suNg_scalar_field *in) {
    memcpy(out->ptr, in->ptr, out->type->gsize_gauge * sizeof(suNg_vector));
}

void copy_suNg_field_flt_cpu(suNg_field_flt *out, suNg_field_flt *in) {
    memcpy(out->ptr, in->ptr, 4 * out->type->gsize_gauge * sizeof(suNg_flt));
}

void copy_suNf_field_cpu(suNf_field *out, suNf_field *in) {
    memcpy(out->ptr, in->ptr, 4 * out->type->gsize_gauge * sizeof(suNf));
}

void copy_suNf_field_flt_cpu(suNf_field_flt *out, suNf_field_flt *in) {
    memcpy(out->ptr, in->ptr, 4 * out->type->gsize_gauge * sizeof(suNf_flt));
}

void copy_suNg_av_field_cpu(suNg_av_field *out, suNg_av_field *in) {
    memcpy(out->ptr, in->ptr, 4 * out->type->gsize_gauge * sizeof(suNg_algebra_vector));
}

void copy_scalar_field_cpu(scalar_field *out, scalar_field *in) {
    memcpy(out->ptr, in->ptr, out->type->gsize_spinor * sizeof(double));
}

void copy_gtransf_cpu(gtransf *out, gtransf *in) {
    memcpy(out->ptr, in->ptr, out->type->gsize_spinor * sizeof(suNg));
}

void copy_ldl_field_cpu(ldl_field *out, ldl_field *in) {
    memcpy(out->ptr, in->ptr, out->type->gsize_gauge * sizeof(ldl_t));
}

void copy_clover_term_cpu(clover_term *out, clover_term *in) {
    memcpy(out->ptr, in->ptr, 4 * out->type->gsize_gauge * sizeof(suNfc));
}

void copy_clover_force_cpu(clover_force *out, clover_force *in) {
    memcpy(out->ptr, in->ptr, 6 * out->type->gsize_gauge * sizeof(suNf));
}

void copy_staple_field_cpu(staple_field *out, staple_field *in) {
    memcpy(out->ptr, in->ptr, 3 * out->type->gsize_gauge * sizeof(suNg));
}

// ** SUB ASSIGN **
void sub_assign_suNg_field_cpu(suNg_field *out, suNg_field *in) {
    suNg *site_out, *site_in;
    _MASTER_FOR(in->type, ix) {
        for (int comp = 0; comp < 4; comp++) {
            site_out = _4FIELD_AT(out, ix, comp);
            site_in = _4FIELD_AT(in, ix, comp);
            _suNg_sub_assign((*site_out), (*site_in));
        }
    }
}

void sub_assign_suNg_field_flt_cpu(suNg_field_flt *out, suNg_field_flt *in) {
    suNg_flt *site_out, *site_in;
    _MASTER_FOR(in->type, ix) {
        for (int comp = 0; comp < 4; comp++) {
            site_out = _4FIELD_AT(out, ix, comp);
            site_in = _4FIELD_AT(in, ix, comp);
            _suNg_sub_assign((*site_out), (*site_in));
        }
    }
}

void sub_assign_suNf_field_cpu(suNf_field *out, suNf_field *in) {
    suNf *site_out, *site_in;
    _MASTER_FOR(in->type, ix) {
        for (int comp = 0; comp < 4; comp++) {
            site_out = _4FIELD_AT(out, ix, comp);
            site_in = _4FIELD_AT(in, ix, comp);
            _suNf_sub_assign((*site_out), (*site_in));
        }
    }
}

void sub_assign_suNf_field_flt_cpu(suNf_field_flt *out, suNf_field_flt *in) {
    suNf_flt *site_out, *site_in;
    _MASTER_FOR(in->type, ix) {
        for (int comp = 0; comp < 4; comp++) {
            site_out = _4FIELD_AT(out, ix, comp);
            site_in = _4FIELD_AT(in, ix, comp);
            _suNf_sub_assign((*site_out), (*site_in));
        }
    }
}

void sub_assign_suNg_scalar_field_cpu(suNg_scalar_field *out, suNg_scalar_field *in) {
    suNg_vector *site_out, *site_in;
    _MASTER_FOR(in->type, ix) {
        site_out = _FIELD_AT(out, ix);
        site_in = _FIELD_AT(in, ix);
        _vector_sub_assign_g((*site_out), (*site_in));
    }
}

void sub_assign_suNg_av_field_cpu(suNg_av_field *out, suNg_av_field *in) {
    suNg_algebra_vector *site_in, *site_out;
    _MASTER_FOR(in->type, ix) {
        for (int comp = 0; comp < 4; comp++) {
            site_out = _4FIELD_AT(out, ix, comp);
            site_in = _4FIELD_AT(in, ix, comp);
            _algebra_vector_sub_assign_g((*site_out), (*site_in));
        }
    }
}

void sub_assign_scalar_field_cpu(scalar_field *out, scalar_field *in) {
    double *site_in, *site_out;
    _MASTER_FOR(in->type, ix) {
        site_out = _FIELD_AT(out, ix);
        site_in = _FIELD_AT(in, ix);
        (*site_out) -= (*site_in);
    }
}

void sub_assign_gtransf_cpu(gtransf *out, gtransf *in) {
    suNg *site_in, *site_out;
    _MASTER_FOR(in->type, ix) {
        site_out = _FIELD_AT(out, ix);
        site_in = _FIELD_AT(in, ix);
        _suNg_sub_assign((*site_out), (*site_in));
    }
}

void sub_assign_ldl_field_cpu(ldl_field *out, ldl_field *in) {
    ldl_t *site_out, *site_in;
    _MASTER_FOR(in->type, ix) {
        site_out = _FIELD_AT(out, ix);
        site_in = _FIELD_AT(in, ix);
        for (int i = 0; i < NF * (2 * NF + 1); ++i) {
            _complex_sub_assign((*site_out).up[i], (*site_in).up[i]);
            _complex_sub_assign((*site_out).dn[i], (*site_in).dn[i]);
        }
    }
}

void sub_assign_clover_term_cpu(clover_term *out, clover_term *in) {
    suNfc *site_out, *site_in;
    _MASTER_FOR(in->type, ix) {
        for (int comp = 0; comp < 4; comp++) {
            site_out = _4FIELD_AT(out, ix, comp);
            site_in = _4FIELD_AT(in, ix, comp);
            _suNf_sub_assign((*site_out), (*site_in));
        }
    }
}

void sub_assign_clover_force_cpu(clover_force *out, clover_force *in) {
    suNf *site_out, *site_in;
    _MASTER_FOR(in->type, ix) {
        for (int comp = 0; comp < 6; comp++) {
            site_out = _6FIELD_AT(out, ix, comp);
            site_in = _6FIELD_AT(in, ix, comp);
            _suNf_sub_assign((*site_out), (*site_in));
        }
    }
}

void sub_assign_staple_field_cpu(staple_field *out, staple_field *in) {
    suNg *site_out, *site_in;
    _MASTER_FOR(in->type, ix) {
        for (int comp = 0; comp < 3; comp++) {
            site_out = _3FIELD_AT(out, ix, comp);
            site_in = _3FIELD_AT(in, ix, comp);
            _suNg_sub_assign((*site_out), (*site_in));
        }
    }
}

// ** SQNORM **
double sqnorm_suNg_field_cpu(suNg_field *f) {
    suNg *site;
    double sqnorm = 0.0;
    _MASTER_FOR(f->type, ix) {
        for (int comp = 0; comp < 4; comp++) {
            double tmp;
            site = _4FIELD_AT(f, ix, comp);
            _suNg_sqnorm(tmp, (*site));
            sqnorm += tmp;
        }
    }

#ifdef WITH_MPI
    global_sum(&sqnorm, 1);
#endif
    return sqnorm;
}

float sqnorm_suNg_field_flt_cpu(suNg_field_flt *f) {
    suNg_flt *site;
    double sqnorm = 0.0;
    _MASTER_FOR(f->type, ix) {
        for (int comp = 0; comp < 4; comp++) {
            float tmp;
            site = _4FIELD_AT(f, ix, comp);
            _suNg_sqnorm(tmp, (*site));
            sqnorm += (double)tmp;
        }
    }

#ifdef WITH_MPI
    global_sum(&sqnorm, 1);
#endif
    return (float)sqnorm;
}

double sqnorm_suNf_field_cpu(suNf_field *f) {
    suNf *site;
    double sqnorm = 0.0;
    _MASTER_FOR(f->type, ix) {
        for (int comp = 0; comp < 4; comp++) {
            double tmp;
            site = _4FIELD_AT(f, ix, comp);
            _suNf_sqnorm(tmp, (*site));
            sqnorm += tmp;
        }
    }

#ifdef WITH_MPI
    global_sum(&sqnorm, 1);
#endif
    return sqnorm;
}

float sqnorm_suNf_field_flt_cpu(suNf_field_flt *f) {
    suNf_flt *site;
    double sqnorm = 0.0;
    _MASTER_FOR(f->type, ix) {
        for (int comp = 0; comp < 4; comp++) {
            float tmp;
            site = _4FIELD_AT(f, ix, comp);
            _suNf_sqnorm(tmp, (*site));
            sqnorm += (double)tmp;
        }
    }

#ifdef WITH_MPI
    global_sum(&sqnorm, 1);
#endif
    return (float)sqnorm;
}

double sqnorm_suNg_scalar_field_cpu(suNg_scalar_field *f) {
    suNg_vector *site;
    double sqnorm = 0.0;
    _MASTER_FOR(f->type, ix) {
        site = _FIELD_AT(f, ix);
        _vector_prod_add_assign_re_g(sqnorm, (*site), (*site));
    }

#ifdef WITH_MPI
    global_sum(&sqnorm, 1);
#endif
    return sqnorm;
}

double sqnorm_suNg_av_field_cpu(suNg_av_field *f) {
    suNg_algebra_vector *site;
    double sqnorm = 0.0;
    _MASTER_FOR(f->type, ix) {
        for (int comp = 0; comp < 4; comp++) {
            double tmp;
            site = _4FIELD_AT(f, ix, comp);
            _algebra_vector_sqnorm_g(tmp, (*site));
            sqnorm += tmp;
        }
    }

#ifdef WITH_MPI
    global_sum(&sqnorm, 1);
#endif
    return sqnorm;
}

double sqnorm_scalar_field_cpu(scalar_field *f) {
    double *site;
    double sqnorm = 0.0;
    _MASTER_FOR(f->type, ix) {
        site = _FIELD_AT(f, ix);
        sqnorm += (*site) * (*site);
    }

#ifdef WITH_MPI
    global_sum(&sqnorm, 1);
#endif
    return sqnorm;
}

double sqnorm_gtransf_cpu(gtransf *f) {
    suNg *site;
    double sqnorm = 0.0;
    _MASTER_FOR(f->type, ix) {
        double tmp;
        site = _FIELD_AT(f, ix);
        _suNg_sqnorm(tmp, (*site));
        sqnorm += tmp;
    }

#ifdef WITH_MPI
    global_sum(&sqnorm, 1);
#endif
    return sqnorm;
}

double sqnorm_ldl_field_cpu(ldl_field *f) {
    ldl_t *site;
    double sqnorm = 0.0;
    _MASTER_FOR(f->type, ix) {
        site = _FIELD_AT(f, ix);
        for (int i = 0; i < NF * (2 * NF + 1); ++i) {
            hr_complex up_comp = (*site).up[i];
            hr_complex dn_comp = (*site).dn[i];
            sqnorm += _complex_re(_complex_prod(up_comp, up_comp));
            sqnorm += _complex_re(_complex_prod(dn_comp, dn_comp));
        }
    }

#ifdef WITH_MPI
    global_sum(&sqnorm, 1);
#endif
    return sqnorm;
}

double sqnorm_clover_term_cpu(clover_term *f) {
    suNfc *site;
    double sqnorm = 0.0;
    _MASTER_FOR(f->type, ix) {
        for (int comp = 0; comp < 4; comp++) {
            double tmp = 0.0;
            site = _4FIELD_AT(f, ix, comp);
            _suNfc_sqnorm(tmp, (*site));
            sqnorm += tmp;
        }
    }

#ifdef WITH_MPI
    global_sum(&sqnorm, 1);
#endif
    return sqnorm;
}

double sqnorm_clover_force_cpu(clover_force *f) {
    suNf *site;
    double sqnorm = 0.0;
    _MASTER_FOR(f->type, ix) {
        for (int comp = 0; comp < 6; comp++) {
            double tmp;
            site = _6FIELD_AT(f, ix, comp);
            _suNf_sqnorm(tmp, (*site));
            sqnorm += tmp;
        }
    }

#ifdef WITH_MPI
    global_sum(&sqnorm, 1);
#endif
    return sqnorm;
}

double sqnorm_spinor_field_cpu(spinor_field *f) {
    suNf_spinor *site;
    double sqnorm = 0.0;
    _MASTER_FOR(f->type, ix) {
        site = _FIELD_AT(f, ix);
        double tmp;
        _spinor_prod_re_f(tmp, (*site), (*site));
        sqnorm += tmp;
    }

#ifdef WITH_MPI
    global_sum(&sqnorm, 1);
#endif
    return sqnorm;
}

float sqnorm_spinor_field_flt_cpu(spinor_field_flt *f) {
    suNf_spinor_flt *site;
    double sqnorm = 0.0;
    _MASTER_FOR(f->type, ix) {
        site = _FIELD_AT(f, ix);
        float tmp;
        _spinor_prod_re_f(tmp, (*site), (*site));
        sqnorm += (double)tmp;
    }

#ifdef WITH_MPI
    global_sum(&sqnorm, 1);
#endif
    return (float)sqnorm;
}

double sqnorm_staple_field_cpu(staple_field *f) {
    suNg *site;
    double sqnorm = 0.0;
    _MASTER_FOR(f->type, ix) {
        for (int comp = 0; comp < 3; comp++) {
            double tmp;
            site = _3FIELD_AT(f, ix, comp);
            _suNg_sqnorm(tmp, (*site));
            sqnorm += tmp;
        }
    }

#ifdef WITH_MPI
    global_sum(&sqnorm, 1);
#endif
    return sqnorm;
}

// Set field to zero
void zero_suNg_field_cpu(suNg_field *f) {
    int len = 4 * f->type->gsize_gauge * sizeof(suNg) / sizeof(double);
    double *dbl_ptr = (double *)(f->ptr);
    for (int i = 0; i < len; ++i) {
        dbl_ptr[i] = 0.0;
    }
}

void zero_scalar_field_cpu(scalar_field *f) {
    int len = f->type->gsize_gauge;
    double *dbl_ptr = f->ptr;
    for (int i = 0; i < len; ++i) {
        dbl_ptr[i] = 0.0;
    }
}

void zero_suNf_field_cpu(suNf_field *f) {
    int len = 4 * f->type->gsize_gauge * sizeof(suNf) / sizeof(double);
    double *dbl_ptr = (double *)(f->ptr);
    for (int i = 0; i < len; ++i) {
        dbl_ptr[i] = 0.0;
    }
}

void zero_suNg_field_flt_cpu(suNg_field_flt *f) {
    int len = 4 * f->type->gsize_gauge * sizeof(suNg_flt) / sizeof(float);
    float *flt_ptr = (float *)(f->ptr);
    for (int i = 0; i < len; ++i) {
        flt_ptr[i] = 0.0f;
    }
}

void zero_suNf_field_flt_cpu(suNf_field_flt *f) {
    int len = 4 * f->type->gsize_gauge * sizeof(suNf_flt) / sizeof(float);
    float *flt_ptr = (float *)(f->ptr);
    for (int i = 0; i < len; ++i) {
        flt_ptr[i] = 0.0f;
    }
}

void zero_suNg_scalar_field_cpu(suNg_scalar_field *f) {
    int len = f->type->gsize_gauge * sizeof(suNg_vector) / sizeof(double);
    double *dbl_ptr = (double *)(f->ptr);
    for (int i = 0; i < len; ++i) {
        dbl_ptr[i] = 0.0;
    }
}

void zero_suNg_av_field_cpu(suNg_av_field *f) {
    int len = 4 * f->type->gsize_gauge * sizeof(suNg_algebra_vector) / sizeof(double);
    double *dbl_ptr = (double *)(f->ptr);
    for (int i = 0; i < len; ++i) {
        dbl_ptr[i] = 0.0;
    }
}

void zero_gtransf_cpu(gtransf *f) {
    int len = f->type->gsize_gauge * sizeof(suNg_field) / sizeof(double);
    double *dbl_ptr = (double *)(f->ptr);
    for (int i = 0; i < len; ++i) {
        dbl_ptr[i] = 0.0;
    }
}

void zero_ldl_field_cpu(ldl_field *f) {
    int len = f->type->gsize_gauge * sizeof(ldl_t) / sizeof(double);
    double *dbl_ptr = (double *)(f->ptr);
    for (int i = 0; i < len; ++i) {
        dbl_ptr[i] = 0.0;
    }
}

void zero_clover_term_cpu(clover_term *f) {
    int len = 4 * f->type->gsize_gauge * sizeof(suNfc) / sizeof(double);
    double *dbl_ptr = (double *)(f->ptr);
    for (int i = 0; i < len; ++i) {
        dbl_ptr[i] = 0.0;
    }
}

void zero_clover_force_cpu(clover_force *f) {
    int len = 6 * f->type->gsize_gauge * sizeof(suNf) / sizeof(double);
    double *dbl_ptr = (double *)(f->ptr);
    for (int i = 0; i < len; ++i) {
        dbl_ptr[i] = 0.0;
    }
}

void zero_staple_field_cpu(staple_field *f) {
    int len = 3 * f->type->gsize_gauge * sizeof(suNg) / sizeof(double);
    double *dbl_ptr = (double *)(f->ptr);
    for (int i = 0; i < len; ++i) {
        dbl_ptr[i] = 0.0;
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

void random_suNf_field_cpu(suNf_field *f) {
    int n = 4 * f->type->gsize_gauge * sizeof(suNf) / sizeof(double);
    ranlxd((double *)(f->ptr), n);
}

void random_suNf_field_flt_cpu(suNf_field_flt *f) {
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
