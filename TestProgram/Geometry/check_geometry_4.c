/*******************************************************************************
*
* Check communication of gauge field in T direction
*
*******************************************************************************/

#include "libhr.h"

int isNotZero(double zero, double tolerance) {
    return (fabs(zero) > tolerance);
}

int checkNorm(char *name, double norm_diff) {
    const double tolerance = (1.e-15 * (GLB_VOLUME) * sizeof(suNf_spinor) / sizeof(double));
    if (isNotZero(norm_diff, tolerance)) {
        lprintf("ERROR", 1, "%s SPINOR norm is wrong [diff=%e >%e]\n", name, norm_diff, tolerance);
        return 1;
    }
    return 0;
}

void spinor_field_one_f(spinor_field *s1) {
    double *p = (double *)s1->ptr;
    int size = s1->type->gsize_spinor * sizeof(*s1->ptr) / sizeof(double);
    for (int i = 0; i < size; i++) {
        p[i] = 1.0;
    }
#ifdef WITH_GPU
    copy_to_gpu_spinor_field_f(s1);
#endif
}

int main(int argc, char *argv[]) {
    int errors = 0;
    setup_process(&argc, &argv);
    atexit(&finalize_process);

#ifndef NDEBUG
    lprintf("test", 1, "GLATTICE\n");
    print_gd(&glattice);
    lprintf("test", 1, "GLAT_EVEN\n");
    print_gd(&glat_even);
    lprintf("test", 1, "GLAT_ODD\n");
    print_gd(&glat_odd);
#endif

    spinor_field *in = alloc_spinor_field_f(1, &glattice);
    spinor_field *even = alloc_spinor_field_f(1, &glat_even);
    spinor_field *odd = alloc_spinor_field_f(1, &glat_odd);

    spinor_field_one_f(in);
    spinor_field_one_f(even);
    spinor_field_one_f(odd);

    start_sendrecv_spinor_field_f(in);
    complete_sendrecv_spinor_field_f(in);
    start_sendrecv_spinor_field_f(even);
    complete_sendrecv_spinor_field_f(even);
    start_sendrecv_spinor_field_f(odd);
    complete_sendrecv_spinor_field_f(odd);

    double in_norm = spinor_field_sqnorm_f(in);
    double even_norm = spinor_field_sqnorm_f(even);
    double odd_norm = spinor_field_sqnorm_f(odd);
    lprintf("TEST", 1, "[START] Norm^2: global=%lf even=%lf odd=%lf\n", in_norm, even_norm, odd_norm);
    double expected_global = ((double)GLB_VOLUME) * sizeof(suNf_spinor) / sizeof(double);
    double expected_eo = expected_global / 2;
    const double tolerance = (1.e-15 * (GLB_VOLUME) * sizeof(suNf_spinor) / sizeof(double));
    if (isNotZero(in_norm - expected_global, tolerance)) {
        lprintf("ERROR", 1, "[GLOBAL] Norm^2 do not match diff=%e\n", in_norm - expected_global);
        errors++;
    }
    if (isNotZero(even_norm - expected_eo, tolerance)) {
        lprintf("ERROR", 1, "[EVEN] Norm^2 do not match diff=%e\n", even_norm - expected_eo);
        errors++;
    }
    if (isNotZero(odd_norm - expected_eo, tolerance)) {
        lprintf("ERROR", 1, "[ODD] Norm^2 do not match diff=%e\n", odd_norm - expected_eo);
        errors++;
    }

    lprintf("TEST", 1, "Errors = %d\n", errors);
    return errors;
}
