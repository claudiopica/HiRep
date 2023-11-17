/*******************************************************************************
*
* NOCOMPILE= !WITH_MPI
* NOCOMPILE= WITH_GPU
*
* Check communication of gauge field in T direction
*
*******************************************************************************/

#include "libhr.h"

int isNotZero(double zero, double tolerance) {
    return (fabs(zero) > tolerance);
}

int checkNorm(char *name, double norm_diff) {
    const double tolerance = 1.e-15;
    if (isNotZero(norm_diff, tolerance)) {
        lprintf("ERROR", 1, "%s SPINOR norm is wrong [rel diff=%e >%e]\n", name, norm_diff, tolerance);
        return 1;
    }
    return 0;
}

#define checkForErrors(_LABEL)                                                                                                 \
    start_sendrecv_spinor_field(in);                                                                                           \
    complete_sendrecv_spinor_field(in);                                                                                        \
    start_sendrecv_spinor_field(even);                                                                                         \
    complete_sendrecv_spinor_field(even);                                                                                      \
    start_sendrecv_spinor_field(odd);                                                                                          \
    complete_sendrecv_spinor_field(odd);                                                                                       \
    in_diff = spinor_field_sqnorm_f(in) - in_norm;                                                                             \
    even_diff = spinor_field_sqnorm_f(even) - even_norm;                                                                       \
    odd_diff = spinor_field_sqnorm_f(odd) - odd_norm;                                                                          \
    reldiff = in_diff / in_norm;                                                                                               \
    reldiff_even = even_diff / even_norm;                                                                                      \
    reldiff_odd = odd_diff / odd_norm;                                                                                         \
    lprintf("TEST", 1, _LABEL " Diff Norm^2=global=%.10e even=%.10e odd=%.10e Rel Diff: global=%.10e even=%.10e odd=%.10e \n", \
            in_diff, even_diff, odd_diff, reldiff, reldiff_even, reldiff_odd);                                                 \
    errors += checkNorm("GLOBAL", reldiff);                                                                                    \
    errors += checkNorm("EVEN", reldiff_even);                                                                                 \
    errors += checkNorm("ODD", reldiff_odd)

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

    spinor_field *in = alloc_spinor_field(1, &glattice);
    spinor_field *even = alloc_spinor_field(1, &glat_even);
    spinor_field *odd = alloc_spinor_field(1, &glat_odd);

    gaussian_spinor_field(in);

    in->type = &glat_even;
    spinor_field_copy_f(even, in);

    in->type = &glat_odd;
    in->ptr += glat_odd.master_shift;
    spinor_field_copy_f(odd, in);

    in->type = &glattice;
    in->ptr -= glat_odd.master_shift;

    double in_norm = spinor_field_sqnorm_f(in);
    double even_norm = spinor_field_sqnorm_f(even);
    double odd_norm = spinor_field_sqnorm_f(odd);
    double diff = (in_norm - (even_norm + odd_norm));
    double reldiff = diff / in_norm;
    lprintf("TEST", 1, "[START] Norm^2: global=%.10e even=%.10e odd=%.10e TEST: relative diff=%.10e  diff=%.10e\n", in_norm,
            even_norm, odd_norm, reldiff, diff);
    const double tolerance = 1.e-15;
    if (isNotZero(reldiff, tolerance)) {
        lprintf("ERROR", 1, "[START] Norm^2 do not match rel=%.10e diff=%.10e\n", reldiff, diff);
        errors++;
    }

    double in_diff, even_diff, odd_diff, reldiff_even, reldiff_odd;

    checkForErrors("[SENDRECV]");

    if (NP_T > 1) {
        for (int rep = 0; rep < GLB_T; rep++) {
            //shift field down by one in T direction
            for (int t = 0; t < T; t++) {
                for (int x = 0; x < X; x++) {
                    for (int y = 0; y < Y; y++) {
                        for (int z = 0; z < Z; z++) {
                            int parity = (t + x + y + z + PSIGN) & 1; //parity of dst
                            int dst = ipt(t, x, y, z);
                            int src = ipt(t + 1, x, y, z);
                            *(_FIELD_AT(in, dst)) = *(_FIELD_AT(in, src));
                            if (parity) {
                                *(_FIELD_AT(odd, dst)) = *(_FIELD_AT(even, src));
                            } else {
                                *(_FIELD_AT(even, dst)) = *(_FIELD_AT(odd, src));
                            }
                        }
                    }
                }
            }
            //switch even and odd norms
            double tmp = even_norm;
            even_norm = odd_norm;
            odd_norm = tmp;
            checkForErrors("[T SHIFT]");
        }
    }

    if (NP_X > 1) {
        for (int rep = 0; rep < GLB_X; rep++) {
            //shift field down by one in X direction
            for (int t = 0; t < T; t++) {
                for (int x = 0; x < X; x++) {
                    for (int y = 0; y < Y; y++) {
                        for (int z = 0; z < Z; z++) {
                            int parity = (t + x + y + z + PSIGN) & 1; //parity of dst
                            int dst = ipt(t, x, y, z);
                            int src = ipt(t, x + 1, y, z);
                            *(_FIELD_AT(in, dst)) = *(_FIELD_AT(in, src));
                            if (parity) {
                                *(_FIELD_AT(odd, dst)) = *(_FIELD_AT(even, src));
                            } else {
                                *(_FIELD_AT(even, dst)) = *(_FIELD_AT(odd, src));
                            }
                        }
                    }
                }
            }
            //switch even and odd norms
            double tmp = even_norm;
            even_norm = odd_norm;
            odd_norm = tmp;
            checkForErrors("[X SHIFT]");
        }
    }

    if (NP_Y > 1) {
        for (int rep = 0; rep < GLB_Y; rep++) {
            //shift field down by one in Y direction
            for (int t = 0; t < T; t++) {
                for (int x = 0; x < X; x++) {
                    for (int y = 0; y < Y; y++) {
                        for (int z = 0; z < Z; z++) {
                            int parity = (t + x + y + z + PSIGN) & 1; //parity of dst
                            int dst = ipt(t, x, y, z);
                            int src = ipt(t, x, y + 1, z);
                            *(_FIELD_AT(in, dst)) = *(_FIELD_AT(in, src));
                            if (parity) {
                                *(_FIELD_AT(odd, dst)) = *(_FIELD_AT(even, src));
                            } else {
                                *(_FIELD_AT(even, dst)) = *(_FIELD_AT(odd, src));
                            }
                        }
                    }
                }
            }
            //switch even and odd norms
            double tmp = even_norm;
            even_norm = odd_norm;
            odd_norm = tmp;
            checkForErrors("[Y SHIFT]");
        }
    }

    if (NP_Z > 1) {
        for (int rep = 0; rep < GLB_Z; rep++) {
            //shift field down by one in Z direction
            for (int t = 0; t < T; t++) {
                for (int x = 0; x < X; x++) {
                    for (int y = 0; y < Y; y++) {
                        for (int z = 0; z < Z; z++) {
                            int parity = (t + x + y + z + PSIGN) & 1; //parity of dst
                            int dst = ipt(t, x, y, z);
                            int src = ipt(t, x, y, z + 1);
                            *(_FIELD_AT(in, dst)) = *(_FIELD_AT(in, src));
                            if (parity) {
                                *(_FIELD_AT(odd, dst)) = *(_FIELD_AT(even, src));
                            } else {
                                *(_FIELD_AT(even, dst)) = *(_FIELD_AT(odd, src));
                            }
                        }
                    }
                }
            }
            //switch even and odd norms
            double tmp = even_norm;
            even_norm = odd_norm;
            odd_norm = tmp;
            checkForErrors("[Z SHIFT]");
        }
    }

    lprintf("TEST", 1, "Errors = %d\n", errors);
    return errors;
}
