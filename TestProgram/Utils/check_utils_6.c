/*******************************************************************************
 *
 * Check of SIMD routines
 * NOCOMPILE= !AVX2_HIREP && !SIMD_VECTOR_HIREP
 *
 *******************************************************************************/

#include "libhr.h"

int main(int argc, char *argv[])
{
    int return_value = 0;

    spinor_field *ss1;
    suNf_spinor *s1, *s2;
    suNf_vector v3, v4, v5, v6;
    suNg_vector w3, w4, w5, w6, *f1, *f2;
    suNg_scalar_field *ff1, *ff2;

    setup_process(&argc, &argv);

    setup_gauge_fields();
    ss1 = alloc_spinor_field_f(2, &glattice);
    ff1 = alloc_suNg_scalar_field(&glattice);
    ff2 = alloc_suNg_scalar_field(&glattice);

    lprintf("MAIN", 0, "Generating a random gauge field... ");
    random_u(u_gauge);
    start_gf_sendrecv(u_gauge);
    represent_gauge_field();
    lprintf("MAIN", 0, "done.\n\n");

    lprintf("MAIN", 0, "Generating two random spinor fields... ");
    gaussian_spinor_field(ss1);
    gaussian_spinor_field(ss1 + 1);
    lprintf("MAIN", 0, "done.\n\n");

    lprintf("MAIN", 0, "Generating two random suNg scalar fields... ");
    random_s(ff1);
    random_s(ff2);
    lprintf("MAIN", 0, "done.\n\n");

    lprintf("MAIN", 0, "Evaluating the difference of the AVX operators versus the default ones over the entire lattice.\n\n");

    double err[8] = { 0 }, diff;
    _MASTER_FOR(&glattice, ix) {
        s1 = _FIELD_AT(ss1, ix);
        s2 = _FIELD_AT(ss1 + 1, ix);
        f1 = _FIELD_AT(ff1, ix);
        f2 = _FIELD_AT(ff2, ix);

        for (int mu = 0; mu < 4; mu++) {
            for (int i = 0; i < 4; i++) {
                _suNf_multiply(v3, *pu_gauge_f(ix, mu), s1->c[i]);

                _suNf_multiply_default(v5, *pu_gauge_f(ix, mu), s1->c[i]);

                for (int j = 0; j < NF; j++) {
                    diff = cabs(v3.c[j] - v5.c[j]);
                    if (diff > err[0]) err[0] = diff;
                }

                _suNf_inverse_multiply(v3, *pu_gauge_f(ix, mu), s1->c[i]);

                _suNf_inverse_multiply_default(v5, *pu_gauge_f(ix, mu), s1->c[i]);

                for (int j = 0; j < NF; j++) {
                    diff = cabs(v3.c[j] - v5.c[j]);
                    if (diff > err[1]) err[1] = diff;
                }

                _suNf_double_multiply(v3, v4, *pu_gauge_f(ix, mu), s1->c[i], s2->c[i]);

                _suNf_double_multiply_default(v5, v6, *pu_gauge_f(ix, mu), s1->c[i], s2->c[i]);

                for (int j = 0; j < NF; j++) {
                    diff = cabs(v3.c[j] - v5.c[j]);
                    if (diff > err[2]) err[2] = diff;

                    diff = cabs(v4.c[j] - v6.c[j]);
                    if (diff > err[2]) err[2] = diff;
                }

                _suNf_double_inverse_multiply(v3, v4, *pu_gauge_f(ix, mu), s1->c[i], s2->c[i]);

                _suNf_double_inverse_multiply_default(v5, v6, *pu_gauge_f(ix, mu), s1->c[i], s2->c[i]);

                for (int j = 0; j < NF; j++) {
                    diff = cabs(v3.c[j] - v5.c[j]);
                    if (diff > err[3]) err[3] = diff;

                    diff = cabs(v4.c[j] - v6.c[j]);
                    if (diff > err[3]) err[3] = diff;
                }
            }

            _suNg_multiply(w3, *pu_gauge_f(ix, mu), *f1);

            _suNg_multiply_default(w5, *pu_gauge_f(ix, mu), *f1);

            for (int j = 0; j < NG; j++) {
                diff = cabs(w3.c[j] - w5.c[j]);
                if (diff > err[4]) err[4] = diff;
            }

            _suNg_inverse_multiply(w3, *pu_gauge_f(ix, mu), *f1);

            _suNg_inverse_multiply_default(w5, *pu_gauge_f(ix, mu), *f1);

            for (int j = 0; j < NG; j++) {
                diff = cabs(w3.c[j] - w5.c[j]);
                if (diff > err[5]) err[5] = diff;
            }

            _suNg_double_multiply(w3, w4, *pu_gauge_f(ix, mu), *f1, *f2);

            _suNg_double_multiply_default(w5, w6, *pu_gauge_f(ix, mu), *f1, *f2);

            for (int j = 0; j < NG; j++) {
                diff = cabs(w3.c[j] - w5.c[j]);
                if (diff > err[6]) err[6] = diff;

                diff = cabs(w4.c[j] - w6.c[j]);
                if (diff > err[6]) err[6] = diff;
            }

            _suNg_double_inverse_multiply(w3, w4, *pu_gauge_f(ix, mu), *f1, *f2);

            _suNg_double_inverse_multiply_default(w5, w6, *pu_gauge_f(ix, mu), *f1, *f2);

            for (int j = 0; j < NG; j++) {
                diff = cabs(w3.c[j] - w5.c[j]);
                if (diff > err[7]) err[7] = diff;

                diff = cabs(w4.c[j] - w6.c[j]);
                if (diff > err[7]) err[7] = diff;
            }
        }
    }
    global_max(err, 8);

    lprintf("MAIN", 0, "Checking difference between _suNf_multiply_default and the SIMD _suNf_multiply.\n ");
    lprintf("MAIN", 0, "Maximal difference = %.14e\n", fabs(err[0]));
    lprintf("MAIN", 0, "(should be around 1*10^(-15) or so)\n\n");
    if (fabs(err[0]) > 10.e-14) return_value++;

    lprintf("MAIN", 0, "Checking difference between _suNf_inverse_multiply_default and SIMD _suNf_inverse_multiply.\n ");
    lprintf("MAIN", 0, "Maximal difference = %.14e\n", fabs(err[1]));
    lprintf("MAIN", 0, "(should be around 1*10^(-15) or so)\n\n");
    if (fabs(err[1]) > 10.e-14) return_value++;

    lprintf("MAIN", 0, "Checking difference between _suNf_double_multiply_default and SIMD _suNf_double_multiply.\n ");
    lprintf("MAIN", 0, "Maximal difference = %.14e\n", fabs(err[2]));
    lprintf("MAIN", 0, "(should be around 1*10^(-15) or so)\n\n");
    if (fabs(err[2]) > 10.e-14) return_value++;

    lprintf("MAIN", 0,
            "Checking difference between _suNf_double_inverse_multiply_default and SIMD _suNf_double_inverse_multiply.\n ");
    lprintf("MAIN", 0, "Maximal difference = %.14e\n", fabs(err[3]));
    lprintf("MAIN", 0, "(should be around 1*10^(-15) or so)\n\n");
    if (fabs(err[3]) > 10.e-14) return_value++;

    lprintf("MAIN", 0, "Checking difference between _suNg_multiply_default and SIMD _suNg_multiply.\n ");
    lprintf("MAIN", 0, "Maximal difference = %.14e\n", fabs(err[4]));
    lprintf("MAIN", 0, "(should be around 1*10^(-15) or so)\n\n");
    if (fabs(err[4]) > 10.e-14) return_value++;

    lprintf("MAIN", 0, "Checking difference between _suNg_inverse_multiply_default and SIMD _suNg_inverse_multiply.\n ");
    lprintf("MAIN", 0, "Maximal difference = %.14e\n", fabs(err[5]));
    lprintf("MAIN", 0, "(should be around 1*10^(-15) or so)\n\n");
    if (fabs(err[1]) > 10.e-14) return_value++;

    lprintf("MAIN", 0, "Checking difference between _suNg_double_multiply_default and SIMD _suNg_double_multiply.\n ");
    lprintf("MAIN", 0, "Maximal difference = %.14e\n", fabs(err[6]));
    lprintf("MAIN", 0, "(should be around 1*10^(-15) or so)\n\n");
    if (fabs(err[2]) > 10.e-14) return_value++;

    lprintf("MAIN", 0,
            "Checking difference between _suNg_double_inverse_multiply_default and SIMD _suNg_double_inverse_multiply.\n ");
    lprintf("MAIN", 0, "Maximal difference = %.14e\n", fabs(err[7]));
    lprintf("MAIN", 0, "(should be around 1*10^(-15) or so)\n\n");
    if (fabs(err[3]) > 10.e-14) return_value++;
    finalize_process();
    return return_value;
}
