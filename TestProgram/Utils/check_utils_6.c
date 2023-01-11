/*******************************************************************************
 *
 * Check of AVX2 routine
 * NOCOMPILE= !AVX2_HIREP
 *
 *******************************************************************************/

#include "libhr.h"

int main(int argc, char *argv[])
{

    int return_value = 0;

    spinor_field *ss1;
    suNf_spinor *s1, *s2;
    suNf_vector v3, v4, v5, v6;

    setup_process(&argc, &argv);

    setup_gauge_fields();
    ss1 = alloc_spinor_field_f(2, &glattice);

    lprintf("MAIN", 0, "Generating a random gauge field... ");
    random_u(u_gauge);
    start_gf_sendrecv(u_gauge);
    represent_gauge_field();
    lprintf("MAIN", 0, "done.\n\n");

    lprintf("MAIN", 0, "Generating a random spinor fields... ");
    gaussian_spinor_field(ss1);
    gaussian_spinor_field(ss1 + 1);
    lprintf("MAIN", 0, "done.\n\n");

    double err = 0, diff;
    _MASTER_FOR(&glattice, ix)
    {
        s1 = _FIELD_AT(ss1, ix);
        s2 = _FIELD_AT(ss1 + 1, ix);
        for (int mu = 0; mu < 4; mu++)
        {
            for (int i = 0; i < 4; i++)
            {

                _suNf_double_multiply(v3, v4, *pu_gauge_f(ix, mu), s1->c[i], s2->c[i]);

                _suNf_double_multiply_default(v5, v6, *pu_gauge_f(ix, mu), s1->c[i], s2->c[i]);

                for (int j = 0; j < NF; j++)
                {
                    diff = cabs(v3.c[j] - v5.c[j]);
                    if (diff > err)
                        err = diff;

                    diff = cabs(v4.c[j] - v6.c[j]);
                    if (diff > err)
                        err = diff;
                }

            }
        }
    }
    global_max(&err, 1);

    lprintf("MAIN", 0, "Checking difference between _suNf_double_multiply_default and _suNf_double_multiply_AVX2.\n ");
    lprintf("MAIN", 0, "Maximal difference = %.14e\n", cabs(err));
    lprintf("MAIN", 0, "(should be around 1*10^(-15) or so)\n\n");
    if (cabs(err) > 10.e-14)
        return_value++;

    finalize_process();
    return return_value;
}
