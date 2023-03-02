/*******************************************************************************
 *
 * test generation of the wilson lines
 *
 *******************************************************************************/

#include "libhr.h"

int main(int argc, char *argv[]) {
    wilson_lines *pl;

    int return_value = 0, in, in1, pt, dir, Ldir;
    suNg *w1, *w2, *w3, s1, s2;
    w2 = &s1;
    w3 = &s2;

    hr_complex dop, dop1;
    double max_diff[2], min_val;
    double cop;

    setup_process(&argc, &argv);

    setup_gauge_fields();

    lprintf("MAIN", 0, "Generating a random gauge field... ");
    random_u(u_gauge);
    start_sendrecv_gfield(u_gauge);
    represent_gauge_field();
    lprintf("MAIN", 0, "done.\n\n");
    error(NP_X != 1, 0, "MAIN", "NP_X direction must be 1 for poly_legs measures\n");
    error(NP_Y != 1, 0, "MAIN", "NP_Y direction must be 1 for poly_legs measures\n");
    error(NP_Z != 1, 0, "MAIN", "NP_Z direction must be 1 for poly_legs measures\n");

    max_diff[0] = max_diff[1] = 0.0;

    dir = 1;
    Ldir = X;

    min_val = 10.;

    for (int t = 0; t < T; t++) {
        for (int y = 0; y < Y; y++) {
            for (int z = 0; z < Z; z++) {
                for (int x = 0; x < X; x++) {
                    in1 = ipt(t, x, y, z);
                    pl = polyleg(in1, dir);

                    if (fabs(pl->tr) < min_val) { min_val = fabs(pl->tr); }
                    if (fabs(pl->tr) < 10.e-13) { return_value++; }

                    for (int j = npoly_dist - 1; j >= 0; j--) {
                        in = in1;
                        _suNg_unit(*w2);

                        for (pt = 0; pt < Ldir - j - 1; pt++) {
                            in = idn(in, dir);
                            w1 = pu_gauge_wrk(in, dir);
                            _suNg_times_suNg(*w3, *w1, *w2);
                            w1 = w2;
                            w2 = w3;
                            w3 = w1;
                        }
                        _suNg_sub_assign(*w2, pl->p[j]);
                        _suNg_sqnorm(dop, *w2);

                        cop = sqrt(creal(dop));

                        if (cop > max_diff[0]) { max_diff[0] = cop; }

                        if (cop > 10.e-13) { return_value++; }
                    }
                    in = in1;
                    _suNg_unit(*w2);

                    for (pt = 0; pt < Ldir; pt++) {
                        in = idn(in, dir);
                        w1 = pu_gauge_wrk(in, dir);
                        _suNg_times_suNg(*w3, *w1, *w2);
                        w1 = w2;
                        w2 = w3;
                        w3 = w1;
                    }
                    _suNg_trace(dop, *w2);
                    dop -= pl->tr;
                    _complex_mul_star(dop1, dop, dop);
                    cop = sqrt(creal(dop1));

                    if (cop > max_diff[1]) { max_diff[1] = cop; }

                    if (cop > 10.e-13) { return_value++; }
                }
            }
        }
    }
    global_max(max_diff, 2);
    global_min(&min_val, 1);

    lprintf("MAIN", 0, "Maximal trace real difference in X= %.16e\n", max_diff[1]);
    lprintf("MAIN", 0, "Maximal normalized leg difference in X= %.16e\n", max_diff[0]);
    lprintf("MAIN", 0, "Minimal amplitude for the trace of the polyakov line in X= %.16e\n", min_val);
    max_diff[0] = max_diff[1] = 0.0;

    dir = 2;
    Ldir = Y;
    min_val = 10.;
    for (int t = 0; t < T; t++) {
        for (int z = 0; z < Z; z++) {
            for (int x = 0; x < X; x++) {
                for (int y = 0; y < Y; y++) {
                    in1 = ipt(t, x, y, z);
                    pl = polyleg(in1, dir);

                    if (fabs(pl->tr) < min_val) { min_val = fabs(pl->tr); }
                    if (fabs(pl->tr) < 10.e-13) { return_value++; }

                    for (int j = npoly_dist - 1; j >= 0; j--) {
                        in = in1;
                        _suNg_unit(*w2);

                        for (pt = 0; pt < Ldir - j - 1; pt++) {
                            in = idn(in, dir);
                            w1 = pu_gauge_wrk(in, dir);
                            _suNg_times_suNg(*w3, *w1, *w2);
                            w1 = w2;
                            w2 = w3;
                            w3 = w1;
                        }

                        _suNg_sub_assign(*w2, pl->p[j]);
                        _suNg_sqnorm(dop, *w2);

                        cop = sqrt(creal(dop));

                        if (cop > max_diff[0]) { max_diff[0] = cop; }

                        if (cop > 10.e-13) { return_value++; }
                    }
                    in = in1;
                    _suNg_unit(*w2);

                    for (pt = 0; pt < Ldir; pt++) {
                        in = idn(in, dir);
                        w1 = pu_gauge_wrk(in, dir);
                        _suNg_times_suNg(*w3, *w1, *w2);
                        w1 = w2;
                        w2 = w3;
                        w3 = w1;
                    }
                    _suNgc_trace(dop, *w2);
                    dop -= pl->tr;
                    _complex_mul_star(dop1, dop, dop);
                    cop = sqrt(creal(dop1));

                    if (cop > max_diff[1]) { max_diff[1] = cop; }

                    if (cop > 10.e-13) { return_value++; }
                }
            }
        }
    }
    global_max(max_diff, 2);
    global_min(&min_val, 1);

    lprintf("MAIN", 0, "Maximal trace real difference in Y= %.16e\n", max_diff[1]);
    lprintf("MAIN", 0, "Maximal normalized leg difference in Y= %.16e\n", max_diff[0]);
    lprintf("MAIN", 0, "Minimal amplitude for the trace of the polyakov line in Y= %.16e\n", min_val);

    max_diff[0] = max_diff[1] = 0.0;

    dir = 3;
    Ldir = Z;

    for (int t = 0; t < T; t++) {
        for (int x = 0; x < X; x++) {
            for (int y = 0; y < Y; y++) {
                for (int z = 0; z < Z; z++) {
                    in1 = ipt(t, x, y, z);
                    pl = polyleg(in1, dir);

                    if (fabs(pl->tr) < min_val) { min_val = fabs(pl->tr); }
                    if (fabs(pl->tr) < 10.e-13) { return_value++; }

                    for (int j = npoly_dist - 1; j >= 0; j--) {
                        in = in1;
                        _suNg_unit(*w2);

                        for (pt = 0; pt < Ldir - j - 1; pt++) {
                            in = idn(in, dir);
                            w1 = pu_gauge_wrk(in, dir);
                            _suNg_times_suNg(*w3, *w1, *w2);
                            w1 = w2;
                            w2 = w3;
                            w3 = w1;
                        }

                        _suNg_sub_assign(*w2, pl->p[j]);
                        _suNg_sqnorm(dop, *w2);

                        cop = sqrt(creal(dop));

                        if (cop > max_diff[0]) { max_diff[0] = cop; }

                        if (cop > 10.e-13) { return_value++; }
                    }
                    in = in1;
                    _suNg_unit(*w2);

                    for (pt = 0; pt < Ldir; pt++) {
                        in = idn(in, dir);
                        w1 = pu_gauge_wrk(in, dir);
                        _suNg_times_suNg(*w3, *w1, *w2);
                        w1 = w2;
                        w2 = w3;
                        w3 = w1;
                    }
                    _suNgc_trace(dop, *w2);
                    dop -= pl->tr;
                    _complex_mul_star(dop1, dop, dop);
                    cop = sqrt(creal(dop1));

                    if (cop > max_diff[1]) { max_diff[1] = cop; }

                    if (cop > 10.e-13) { return_value++; }
                }
            }
        }
    }

    global_max(max_diff, 2);
    global_min(&min_val, 1);

    lprintf("MAIN", 0, "Maximal trace real difference in Z= %.16e\n", max_diff[1]);
    lprintf("MAIN", 0, "Maximal normalized leg difference in Z= %.16e\n", max_diff[0]);
    lprintf("MAIN", 0, "Minimal amplitude for the trace of the polyakov line in Z= %.16e\n", min_val);

    global_sum_int(&return_value, 1);

    finalize_process();

    return return_value;
}
