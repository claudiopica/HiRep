/*******************************************************************************
 *
 * Test of luscherweisz action and forces
 * NCOMPILE= !WITH_MPI
 *******************************************************************************/

#include "libhr.h"
#include <string.h>
#define COMM (1 == 1)
#define NOCOMM (1 == 0)

int test_wilson_action_and_force(double beta);
int test_lw_force(double beta, double c0, double c1);
int test_ginv_lw_action(double beta, double c0, double c1);
int test_gcov_lw_force(double beta, double c0, double c1);

/*
Test functions:

=== int test_wilson_action_and_force(double beta);
  It tests the LW action and force with c0=1 and c1=0 against Wilson action and force already present in the code.

=== int test_ginv_lw_action(double beta, double c0, double c1);
  It tests the gauge invariance of the LW action.

=== int test_gcov_lw_force(double beta, double c0, double c1);
  It tests the gauge covariance of the LW force.

=== int test_lw_force(double beta, double c0, double c1) {
  It tests the LW force against the numerical derivative of the LW action, point by point in the global lattice.

===  static void random_g(suNg_field* g);
===  static void transform_gauge(suNg_field* gtransf, suNg_field* suNg_field);
===  static void transform_force(suNg_field* gtransf, suNg_av_field* force);
  Utilities for gauge transformations.

*/

int test_wilson_action_and_force(double beta) {
    int ret = 0;
    double s1, s2, diff = 0., err;
    suNg_av_field *f1, *f2, **ff1, **ff2;
    suNg_algebra_vector *v1, *v2;

    ff1 = (suNg_av_field **)malloc(sizeof(suNg_av_field *));
    ff2 = (suNg_av_field **)malloc(sizeof(suNg_av_field *));

    f1 = alloc_suNg_av_field(&glattice);
    ff1[0] = f1;

    f2 = alloc_suNg_av_field(&glattice);
    ff2[0] = f2;

    force_gauge_par par;
    par.momenta = ff1;
    par.c0 = 1;
    par.c1 = 0;
    par.beta = beta;

    force_gauge_par gforce_par;
    gforce_par.momenta = ff2;
    gforce_par.beta = beta;

    calculate_stfld(NOCOMM);

    err = 0.;
    _MASTER_FOR(&glattice, ix) {
        s1 = -(beta / ((double)NG)) * local_plaq(ix);
        s2 = lw_action_density(ix, beta, 1.0, 0.0);
        diff = fabs(s1 - s2);
        if (diff > err) { err = diff; }
    }

    global_max(&err, 1);

    lprintf("TEST", 0, "beta=%f :  Deviation from the standard Wilson action density = %e\n(Should be 10^-15 or so)\n\n", beta,
            err);
    if (diff > 1.e-14) { ret++; }

    _MASTER_FOR(&glattice, ix) {
        for (int mu = 0; mu < 4; mu++) {
            _algebra_vector_zero_g(*_4FIELD_AT(f1, ix, mu));
            _algebra_vector_zero_g(*_4FIELD_AT(f2, ix, mu));
        }
    }

    force0(1., (void *)(&gforce_par));

    lw_force(1., &par);

    err = 0.;
    _MASTER_FOR(&glattice, ix) {
        for (int mu = 0; mu < 4; mu++) {
            v1 = _4FIELD_AT(f1, ix, mu);
            v2 = _4FIELD_AT(f2, ix, mu);
            for (int i = 0; i < NG * NG - 1; i++) {
                diff = fabs(v1->c[i] - v2->c[i]);
                if (diff > err) { err = diff; }
            }
        }
    }
    global_max(&err, 1);

    lprintf("TEST", 0, "beta=%f :  Deviation from the standard Wilson force = %e\n(Should be 10^-15 or so)\n\n", beta, err);
    if (err > 1.e-14) { ret++; }

    free_suNg_av_field(f1);
    free_suNg_av_field(f2);

    free(ff1);
    free(ff2);

    return ret;
}

static void random_g(gtransf *g) {
    _MASTER_FOR(&glattice, ix) {
        random_suNg(_FIELD_AT(g, ix));
    }
    start_sendrecv_gtransf(g);
    complete_sendrecv_gtransf(g);
}

void loc_unit_gauge(gtransf *gauge) {
    _MASTER_FOR(&glattice, ix) {
        _suNg_unit(*_FIELD_AT(gauge, ix));
    }
    start_sendrecv_gtransf(gauge);
    complete_sendrecv_gtransf(gauge);
}

static void transform_gauge(gtransf *gtransf, suNg_field *suNg_field) {
    _MASTER_FOR(&glattice, ix) {
        for (int mu = 0; mu < 4; mu++) {
            int iy = iup(ix, mu);
            suNg *u = _4FIELD_AT(suNg_field, ix, mu);
            suNg v;
            _suNg_times_suNg_dagger(v, *u, *_FIELD_AT(gtransf, iy));
            _suNg_times_suNg(*u, *_FIELD_AT(gtransf, ix), v);
        }
    }
    start_sendrecv_suNg_field(suNg_field);
    complete_sendrecv_suNg_field(suNg_field);
}

static void transform_force(gtransf *gtransf, suNg_av_field *force) {
    _MASTER_FOR(&glattice, ix) {
        for (int mu = 0; mu < 4; mu++) {
            suNg v, m;
            _fund_algebra_represent(m, *_4FIELD_AT(force, ix, mu));
            _suNg_times_suNg_dagger(v, m, *_FIELD_AT(gtransf, ix));
            _suNg_times_suNg(m, *_FIELD_AT(gtransf, ix), v);
            _fund_algebra_project(*_4FIELD_AT(force, ix, mu), m);
        }
    }
}

int test_ginv_lw_action(double beta, double c0, double c1) {
    gtransf *g;
    double *s;
    double diff = 0., err;
    int ret = 0;

    g = alloc_gtransf(&glattice);
    s = malloc(glattice.gsize_gauge * sizeof(double));

    calculate_stfld(NOCOMM);
    _MASTER_FOR(&glattice, ix) {
        s[ix] = lw_action_density(ix, beta, c0, c1);
    }

    random_g(g);
    transform_gauge(g, u_gauge);

    calculate_stfld(NOCOMM);
    err = 0.;
    _MASTER_FOR(&glattice, ix) {
        diff = fabs(s[ix] - lw_action_density(ix, beta, c0, c1));
        if (diff > err) { err = diff; }
    }
    global_max(&err, 1);

    lprintf("TEST", 0, "pars=(%f,%f,%f) :  Gauge invariance LW action = %e\n(Should be 10^-15 or so)\n\n", beta, c0, c1, err);
    if (err > 1.e-14) { ret++; }

    free(s);
    free_gtransf(g);

    return ret;
}

int test_gcov_lw_force(double beta, double c0, double c1) {
    gtransf *g;
    suNg_av_field *f1, *f2, *force;
    suNg_algebra_vector *v1, *v2;
    double diff = 0., err;
    suNg_av_field *momenta[1];

    int ret = 0;

    force = alloc_suNg_av_field(&glattice);
    momenta[0] = force;

    force_gauge_par par;
    par.momenta = momenta;
    par.c0 = c0;
    par.c1 = c1;
    par.beta = beta;

    g = alloc_gtransf(&glattice);
    f1 = alloc_suNg_av_field(&glattice);
    f2 = alloc_suNg_av_field(&glattice);
    _MASTER_FOR(&glattice, ix) {
        for (int mu = 0; mu < 4; mu++) {
            _algebra_vector_zero_g(*_4FIELD_AT(f1, ix, mu));
            _algebra_vector_zero_g(*_4FIELD_AT(f2, ix, mu));
            _algebra_vector_zero_g(*_4FIELD_AT(force, ix, mu));
        }
    }
    random_u(u_gauge);
    start_sendrecv_suNg_field(u_gauge);
    complete_sendrecv_suNg_field(u_gauge);

    lw_force(1., &par);

    copy_suNg_av_field_cpu(f1, force);

    random_g(g);
    transform_force(g, f1);
    transform_gauge(g, u_gauge);

    _MASTER_FOR(&glattice, ix) {
        for (int mu = 0; mu < 4; mu++) {
            _algebra_vector_zero_g(*_4FIELD_AT(force, ix, mu));
        }
    }

    lw_force(1., &par);
    copy_suNg_av_field_cpu(f2, force);

    err = 0.;
    _MASTER_FOR(&glattice, ix) {
        for (int mu = 0; mu < 4; mu++) {
            v1 = _4FIELD_AT(f1, ix, mu);
            v2 = _4FIELD_AT(f2, ix, mu);
            for (int i = 0; i < NG * NG - 1; i++) {
                // lprintf("TEST", 0, "%.10e %.10e\n", v1->c[i], v2->c[i]);
                diff = fabs(v1->c[i] - v2->c[i]);
                if (diff > err) { err = diff; }
            }
        }
    }
    global_max(&err, 1);

    lprintf("TEST", 0, "pars=(%f,%f,%f) :  Gauge covariance LW force = %.10e\n(Should be 10^-15 or so)\n\n", beta, c0, c1, err);
    if (err > 1.e-14) { ret++; }
    free_suNg_av_field(f1);
    free_suNg_av_field(f2);
    free_suNg_av_field(force);
    free_gtransf(g);
    return ret;
}

int test_lw_force(double beta, double c0, double c1) {
    int ret = 0;

    suNg_algebra_vector mom;
    suNg_field *u;
    suNg_av_field *f;
    double *s;
    double eps;
    double err, diff;
    int x[4];

    suNg_av_field *momenta[1];

    f = alloc_suNg_av_field(&glattice);
    momenta[0] = f;

    force_gauge_par par;
    par.momenta = momenta;
    par.c0 = c0;
    par.c1 = c1;
    par.beta = beta;

    u = alloc_suNg_field(&glattice);

    _MASTER_FOR(&glattice, ix) {
        for (int mu = 0; mu < 4; mu++) {
            _algebra_vector_zero_g(*_4FIELD_AT(f, ix, mu));
        }
    }
    s = malloc(glattice.gsize_gauge * sizeof(double));

    lw_force(1., &par);

    calculate_stfld(NOCOMM);
    _MASTER_FOR(&glattice, ix) {
        s[ix] = lw_action_density(ix, beta, c0, c1);
    }

    memcpy(u->ptr, u_gauge->ptr, 4 * glattice.gsize_gauge * sizeof(suNg));

    eps = .000001;
    for (int k = 0; k < 2; k++) {
        err = 0.;
        for (x[0] = 0; x[0] < GLB_T; x[0]++) {
            for (x[1] = 0; x[1] < GLB_X; x[1]++) {
                for (x[2] = 0; x[2] < GLB_Y; x[2]++) {
                    for (x[3] = 0; x[3] < GLB_Z; x[3]++) {
                        int local = (1 == 0);
                        if ((x[0] >= zerocoord[0] && x[0] < zerocoord[0] + T) &&
                            (x[1] >= zerocoord[1] && x[1] < zerocoord[1] + X) &&
                            (x[2] >= zerocoord[2] && x[2] < zerocoord[2] + Y) &&
                            (x[3] >= zerocoord[3] && x[3] < zerocoord[3] + Z)) {
                            local = (1 == 1);
                        }

                        int ix = -1;
                        if (local) {
                            ix = ipt(x[0] - zerocoord[0], x[1] - zerocoord[1], x[2] - zerocoord[2], x[3] - zerocoord[3]);
                        }

                        for (int mu = 0; mu < 4; mu++) {
                            if (local) {
                                gauss((double *)(&mom), NG * NG - 1);
                                ExpX(eps, &mom, pu_gauge(ix, mu));
                            }
                            start_sendrecv_suNg_field(u_gauge);
                            complete_sendrecv_suNg_field(u_gauge);

                            double deltaS = 0.;
                            calculate_stfld(NOCOMM);
                            _MASTER_FOR(&glattice, iy) {
                                deltaS += s[iy] - lw_action_density(iy, beta, c0, c1);
                            }
                            global_sum(&deltaS, 1);

                            double Xf = 0.;
                            if (local) {
                                for (int i = 0; i < NG * NG - 1; i++) {
                                    Xf += _FUND_NORM2 * mom.c[i] * _4FIELD_AT(f, ix, mu)->c[i];
                                }
                            }
                            global_sum(&Xf, 1);

                            diff = fabs(Xf - deltaS / eps);
                            if (diff > err) { err = diff; }

                            memcpy(u_gauge->ptr, u->ptr, 4 * glattice.gsize_gauge * sizeof(suNg));
                            start_sendrecv_suNg_field(u_gauge);
                            complete_sendrecv_suNg_field(u_gauge);
                        }
                    }
                }
            }
        }

        lprintf(
            "TEST", 0,
            "pars=(%f,%f,%f) :  Derivative of the action,  eps = %.3e     fabs(DeltaS - eps*X.force)/eps^2 = %e\n(Should of be order 1)\n\n",
            beta, c0, c1, eps, err / eps);
        if (err > 10) { ret++; }

        eps *= .1;
    }

    free_suNg_field(u);
    free(s);
    return ret;
}

int main(int argc, char *argv[]) {
    int return_value = 0;

    setup_process(&argc, &argv);

    setup_gauge_fields();

    random_u(u_gauge);
    start_sendrecv_suNg_field(u_gauge);
    complete_sendrecv_suNg_field(u_gauge);

    return_value += test_wilson_action_and_force(1.);

    return_value += test_wilson_action_and_force(3.);

    return_value += test_ginv_lw_action(1., 1., 0.);

    return_value += test_ginv_lw_action(1., 0., 1.);

    return_value += test_gcov_lw_force(1., 1., 0.);

    return_value += test_gcov_lw_force(1., 0., 1.);

    return_value += test_lw_force(1., 1., 0.);

    return_value += test_lw_force(1., 0., 1.);

    finalize_process();
    return return_value;
}
