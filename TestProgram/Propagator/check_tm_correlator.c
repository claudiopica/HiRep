/***************************************************************************\
* Copyright (c) 2008-2024, Antonio Rago, Fernando Romero Lopez              *
* All rights reserved.                                                      *
\***************************************************************************/

/******************************************************************************
 *
 * Checks of the tm correlator
 * 
 * NOCOMPILE= BC_X_ANTIPERIODIC
 * NOCOMPILE= BC_Y_ANTIPERIODIC
 * NOCOMPILE= BC_Z_ANTIPERIODIC
 * NOCOMPILE= BC_T_SF
 * NOCOMPILE= BC_T_SF_ROTATED
 * NOCOMPILE= FERMION_THETA
 * 
 ******************************************************************************/

#include "libhr.h"
// #include <string.h>

#define SQR(A) ((A) * (A))
#define CMUL(a, b) \
    (a).re *= b;   \
    (a).im *= b

/// \cond
#define BASENAME(filename) (strrchr((filename), '/') ? strrchr((filename), '/') + 1 : filename)

#define corr_ind(px, py, pz, n_mom, tc, nm, cm) \
    ((px) * (n_mom) * (n_mom) * GLB_T * (nm) + (py) * (n_mom) * GLB_T * (nm) + (pz) * GLB_T * (nm) + ((cm) * GLB_T) + (tc))
/// \endcond

/// \cond
#define INDEX(px, py, pz, n_mom, tc)                                                                         \
    ((px + n_mom) * (2 * n_mom + 1) * (2 * n_mom + 1) * (GLB_T) + (py + n_mom) * (2 * n_mom + 1) * (GLB_T) + \
     (pz + n_mom) * (GLB_T) + (tc))
/// \endcond

/**
 * @brief Structure containing data from the input file relevant to scattering.
 */
typedef struct input_cor {
    char mstring[256];
    char mustring[256];
    double precision;

    /* for the reading function */
    input_record_t read[4];

} input_cor;

#define init_input_cor(varname)                                                               \
    {                                                                                         \
        .read = {                                                                             \
            { "Fermion mass", "prop:mass = %s", STRING_T, (varname).mstring },                \
            { "twisted mass", "prop:mu = %s", STRING_T, (varname).mustring },                 \
            { "inverter precision", "prop:precision = %lf", DOUBLE_T, &(varname).precision }, \
            { NULL, NULL, INT_T, NULL }                                                       \
        }                                                                                     \
    }

char input_filename[256] = "input_file";

input_cor cor_var = init_input_cor(cor_var);

/* Random timeslice not previously chosen */
static int random_tau() {
    static int *slices = NULL;
    if (slices == NULL) { slices = (int *)malloc(GLB_T * sizeof(int)); }
    static int counter = 0;
    int itmp, tau, i;
    double ran;

    if (counter == 0) {
        for (i = 0; i < GLB_T; ++i) {
            slices[i] = i;
        }
        counter = GLB_T;
    }
    do {
        ranlxd(&ran, 1);
        itmp = (int)(ran * counter);
    } while (itmp == counter);
    counter--;
    tau = slices[itmp];
    slices[itmp] = slices[counter];
    slices[counter] = tau;
    bcast_int(&tau, 1);
    return tau;
}

typedef struct fourvector {
    double v[4];
} fourvec;

static void do_global_sum(meson_observable *mo, double norm) {
    meson_observable *motmp = mo;
    int i;
    while (motmp != NULL) {
        global_sum(motmp->corr_re, motmp->corr_size);
        global_sum(motmp->corr_im, motmp->corr_size);
        for (i = 0; i < motmp->corr_size; i++) {
            motmp->corr_re[i] *= norm;
            motmp->corr_im[i] *= norm;
        }
        motmp = motmp->next;
    }
}
/**
 * @brief Adds two four-vectors together replacing the first one with the sum, v1 += v2
 */
static void iadd(fourvec *v1, fourvec *v2) {
    for (int i = 0; i < 4; ++i) {
        v1->v[i] += v2->v[i];
    }
}

/**
 * @brief Multiply four-vector by a real number
 */
static void imul(fourvec *v1, double a) {
    for (int i = 0; i < 4; ++i) {
        v1->v[i] *= a;
    }
}

/**
 * @brief Returns sum over phat^2 + m (part of the propagator)
 */
double f1(fourvec p, double m) {
    double tmp = 0.0;
    int i;

    for (i = 0; i < 4; ++i) {
        tmp += sin(p.v[i] / 2) * sin(p.v[i] / 2);
    }
    return m + 2 * tmp;
}

/**
 * @brief Part of the propagator
 */
static double f2(fourvec v1, fourvec v2) {
    int i;
    double result = 0.0;
    for (i = 0; i < 4; ++i) {
        result += sin(v1.v[i]) * sin(v2.v[i]);
    }
    return result;
}

/**
 * @brief Calculates analytic expression for pion 2-point function
 * @param p Momentum at the sink
 * @param m Quark mass
 * @param L spatial size of the box
 * @param LT time extent of the box
 * @param t time slice
 */
hr_complex tw_twopoint(fourvec p, double m, double mu, int L, int LT, int t) {
    fourvec mom1, mom2;
    int q1, q2, q3, q41, q42;
    hr_complex res;
    res = 0.;
    double tmp;

    for (q1 = 0; q1 < L; ++q1) {
        for (q2 = 0; q2 < L; ++q2) {
            for (q3 = 0; q3 < L; ++q3) {
                for (q41 = 0; q41 < LT; ++q41) {
                    for (q42 = 0; q42 < LT; ++q42) {
#ifdef BC_T_PERIODIC
                        mom1 = (fourvec){ { q1, q2, q3, ((double)q41) * L / LT } };
                        iadd(&mom1, &p);
                        imul(&mom1, 2.0 * PI / L);
#elif BC_T_ANTIPERIODIC
                        mom1 = (fourvec){ { q1 * 2.0 * PI / L, q2 * 2.0 * PI / L, q3 * 2.0 * PI / L,
                                            ((double)(2 * q41 + 1)) * PI / LT } };
#endif

#ifdef BC_T_PERIODIC
                        mom2 = (fourvec){ { q1, q2, q3, ((double)q42) * L / LT } };
                        imul(&mom2, 2.0 * PI / L);
#elif BC_T_ANTIPERIODIC
                        mom2 = (fourvec){ { q1 * 2.0 * PI / L, q2 * 2.0 * PI / L, q3 * 2.0 * PI / L,
                                            ((double)(2 * q42 + 1)) * PI / LT } };
#endif

                        tmp = (f1(mom1, m) * f1(mom2, m) + mu * mu + f2(mom1, mom2)) /
                              ((SQR(f1(mom1, m)) + mu * mu + f2(mom1, mom1)) * (SQR(f1(mom2, m)) + mu * mu + f2(mom2, mom2)));

                        res += tmp * cexp(I * (2.0 * PI / LT) * t * (q42 - q41));
                    }
                }
            }
        }
    }

    res = 4 * res / L / L / L / LT / LT;
    return res;
}

static int compare_corr(hr_complex *corr_ex, hr_complex *corr_num, int tstart, char *name, double tol) {
    int retval = 0;
    for (int t = tstart; t < GLB_T; t++) {
        if (cabs(corr_ex[t] - corr_num[t]) / cabs(corr_ex[t]) > tol) {
            lprintf("TEST", 0, "Mismatch %s, t=%d, relative diff: %e, numeric = %e + I*(%e), analytic = %e + I*(%e) \n", name,
                    t, cabs(corr_ex[t] - corr_num[t]) / cabs(corr_ex[t]), creal(corr_num[t]), cimag(corr_num[t]),
                    creal(corr_ex[t]), cimag(corr_ex[t]));
            retval += 1;
        } else {
            lprintf("TEST", 0, "Match %s, t=%d, numeric = %e + I*(%e), analytic = %e + I*(%e) \n", name, t, creal(corr_num[t]),
                    cimag(corr_num[t]), creal(corr_ex[t]), cimag(corr_ex[t]));
        }
    }
    return retval;
}

static int compare_prop(spinor_field *prop1, spinor_field *prop2, double tol) {
    int retval = 0;
    suNf_spinor *p1, *p2;

    for (int l = 0; l < 4; l++) {
        _MASTER_FOR(prop2->type, ix) {
            p1 = _FIELD_AT(prop1 + l, ix);
            p2 = _FIELD_AT(prop2 + l, ix);
            for (int j = 0; j < 4; j++) {
                for (int k = 0; k < NF; k++) {
                    if (cabs(p1->c[j].c[k] - p2->c[j].c[k]) > tol) {
                        lprintf("ERROR", 0, "%0.8e %0.8e %0.8e %0.8e\n ", creal(p1->c[j].c[k]), cimag(p1->c[j].c[k]),
                                creal(p2->c[j].c[k]), cimag(p2->c[j].c[k]));
                        retval++;
                    }
                }
            }
        }
    }

    global_sum_int(&retval, 1);

    return retval;
}

int main(int argc, char *argv[]) {
    int return_value = 0;
    fourvec zero_p = (fourvec){ { 0, 0, 0, 0 } };
    double mu, m;
    int ts;

    error(!(GLB_X == GLB_Y && GLB_X == GLB_Z), 1, "main", "This test works only for GLB_X=GLB_Y=GLB_Z");
    std_comm_t = ALL_COMMS; // Communications of both the CPU and GPU field copy are necessary
    setup_process(&argc, &argv);
    setup_gauge_fields();

    spinor_field *source = alloc_spinor_field(4, &glattice);
    spinor_field *prop = alloc_spinor_field(4, &glattice);
    spinor_field *propmu0 = alloc_spinor_field(4, &glattice);

    read_input(glb_var.read, get_input_filename());
    read_input(cor_var.read, get_input_filename());
    read_input(rlx_var.read, get_input_filename());

    m = atof(cor_var.mstring);
    mu = atof(cor_var.mustring);
    init_propagator_eo(1, &m, cor_var.precision);
    lprintf("MAIN", 0, "mass is : %e\n", m);

    struct timeval start, end, etime;
    gettimeofday(&start, 0);

    unit_gauge(u_gauge);
    represent_gauge_field();
#ifdef REPR_FUNDAMENTAL
    apply_BCs_on_represented_gauge_field(); // This is a trick: the BCs are not applied in the case the REPR is fundamental because represent_gauge field assumes that the right BCs are already applied on the fundamental field!
#endif

    gettimeofday(&start, 0);
    hr_complex Piseq[GLB_T];
    hr_complex Pith[GLB_T];

    for (int t = 0; t < GLB_T; t++) {
        Piseq[t] = 0.0;
        Pith[t] = tw_twopoint(zero_p, m, mu, GLB_X, GLB_T, t);
    }

    meson_observable mo;
    init_mo(&mo, "pi", GLB_T);
    reset_mo(&mo);

    for (int j = 0; j < 4; j++) {
#ifdef WITH_GPU
        zero_spinor_field_cpu(prop + j);
#endif
        zero_spinor_field(prop + j);
    }

    ts = 0 * random_tau();
    lprintf("MAIN", 0, "ts = %d \n", ts);
    create_point_source(source, ts, 1);
    calc_propagator_tw(&m, 0.0, prop, source, 4);
    create_point_source(source, ts, 1);
    calc_propagator(propmu0, source, 4);

#ifdef WITH_GPU
    for (int beta = 0; beta < 4; beta++) {
        copy_from_gpu(prop + beta);
        copy_from_gpu(propmu0 + beta);
    }
#endif

    return_value = compare_prop(propmu0, prop, 1.e-14);

    lprintf("MAIN", 0, "\nComparing propagators at mu=0\n ");
    if (return_value != 0) {
        lprintf("MAIN", 0, "Check failed\n ");
    } else {
        lprintf("MAIN", 0, "Check passed\n ");
    }

    calc_propagator_tw(&m, mu, prop, source, 4);

#ifdef WITH_GPU
    for (int beta = 0; beta < 4; beta++) {
        copy_from_gpu(prop + beta);
    }
#endif

    // "standard" two points : pi and rho
    measure_mesons_core(prop, prop, source, &mo, 1, ts, 1, 0, GLB_T);
    do_global_sum(&mo, 1.0);

    // now copy the content to a placeholder.
    for (int t = 0; t < GLB_T; t++) {
        Piseq[t] = (mo.corr_re[corr_ind(0, 0, 0, 0, t, 1, 0)] + I * mo.corr_im[corr_ind(0, 0, 0, 0, t, 1, 0)]);
    }

    return_value += compare_corr(Pith, Piseq, 0, "Pi", 1.e-8);

    gettimeofday(&end, 0);
    timeval_subtract(&etime, &end, &start);

    global_sum_int(&return_value, 1);
    lprintf("MAIN", 0, "return_value= %d\n ", return_value);

    finalize_process();

    return return_value;
}
