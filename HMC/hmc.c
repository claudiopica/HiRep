/****************************************************************************
* Copyright (c) 2008, Claudio Pica                                          *
* All rights reserved.                                                      *
\***************************************************************************/

/*******************************************************************************
*
* Main HMC program
*
*******************************************************************************/

#include "libhr.h"
#include "hmc_utils.h"
#include <string.h>

/* Mesons parameters */
typedef struct input_mesons {
    char make[256];
    double precision;
    int nhits;
    double mesmass; /* valence mass */

    /* for the reading function */
    input_record_t read[5];

} input_mesons;

#define init_input_mesons(varname)                                                             \
    {                                                                                          \
        .read = {                                                                              \
            { "make mesons", "mes:make = %s", STRING_T, (varname).make },                      \
            { "inverter precision", "mes:precision = %lf", DOUBLE_T, &(varname).precision },   \
            { "number of noisy sources per cnfg", "mes:nhits = %d", INT_T, &(varname).nhits }, \
            { "valence mass", "mes:mass = %lf", DOUBLE_T, &(varname).mesmass },                \
            { NULL, NULL, INT_T, NULL }                                                        \
        }                                                                                      \
    }

input_mesons mes_var = init_input_mesons(mes_var);

/* Polyakov-loop parameters */
typedef struct input_polyakov {
    char make[256];

    /* for the reading function */
    input_record_t read[2];

} input_polyakov;

#define init_input_polyakov(varname)                                                                                   \
    {                                                                                                                  \
        .read = { { "make polyakov loops", "poly:make = %s", STRING_T, (varname).make }, { NULL, NULL, INT_T, NULL } } \
    }

input_polyakov poly_var = init_input_polyakov(poly_var);

typedef struct input_forces {
    char measure[256];
    input_record_t read[2];
} input_forces;

#define init_input_forces(varname)                                                                                        \
    {                                                                                                                     \
        .read = { { "Measure forces", "forces:measure = %s", STRING_T, (varname).measure }, { NULL, NULL, INT_T, NULL } } \
    }

input_forces force_var = init_input_forces(force_var);

/* Lowest-eigenvalue parameters */
typedef struct input_eigval {
    char make[256];
    int nevt; /* search space dimension */
    int nev; /* number of accurate eigenvalues */
    int kmax; /* max degree of polynomial */
    int maxiter; /* max number of subiterations */
    double omega1; /* absolute precision */
    double omega2; /* relative precision */
    double evamass; /* mass to use in Dirac operator */

    /* for the reading function */
    input_record_t read[9];

} input_eigval;

#define init_input_eigval(varname)                                                            \
    {                                                                                         \
        .read = {                                                                             \
            { "make lowest eigenvalues", "eva:make = %s", STRING_T, (varname).make },         \
            { "search space dimension", "eva:nevt = %d", INT_T, &(varname).nevt },            \
            { "number of accurate eigenvalues", "eva:nev = %d", INT_T, &(varname).nev },      \
            { "max degree of polynomial", "eva:kmax = %d", INT_T, &(varname).kmax },          \
            { "max number of subiterations", "eva:maxiter = %d", INT_T, &(varname).maxiter }, \
            { "absolute precision", "eva:omega1 = %lf", DOUBLE_T, &(varname).omega1 },        \
            { "relative precision", "eva:omega2 = %lf", DOUBLE_T, &(varname).omega2 },        \
            { "Dirac op mass", "eva:mass = %lf", DOUBLE_T, &(varname).evamass },              \
            { NULL, NULL, INT_T, NULL }                                                       \
        }                                                                                     \
    }

input_eigval eigval_var = init_input_eigval(eigval_var);

/* flow control variable */
hmc_flow flow = init_hmc_flow(flow);

static void H2eva(spinor_field *out, spinor_field *in) {
    g5Dphi_sq(eigval_var.evamass, out, in);
}

int main(int argc, char *argv[]) {
    int i, acc, rc;

    /* setup process communications */
    setup_process(&argc, &argv);

    setup_gauge_fields();

    /* read input for measures */
    read_input(mes_var.read, get_input_filename());
    read_input(poly_var.read, get_input_filename());
    read_input(eigval_var.read, get_input_filename());
    read_input(force_var.read, get_input_filename());

    /* Init Monte Carlo */

    init_mc_ghmc(&flow, get_input_filename());
    lprintf("MAIN", 0, "MVM during HMC initialization: %ld\n", getMVM());

    lprintf("MAIN", 0, "Initial plaquette: %1.16e\n", avr_plaquette());

    if (strcmp(mes_var.make, "true") == 0) {
        init_meson_correlators(0);
        lprintf("MAIN", 0, "Measuring Gamma Gamma correlators and PCAC-mass\n");
        lprintf("OBSERVABLES", 0, "Inverter precision for mesons = %e\n", mes_var.precision);
        lprintf("OBSERVABLES", 0, "Number of noisy sources for mesons per cnfg = %d\n", mes_var.nhits);
    }

    if (strcmp(eigval_var.make, "true") == 0) {
        lprintf("OBSERVABLES", 0, "EVA Search space dimension  (eva:nevt) = %d\n", eigval_var.nevt);
        lprintf("OBSERVABLES", 0, "EVA Number of accurate eigenvalues (eva:nev) = %d\n", eigval_var.nev);
        lprintf("OBSERVABLES", 0, "EVA Max degree of polynomial (eva:kmax) = %d\n", eigval_var.kmax);
        lprintf("OBSERVABLES", 0, "EVA Max number of subiterations (eva:maxiter) = %d\n", eigval_var.maxiter);
        lprintf("OBSERVABLES", 0, "EVA Absolute precision  (eva:omega1) = %e\n", eigval_var.omega1);
        lprintf("OBSERVABLES", 0, "EVA Relative precision (eva:omega2) = %e\n", eigval_var.omega2);
    }

    double *eva_vals = NULL;
    spinor_field *eva_vecs = NULL;
    if (strcmp(eigval_var.make, "true") == 0) {
        eva_vals = malloc(sizeof(double) * eigval_var.nevt);
        eva_vecs = alloc_spinor_field(eigval_var.nevt, &glattice);
    }
    rc = acc = 0;
    for (i = flow.start; i < flow.end; ++i) {
        int rr;
        double perc;
        lprintf("MAIN", 0, "Trajectory #%d...\n", i);

        Timer clock;
        timer_set(&clock);

        rr = update_ghmc();

        double elapsed_sec = timer_lap(&clock) * 1.e-6; //time in seconds
        lprintf("MAIN", 0, "Trajectory #%d: generated in [%lf sec]\n", i, elapsed_sec);

        if (rr < 0) {
            lprintf("MAIN", 0, "Error in updating the gauge field!!\n");
            return 1;
        } else if (rr != 0) {
            acc++;
        }
        rc++;
        perc = (acc == 0) ? 0. : (float)(100 * acc) / (float)(rc);

        lprintf("MAIN", 0, "Trajectory #%d: %d/%d (%3.4f%%) MVM (f;d) = %ld ; %ld\n", i, acc, rc, perc, getMVM_flt(), getMVM());

        if ((i % flow.save_freq) == 0) {
            save_conf(&flow, i);
            if (u_scalar != NULL) { save_scalar_conf(&flow, i); }
            /* Only save state if we have a file to save to */
            if (rlx_var.rlxd_state[0] != '\0') {
                lprintf("MAIN", 0, "Saving rlxd state to file %s\n", rlx_var.rlxd_state);
                write_ranlxd_state(rlx_var.rlxd_state);
            }
        }

#ifdef WITH_GPU
        copy_from_gpu(u_gauge);
#ifdef ALLOCATE_REPR_GAUGE_FIELD
        copy_from_gpu(u_gauge_f);
#endif
#endif

        if ((i % flow.meas_freq) == 0) {
            /* plaquette */
#ifdef WITH_SMEARING
            lprintf("MAIN", 0, "Plaquette: %1.8e, Smeared: %1.8e\n", avr_plaquette(), avr_smeared_plaquette());
#else
            lprintf("MAIN", 0, "Plaquette: %1.16e\n", avr_plaquette());
            /*avr_ts_plaquette();*/
#endif

            /* Mesons */
            if (strcmp(mes_var.make, "true") == 0) {
                measure_spectrum_semwall(1, &mes_var.mesmass, mes_var.nhits, i, mes_var.precision, DONTSTORE, NULL);
            }

            /* Four fermion observables */
            if (four_fermion_active == 1) { ff_observables(); }

            /* Polyakov loops */
            if (strcmp(poly_var.make, "true") == 0) { polyakov(); }

            if (strcmp(force_var.measure, "true") == 0) { print_force_summary(); }

            /* Lowest eigenvalues */
            if (strcmp(eigval_var.make, "true") == 0) {
                double max;
                max_eigval(&H2eva, &glattice, &max);
                max *= 1.1;
                int status;
                int ie = eva(eigval_var.nev, eigval_var.nevt, 0, eigval_var.kmax, eigval_var.maxiter, max, eigval_var.omega1,
                             eigval_var.omega2, &H2eva, eva_vecs, eva_vals, &status);
                while (ie != 0) { /* if failed restart EVA */
                    lprintf("MAIN", 0, "Restarting EVA!\n");
                    ie = eva(eigval_var.nev, eigval_var.nevt, 2, eigval_var.kmax, eigval_var.maxiter, max, eigval_var.omega1,
                             eigval_var.omega2, &H2eva, eva_vecs, eva_vals, &status);
                }

                for (int n = 0; n < eigval_var.nev; ++n) {
                    lprintf("LOWEIG", 0, "Eig %d = %1.15e\n", n, eva_vals[n]);
                }
            }
        }
    }

    /* save final configuration */
    if (((--i) % flow.save_freq) != 0) {
        save_conf(&flow, i);
        if (u_scalar != NULL) { save_scalar_conf(&flow, i); }
        /* Only save state if we have a file to save to */
        if (rlx_var.rlxd_state[0] != '\0') {
            lprintf("MAIN", 0, "Saving rlxd state to file %s\n", rlx_var.rlxd_state);
            write_ranlxd_state(rlx_var.rlxd_state);
        }
    }

    /* inalize Monte Carlo & close communications */
    finalize_process();

    return 0;
}
