/***************************************************************************\
* Copyright (c) 2025, Antonio Rago
* All rights reserved.                                                      * 
\***************************************************************************/

/*******************************************************************************
*
* Main program to measure minimal and maximal Ev of H2
*
*******************************************************************************/

#include "libhr.h"
#include <unistd.h>

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
#if defined(WITH_CLOVER) || defined(WITH_EXPCLOVER)
    double evacsw;
    /* for the reading function */
    input_record_t read[10];
#else
    /* for the reading function */
    input_record_t read[9];
#endif
} input_eigval;

#if defined(WITH_CLOVER) || defined(WITH_EXPCLOVER)
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
            { "Csw value", "eva:csw = %lf", DOUBLE_T, &(varname).evacsw },                    \
            { NULL, NULL, INT_T, NULL }                                                       \
        }                                                                                     \
    }
#else
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
#endif

typedef struct input_poly {
    char make[256];

    /* for the reading function */
    input_record_t read[2];
} input_poly;

#define init_input_poly(varname)                                                                                 \
    {                                                                                                            \
        .read = { { "make Polyakov", "polyakov:make = %s", STRING_T, (varname).make }, { NULL, NULL, 0, NULL } } \
    }

static input_poly poly_var = init_input_poly(poly_var);
static input_eigval eigval_var = init_input_eigval(eigval_var);

typedef struct obs_measure {
    char configlist[256]; /* directory to store gconfs */

    input_eigval *evs;

    input_poly *poly;

    /* for the reading function */
    input_record_t read[2];

} flow_obs_measure;

#define init_flow_obs_measure(varname)                                                                                        \
    {                                                                                                                         \
        .read = { { "Configuration list", "configlist = %s", STRING_T, &(varname).configlist }, { NULL, NULL, INT_T, NULL } } \
    }
static flow_obs_measure var_obs = init_flow_obs_measure(var_obs);

static void H2eva(spinor_field *out, spinor_field *in) {
    g5Dphi_sq(eigval_var.evamass, out, in);
}

int init_mk_obs(flow_obs_measure *gf, char *ifile) {
    gf->evs = &eigval_var;
    gf->poly = &poly_var;
    read_input(var_obs.read, ifile);

    read_input(eigval_var.read, ifile);

    read_input(poly_var.read, ifile);
    lprintf("INIT WF", 0, "Polyakov make=%s\n", poly_var.make);

    BCs_pars_t BCs_pars = { .fermion_twisting_theta = { 0., 0., 0., 0. },
                            .gauge_boundary_improvement_cs = 1.,
                            .gauge_boundary_improvement_ct = 1.,
                            .chiSF_boundary_improvement_ds = 1.,
                            .SF_BCs = 0 };
    init_BCs(&BCs_pars);

    read_input(eigval_var.read, ifile);

#if defined(WITH_CLOVER) || defined(WITH_EXPCLOVER)
    set_csw(&eigval_var.evacsw);
#endif

    return 0;
}

int main(int argc, char *argv[]) {
    FILE *list = NULL;
    char cnfg_filename[256];

    struct timeval start, end, etime; /* //for measurment timing */

    setup_process(&argc, &argv);

    setup_gauge_fields();

    init_mk_obs(&var_obs, get_input_filename());

    /* Measures */
    lprintf("MAIN", 0, "Configurations list from %s\n", var_obs.configlist);

    error((list = fopen(var_obs.configlist, "r")) == NULL, 1, "main [suN_multilevel_measure.c]",
          "Failed to open config list file\n");

    double *eva_vals = NULL;
    spinor_field *eva_vecs = NULL;
    if (strcmp(eigval_var.make, "true") == 0) {
        eva_vals = malloc(sizeof(double) * eigval_var.nevt);
        eva_vecs = alloc_spinor_field(eigval_var.nevt, &glattice);
    }

    while (1) {
        if (fscanf(list, "%s", cnfg_filename) == 0 || feof(list)) { break; }

        char *str;
        int confid;
        str = strrchr(cnfg_filename, 'n');
        error(sscanf(str, "n%d", &confid) != 1, 1, "main [suN_multilevel_measure.c]",
              "Malformed configuration name (not ending by ...n<number>) \n");

        lprintf("MAIN", 0, "\n\nConfiguration %d from %s\n", confid, cnfg_filename);

        read_gauge_field(cnfg_filename);

        apply_BCs_on_fundamental_gauge_field();

        gettimeofday(&start, 0);

        /* Polyakov loops */
        if (strcmp(poly_var.make, "true") == 0) { polyakov(); }

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
        gettimeofday(&end, 0);
        timeval_subtract(&etime, &end, &start);
        lprintf("MAIN", 0, "Obs for conf #%d: generated in [%ld sec %ld usec]\n", confid, etime.tv_sec, etime.tv_usec);
    }

    /* close communications */
    finalize_process();

    return 0;
}
