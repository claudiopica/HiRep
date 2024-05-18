/***************************************************************************\
 * Copyright (c) 2009-2024, Claudio Pica, Sofie Martins                      *
 * All rights reserved.                                                      *
 \***************************************************************************/

/*******************************************************************************
 *
 * Computation of the lowest eigenvalues of H^2
 *
 * NOCOMPILE= BC_T_SF_ROTATED || BC_T_SF
 * NOCOMPILE= BC_T_THETA || BC_X_THETA || BC_Y_THETA || BC_Z_THETA
 *
 *******************************************************************************/

#include "libhr.h"
#include <string.h>

#if defined(BC_T_SF_ROTATED) && defined(BC_T_SF)
#error This code does not work with the Schroedinger functional !!!
#endif

#ifdef FERMION_THETA
#error This code does not work with the fermion twisting !!!
#endif

static double hevamass;

/* Mesons parameters */
typedef struct input_eigval {
    /* EVA parameters */
    int nevt; /* search space dimension */
    int nev; /* number of accurate eigenvalues */
    int kmax; /* max degree of polynomial */
    int maxiter; /* max number of subiterations */
    double omega1; /* absolute precision */
    double omega2; /* relative precision */
    double mass; /* quenched mass */
    char configlist[256];

    /* for the reading function */
    input_record_t read[9];

} input_eigval;

#define init_input_eigval(varname)                                                            \
    {                                                                                         \
        .read = {                                                                             \
            { "search space dimension", "eva:nevt = %d", INT_T, &(varname).nevt },            \
            { "number of accurate eigenvalues", "eva:nev = %d", INT_T, &(varname).nev },      \
            { "max degree of polynomial", "eva:kmax = %d", INT_T, &(varname).kmax },          \
            { "max number of subiterations", "eva:maxiter = %d", INT_T, &(varname).maxiter }, \
            { "absolute precision", "eva:omega1 = %lf", DOUBLE_T, &(varname).omega1 },        \
            { "relative precision", "eva:omega2 = %lf", DOUBLE_T, &(varname).omega2 },        \
            { "quark quenched mass", "eva:mass = %lf", DOUBLE_T, &(varname).mass },           \
	    { "Configuration list", "eva:configlist = %s", STRING_T, &(varname).configlist }, \
            { NULL, NULL, INT_T, NULL }                                                       \
        }                                                                                     \
    }

input_eigval eig_var = init_input_eigval(eig_var);

int main(int argc, char *argv[]) {
    int i, n;
    FILE *list = NULL;

    /* setup process id and communications */
    setup_process(&argc, &argv);
    setup_gauge_fields();
    read_input(glb_var.read, get_input_filename());
    read_input(eig_var.read, get_input_filename());
    lprintf("MAIN", 0, "list file: [%s]\n", eig_var.configlist);
    if (strcmp(eig_var.configlist, "") != 0) {
    	error((list = fopen(eig_var.configlist, "r")) == NULL, 1, "main [mk_eigvals.c]", "Failed to open list file\n");
    }
    init_BCs(NULL);
    hevamass = eig_var.mass;

    lprintf("MAIN", 0, "EVA Parameters:\n");
    lprintf("MAIN", 0, "search space dimension  (eva:nevt) = %d\n", eig_var.nevt);
    lprintf("MAIN", 0, "number of accurate eigenvalues (eva:nev) = %d\n", eig_var.nev);
    lprintf("MAIN", 0, "max degree of polynomial (eva:kmax) = %d\n", eig_var.kmax);
    lprintf("MAIN", 0, "max number of subiterations (eva:maxiter) = %d\n", eig_var.maxiter);
    lprintf("MAIN", 0, "absolute precision  (eva:omega1) = %e\n", eig_var.omega1);
    lprintf("MAIN", 0, "relative precision (eva:omega2) = %e\n", eig_var.omega2);
    lprintf("MAIN", 0, "mass = %f\n", hevamass);

    /* EVA parameters */
    double max, mupp;
    double *eva_val;
    int status, ie;
    spinor_field *eva_vec, *eva_ws;
    eva_val = malloc(sizeof(double) * eig_var.nevt);
#ifndef UPDATE_EO
    eva_vec = alloc_spinor_field(eig_var.nevt + 1, &glattice);
#else
    eva_vec = alloc_spinor_field(eig_var.nevt + 1, &glat_even);
#endif
    eva_ws = eva_vec + eig_var.nevt;
    mupp = fabs(hevamass + 4) + 4;
    mupp *= mupp;

    char cnfg_filename[256];
    Timer clock;
    timer_set(&clock);

    i = 0;
    while (1) {
        if (list != NULL) {
            if (fscanf(list, "%s", cnfg_filename) == 0 || feof(list)) { break; }
        }
        i++;
        lprintf("MAIN", 0, "Configuration from %s\n", cnfg_filename);
        read_gauge_field(cnfg_filename);
        represent_gauge_field();

	timer_lap(&clock);
#ifdef UPDATE_EO
        max_eigval(&H2, &glat_even, &max);
#else
	max_eigval(&H2, &glattice, &max);
#endif
        lprintf("MAIN", 0, "MAXCHECK: cnfg=%e  uppbound=%e diff=%e %s\n", 
			max, mupp, mupp - max, (mupp - max) < 0 ? "[FAILED]" : "[OK]");
        max *= 1.1;

        ie = eva(eig_var.nev, eig_var.nevt, 0, eig_var.kmax, 
			eig_var.maxiter, max, eig_var.omega1, eig_var.omega2, &H2,
                 eva_vec, eva_val, &status);
        while (ie != 0) { /* if failed restart EVA */
            lprintf("MAIN", 0, "Restarting EVA!\n");
            ie = eva(eig_var.nev, eig_var.nevt, 2, eig_var.kmax, eig_var.maxiter, max, eig_var.omega1, eig_var.omega2, &H2,
                     eva_vec, eva_val, &status);
        }

        for (n = 0; n < eig_var.nev; ++n) {
            H(&eva_ws[0], &eva_vec[n]);
            lprintf("RESULT", 0, "Eig %d = %.15e %.15e\n", n, eva_val[n],
                    prod_re_spinor_field(&eva_ws[0], &eva_vec[n]) / sqnorm_spinor_field(&eva_vec[n]));
        }

	double elapsed = timer_lap(&clock) * 1e-6;
	lprintf("TIMING", 0, "Eigenvalue determination for configuration [%s] done [%lf sec]", cnfg_filename, elapsed);
        if (list == NULL) { break; }
    }

    if (list != NULL) { fclose(list); }
    free_BCs();
    free(eva_val);
    free_spinor_field(eva_vec);

    finalize_process();
    return 0;
}
