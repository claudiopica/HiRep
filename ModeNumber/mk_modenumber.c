/***************************************************************************\
 * Copyright (c) 2009-2024, Claudio Pica, Sofie Martins                      *
 * All rights reserved.                                                      *
 \***************************************************************************/

/*******************************************************************************
 *
 * NOCOMPILE= BC_T_SF_ROTATED || BC_T_SF
 * NOCOMPILE= BC_T_THETA || BC_X_THETA || BC_Y_THETA || BC_Z_THETA
 *
 *******************************************************************************/

#include "libhr.h"
#include "modenumber.h"
#include <string.h>

#if defined(BC_T_SF_ROTATED) && defined(BC_T_SF)
#error This code does not work with the Schroedinger functional !!!
#endif

#ifdef FERMION_THETA
#error This code does not work with the fermion twisting !!!
#endif

static double hevamass;

typedef struct input_nu {
    double inverr2;
    char approx[512];
    int nhits;
    double mass;
    double csw; /* clover coefficient: required with WITH_CLOVER / WITH_EXPCLOVER */
    char list[1024];
    char configlist[256];

    /* for the reading function */
    input_record_t read[8];

} input_nu;

#define init_input_nu(varname)                                                                      \
    {                                                                                               \
        .csw = NAN,                                                                                 \
        .read = {                                                                                   \
            { "squared error for inverter", "nu:inverr2 = %lf", DOUBLE_T, &(varname).inverr2 },     \
            { "Chebyshev approximation file", "nu:approx = %s", STRING_T, (varname).approx },       \
            { "number of stochastic spinors", "nu:nhits = %d", INT_T, &(varname).nhits },           \
            { "bare quark mass am0", "nu:mass = %lf", DOUBLE_T, &(varname).mass },                   \
            { "clover coefficient", "nu:csw = %lf", DOUBLE_T, &(varname).csw },                      \
            { "list of eigenvalues", "nu:list = %s", STRING_T, (varname).list },                    \
            { "Configuration list", "nu:configlist = %s", STRING_T, &(varname).configlist },        \
            { NULL, NULL, INT_T, NULL }                                                             \
        }                                                                                           \
    }

input_nu nu_var = init_input_nu(nu_var);

int main(int argc, char *argv[]) {
    char tmp[1024];
    FILE *list = NULL;
    char *cptr;
    int neig;
    double M[1024];

    /* setup process id and communications */
    setup_process(&argc, &argv);
    setup_gauge_fields();
    read_input(nu_var.read, get_input_filename());
    hevamass = nu_var.mass;
#if defined(WITH_CLOVER) || defined(WITH_EXPCLOVER)
    error(isnan(nu_var.csw), 1, "main [mk_modenumber.c]", "nu:csw missing from the input file (required for clover)\n");
    set_csw(&nu_var.csw);
    lprintf("MAIN", 0, "csw = %.10f\n", get_csw());
#endif
    strcpy(tmp, nu_var.list);
    lprintf("MAIN", 0, "list file: [%s]\n", nu_var.configlist);
    if (strcmp(nu_var.configlist, "") != 0) {
        error((list = fopen(nu_var.configlist, "r")) == NULL, 1, "main [mk_modenumber.c]", "Failed to open list file\n");
    }
    cptr = strtok(tmp, ";");
    neig = 0;
    while (cptr != NULL) {
        M[neig] = atof(cptr);
        neig++;
        cptr = strtok(NULL, ";");
    }
    error(neig == 0, 1, "mk_modenumber.c", "neig == 0 !!!");
    init_BCs(NULL);
    init_modenumber(nu_var.mass, nu_var.inverr2, nu_var.nhits, nu_var.approx);
    for (int k = 0; k < neig; k++) {
        lprintf("MODENUMBER", 0, "M[%d] = %e\n", k, M[k]);
    }

    char cnfg_filename[256];
    Timer clock;
    timer_set(&clock);

    int i = 0;
    while (1) {
        if (list != NULL) {
            if (fscanf(list, "%s", cnfg_filename) == 0 || feof(list)) { break; }
        }
        i++;
        lprintf("MAIN", 0, "Configuration from %s\n", cnfg_filename);
        read_gauge_field(cnfg_filename);
        represent_gauge_field();

        timer_lap(&clock);
        for (int k = 0; k < neig; k++) {
            double number = ModeNumber(M[k] * M[k]);
            int mvm = getMVM();
            lprintf("MODENUMBER", 0, "nu[ %e ] = %.2f\n", M[k], number);
            lprintf("MODENUMBER", 0, "MVM = %d\n", mvm);
        }

        double elapsed = timer_lap(&clock) * 1e-6;
        lprintf("TIMING", 0, "Mode number determination for configuration [%s] done [%lf sec]", cnfg_filename, elapsed);
        if (list == NULL) { break; }
    }

    if (list != NULL) { fclose(list); }
    free_BCs();
    free_modenumber();

    finalize_process();
    return 0;
}
