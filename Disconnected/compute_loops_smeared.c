/*******************************************************************************
 *
 * Compute some disconnected loops
 * Copyright (c) 2014, R. Arthur, V. Drach, A. Hietanen 
 * All rights reserved.
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
/* Disonnected parameters */
typedef struct input_loops {
    char mstring[256];
    double precision;
    int nhits;
    int source_type;
    int n_mom;
    double csw;
    int n_smr;
    double alpha;
    char configlist[256];
    /* for the reading function */
    input_record_t read[10];

} input_loops;

#define init_input_loops(varname)                                                               \
    {                                                                                           \
        .read = {                                                                               \
            { "Fermion mass", "disc:mass = %s", STRING_T, (varname).mstring },                  \
            { "inverter precision", "disc:precision = %lf", DOUBLE_T, &(varname).precision },   \
            { "number of inversions per cnfg", "disc:nhits = %d", INT_T, &(varname).nhits },    \
            { "Source type ", "disc:source_type = %d", INT_T, &(varname).source_type },         \
            { "maximum component of momentum", "disc:n_mom = %d", INT_T, &(varname).n_mom },    \
            { "Configuration list:", "disc:configlist = %s", STRING_T, &(varname).configlist }, \
            { "csw", "disc:csw = %lf", DOUBLE_T, &(varname).csw },                              \
            { "n_smr", "disc:n_smr = %d", INT_T, &(varname).n_smr },                            \
            { "alpha", "disc:alpha = %lf", DOUBLE_T, &(varname).alpha },                        \
            { NULL, NULL, INT_T, NULL }                                                         \
        }                                                                                       \
    }

typedef struct input_HYP {
    /*  int nsteps;*/
    double weight[3];

    /* for the reading function */
    input_record_t read[4];

} input_HYP;

#define init_input_HYP(varname)                                                                  \
    {                                                                                            \
        .read = {                                                                                \
            { "HYP smearing weight[0]", "HYP:weight0 = %lf", DOUBLE_T, &((varname).weight[0]) }, \
            { "HYP smearing weight[1]", "HYP:weight1 = %lf", DOUBLE_T, &((varname).weight[1]) }, \
            { "HYP smearing weight[2]", "HYP:weight2 = %lf", DOUBLE_T, &((varname).weight[2]) }, \
            { NULL, NULL, INT_T, NULL }                                                          \
        }                                                                                        \
    }

char input_filename[256] = "input_file_smeared";
int Nsource;
double M;

enum { UNKNOWN_CNFG, DYNAMICAL_CNFG, QUENCHED_CNFG };

input_loops disc_var = init_input_loops(disc_var);
input_HYP HYP_var = init_input_HYP(HYP_var);

typedef struct {
    char string[256];
    int t, x, y, z;
    int nc, nf;
    double b, m;
    int n;
    int type;
} filename_t;

int main(int argc, char *argv[]) {
    int i;
    FILE *list;
    
    double m[256];
    char list_filename[256] = "";
    char cnfg_filename[256] = "";
    
    static suNg_field *HYP;

    Timer clock;
    timer_set(&clock);
    /* setup process id and communications */

    setup_process(&argc, &argv);
    setup_gauge_fields();

    read_input(disc_var.read, get_input_filename());

    HYP_var.weight[0] = HYP_var.weight[1] = HYP_var.weight[2] = 0.;
    read_input(HYP_var.read, get_input_filename());

    strcpy(list_filename, disc_var.configlist);
    lprintf("MAIN", 0, "list_filename = %s \n", list_filename, disc_var.configlist);

    error((list = fopen(list_filename, "r")) == NULL, 1, "main [measure_spectrum.c]", "Failed to open list file\n");

#if defined(WITH_CLOVER) || defined(WITH_EXPCLOVER)
    set_csw(&disc_var.csw);
#endif

    init_BCs(NULL);

    m[0] = atof(disc_var.mstring); //

    lprintf("MAIN", 0, "HYP parameter smearing: %f, %f, %f\n", HYP_var.weight[0],HYP_var.weight[1],HYP_var.weight[2]);
    lprintf("MAIN", 0, "Inverter precision = %e\n", disc_var.precision);
    lprintf("MAIN", 0, "Mass[%d] = %f\n", 0, m[0]);
    lprintf("MAIN", 0, "nhits = %d\n", disc_var.nhits);
    lprintf("MAIN", 0, "Number of Gaussian smearing levels = %d\n", disc_var.n_smr);
    lprintf("MAIN", 0, "Smearing parameter alpha = %f\n", disc_var.alpha);

    HYP = alloc_suNg_field(&glattice);

    i = 0;
    while (++i) {
        if (list != NULL) {
            if (fscanf(list, "%s", cnfg_filename) == 0 || feof(list)) { break; }
        }

        lprintf("MAIN", 0, "Configuration from %s\n", cnfg_filename);

        // read_gauge_field(cnfg_filename);
        unit_u(u_gauge);

        represent_gauge_field();

        lprintf("TEST", 0, "<p> %1.6f\n", avr_plaquette());
        full_plaquette();

        HYP_smearing(HYP, u_gauge, HYP_var.weight);
        u_gauge = HYP;

        lprintf("TEST", 0, "<p> %1.6f\n", avr_plaquette());
        full_plaquette();

        lprintf("CORR", 0, "Number of noise vector : nhits = %i \n", disc_var.nhits);
        lprintf("CORR", 0, "Number of Gaussian smearing levels = %d\n", 0, disc_var.n_smr);
        lprintf("CORR", 0, "Smearing parameter alpha = %d\n", 0, disc_var.alpha);
        measure_loops_smeared(m, disc_var.nhits, i, disc_var.precision, disc_var.source_type, disc_var.n_mom, disc_var.n_smr,
                              disc_var.alpha, DONTSTORE, NULL);

        if (list == NULL) { break; }
    }

    if (list != NULL) { fclose(list); }

    double elapsed_sec = timer_lap(&clock) * 1.e-6; //time in seconds
    lprintf("TIMING", 0, "Inversions and contractions for configuration  [%s] done [%lf sec]\n", cnfg_filename, elapsed_sec);

    /* close communications */
    free_suNg_field(HYP);
    finalize_process();

    return 0;
}
