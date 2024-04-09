/*******************************************************************************
 *
 * Computation of the Pseudo scalar scattering lengths
 * VD 2014
 * 
 * NOCOMPILE= BC_T_SF_ROTATED || BC_T_SF
 * NOCOMPILE= BC_T_THETA || BC_X_THETA || BC_Y_THETA || BC_Z_THETA
 *******************************************************************************/
#include "libhr.h"
#include <string.h>

#if defined(BC_T_SF_ROTATED) && defined(BC_T_SF)
#error This code does not work with the Schroedinger functional !!!
#endif

#ifdef FERMION_THETA
#error This code does not work with the fermion twisting !!!
#endif

/* Mesons parameters */
typedef struct input_scatt {
    char mstring[256], configlist[256], outpath[256], seq_prop[256];
    double csw;
    double precision;
    int nhits;
    int isospin;

    /* for the reading function */
    input_record_t read[10];

} input_scatt;

#define init_input_scatt(varname)                                                                 \
    {                                                                                             \
        .read = {                                                                                 \
            { "quark quenched masses", "mes:masses = %s", STRING_T, (varname).mstring },          \
            { "csw", "mes:csw = %lf", DOUBLE_T, &(varname).csw },                                 \
            { "inverter precision", "mes:precision = %lf", DOUBLE_T, &(varname).precision },      \
            { "number of inversions per cnfg", "mes:nhits = %d", INT_T, &(varname).nhits },       \
            { "Configuration list:", "mes:configlist = %s", STRING_T, &(varname).configlist },    \
            { "outpath:", "mes:outpath = %s", STRING_T, &(varname).outpath },                     \
            { "Make sequential prop or not", "mes:seq_prop = %s", STRING_T, (varname).seq_prop }, \
            { "Isospin channel", "mes:isospin = %d", INT_T, &(varname).isospin },                 \
            { NULL, NULL, INT_T, NULL }                                                           \
        }                                                                                         \
    }

char cnfg_filename[256] = "";
char list_filename[256] = "";
char source_filename[256] = "";
char input_filename[256] = "input_file";
char output_dir[256] = "./output/";
int Nsource;
double M;

enum { UNKNOWN_CNFG, DYNAMICAL_CNFG, QUENCHED_CNFG };

input_scatt mes_var = init_input_scatt(mes_var);

typedef struct {
    char string[256];
    char configlist[256];
    char outpath[256];
    int t, x, y, z;
    int nc, nf;
    double b, m;
    int n;
    int type;
} filename_t;

int main(int argc, char *argv[]) {
    int k;
    FILE *list;
    int nm;
    double m[256];
    int seq_prop = 0;
    /* setup process communications */
    setup_process(&argc, &argv);

    setup_gauge_fields();

    read_input(glb_var.read, get_input_filename());
    read_input(mes_var.read, get_input_filename());
    read_input(rlx_var.read, get_input_filename());
    if (strcmp(mes_var.seq_prop, "true") == 0) { seq_prop = 1; }
    if (strcmp(mes_var.seq_prop, "none") == 0) { seq_prop = 2; }
    if (strcmp(mes_var.seq_prop, "false") == 0) { seq_prop = 0; }
    lprintf("MAIN", 0, "mes_var.seq_prop = %s, seq_prop=%d (0 = false, 1 = true, 2 = none)\n", mes_var.seq_prop, seq_prop);

    strcpy(list_filename, mes_var.configlist);
    strcpy(output_dir, mes_var.outpath);

    lprintf("MAIN", 0, "list_filename = %s %s\n", list_filename, mes_var.configlist);
    lprintf("MAIN", 0, "output_dir : %s \n", mes_var.outpath);
    lprintf("MAIN", 0, "Isospin channel: %d \n", mes_var.isospin);

    if (strcmp(list_filename, "") != 0) {
        error((list = fopen(list_filename, "r")) == NULL, 1, "main [scattering_lengths.c]", "Failed to open list file\n");
    }

#if defined(WITH_CLOVER) || defined(WITH_EXPCLOVER)
    set_csw(&mes_var.csw);
#endif

    nm = 0;
    m[0] = atof(mes_var.mstring); //
    init_propagator_eo(1, m, mes_var.precision);

    lprintf("MAIN", 0, "Inverter precision = %e\n", mes_var.precision);
    for (k = 0; k < nm; k++) {
        lprintf("MAIN", 0, "Mass[%d] = %f\n", k, m[k]);
    }

    lprintf("CORR", 0, "Number of noise vector : nhits = %i \n", mes_var.nhits);
    while (1) {
        struct timeval start, end, etime;
        gettimeofday(&start, 0);
        if (list != NULL) {
            if (fscanf(list, "%s", cnfg_filename) == 0 || feof(list)) { break; }
        }

        lprintf("MAIN", 0, "Configuration from %s\n", cnfg_filename);
        read_gauge_field(cnfg_filename);
        represent_gauge_field();
        lprintf("TEST", 0, "<p> %1.6f\n", avr_plaquette());
        full_plaquette();
        if (mes_var.isospin == 2) {
            measure_pion_scattering_I2(m, mes_var.nhits, mes_var.precision, mes_var.outpath, cnfg_filename, NULL);
        }
        if (mes_var.isospin == 0) {
            measure_pion_scattering_I0(m, mes_var.nhits, mes_var.precision, mes_var.outpath, cnfg_filename, seq_prop, NULL);
        }
        if (mes_var.isospin != 2 && mes_var.isospin != 0) { lprintf("MAIN", 0, "Isospin channel can be 0 or 2!"); }

        gettimeofday(&end, 0);
        timeval_subtract(&etime, &end, &start);
        lprintf("MAIN", 0, "Configuration : analysed in [%ld sec %ld usec]\n", etime.tv_sec, etime.tv_usec);

        if (list == NULL) { break; }
    }
    lprintf("DEBUG", 0, "ALL done, deallocating\n");

    if (list != NULL) { fclose(list); }

    finalize_process();
    return 0;
}
