/*******************************************************************************
*
* This code takes a list of configuraiton, fix the gauge according to input parameter, and meausre plaq and action.
* Author: R. Arthur ?
*******************************************************************************/

#include "libhr.h"
#include <string.h>

char cnfg_filename[256] = "";
char list_filename[256] = "";
char input_filename[256] = "input_file";
char output_filename[256] = "gaugefix.out";
enum { UNKNOWN_CNFG = 0, DYNAMICAL_CNFG, QUENCHED_CNFG };

/* Gauge fixing parameters */
typedef struct input_gaugefix {
    char make[256];
    int maxit; /* max number of iteration */
    int fixdir; /* fixed direction ? 0, 1, 2, 3 for Coulomb gauge else Landau  */
    double overrelax; /* overrelaxation */
    double tol; /* tolerance */
    char configlist[256]; /* list of configuration */

    /* for the reading function */
    input_record_t read[7];

} input_gaugefix;

#define init_input_gaugefix(varname)                                                                \
    {                                                                                               \
        .read = {                                                                                   \
            { "Max iteration ", "gaugefix:maxit = %d", INT_T, &(varname).maxit },                   \
            { "Fix dir ", "gaugefix:fixdir = %d", INT_T, &(varname).fixdir },                       \
            { "Overrelaxation", "gaugefix:overrelax = %lf", DOUBLE_T, &(varname).overrelax },       \
            { "Tolerance", "gaugefix:tol = %lf", DOUBLE_T, &(varname).tol },                        \
            { "Configuration list:", "gaugefix:configlist = %s", STRING_T, &(varname).configlist }, \
            { "make random gauge transform", "gaugefix:make = %s", STRING_T, (varname).make },      \
            { NULL, NULL, INT_T, NULL }                                                             \
        }                                                                                           \
    }

input_gaugefix gaugefix_var = init_input_gaugefix(gaugefix_var);

typedef struct {
    char string[256];
    int t, x, y, z;
    int nc, nf;
    double b, m;
    int n;
    int type;
} filename_t;

int parse_cnfg_filename(char *filename, filename_t *fn) {
    int hm;
    char *tmp = NULL;
    char *basename;

    basename = filename;
    while ((tmp = strchr(basename, '/')) != NULL) {
        basename = tmp + 1;
    }

    /*#ifdef REPR_FUNDAMENTAL*/
    /*#define repr_name "FUN"*/
    /*#elif defined REPR_SYMMETRIC*/
    /*#define repr_name "SYM"*/
    /*#elif defined REPR_ANTISYMMETRIC*/
    /*#define repr_name "ASY"*/
    /*#elif defined REPR_ADJOINT*/
    /*#define repr_name "ADJ"*/
    /*#endif*/
    hm = sscanf(basename, "%*[^_]_%dx%dx%dx%d%*[Nn]c%dr%*[FSA]%*[UYSD]%*[NMYJ]%*[Nn]f%db%lfm%lfn%d", &(fn->t), &(fn->x),
                &(fn->y), &(fn->z), &(fn->nc), &(fn->nf), &(fn->b), &(fn->m), &(fn->n));
    if (hm == 9) {
        fn->m = -fn->m; /* invert sign of mass */
        fn->type = DYNAMICAL_CNFG;
        return DYNAMICAL_CNFG;
    }
    /*#undef repr_name*/

    double kappa;
    hm = sscanf(basename, "%dx%dx%dx%d%*[Nn]c%d%*[Nn]f%db%lfk%lfn%d", &(fn->t), &(fn->x), &(fn->y), &(fn->z), &(fn->nc),
                &(fn->nf), &(fn->b), &kappa, &(fn->n));
    if (hm == 9) {
        fn->m = .5 / kappa - 4.;
        fn->type = DYNAMICAL_CNFG;
        return DYNAMICAL_CNFG;
    }

    hm = sscanf(basename, "%dx%dx%dx%d%*[Nn]c%db%lfn%d", &(fn->t), &(fn->x), &(fn->y), &(fn->z), &(fn->nc), &(fn->b), &(fn->n));
    if (hm == 7) {
        fn->type = QUENCHED_CNFG;
        return QUENCHED_CNFG;
    }

    hm = sscanf(basename, "%*[^_]_%dx%dx%dx%d%*[Nn]c%db%lfn%d", &(fn->t), &(fn->x), &(fn->y), &(fn->z), &(fn->nc), &(fn->b),
                &(fn->n));
    if (hm == 7) {
        fn->type = QUENCHED_CNFG;
        return QUENCHED_CNFG;
    }

    fn->type = UNKNOWN_CNFG;
    return UNKNOWN_CNFG;
}

int main(int argc, char *argv[]) {
    int i;
    FILE *list;
    filename_t fpars;

    /* setup process communications */
    setup_process(&argc, &argv);

    setup_gauge_fields();

    read_input(gaugefix_var.read, get_input_filename());

    strcpy(list_filename, gaugefix_var.configlist);

    print_compiling_info_short();
    lprintf("MAIN", 0, "PId =  %d [world_size: %d]\n\n", PID, WORLD_SIZE);
    lprintf("MAIN", 0, "input file [%s]\n", input_filename);
    if (strcmp(list_filename, "") != 0) {
        lprintf("MAIN", 0, "list file [%s]\n", list_filename);
    } else {
        error(1, 1, "main [gaugefix_measure.c]", "No list provided");
    }

    lprintf("MAIN", 0, "Gauge Fixing maxit %d \n", gaugefix_var.maxit);
    lprintf("MAIN", 0, "Gauge Fixing fixdir %d\n", gaugefix_var.fixdir); //fixdir= 0, 1, 2, 3 for Coulomb gauge else Landau
    lprintf("MAIN", 0, "Gauge Fixing overrelax %e \n", gaugefix_var.overrelax);
    lprintf("MAIN", 0, "Gauge Fixing tol %e \n", gaugefix_var.tol);
    lprintf("MAIN", 0, "Gauge Fixing make %s\n", gaugefix_var.make);

    suNg_field *fixed_gauge = NULL;
    fixed_gauge = alloc_suNg_field(&glattice);

    list = NULL;
    if (strcmp(list_filename, "") != 0) {
        error((list = fopen(list_filename, "r")) == NULL, 1, "main [gaugefix_measure.c]", "Failed to open list file\n");
    }

    i = 0;
    while (1) {
        if (list != NULL) {
            if (fscanf(list, "%s", cnfg_filename) == 0 || feof(list)) { break; }
        }

        i++;
        /* read & broadcast parameters */
        parse_cnfg_filename(cnfg_filename, &fpars);
        error(fpars.type == UNKNOWN_CNFG, 1, "WF_measure.c", "Bad name for a configuration file");
        error(fpars.nc != NG, 1, "WF_measure.c", "Bad NG");
        error(fpars.t != GLB_T || fpars.x != GLB_X || fpars.y != GLB_Y || fpars.z != GLB_Z, 1, "WF_measure.c",
              "Bad lattice size");

        lprintf("MAIN", 0, "Configuration from %s\n", cnfg_filename);
        read_gauge_field(cnfg_filename);

        lprintf("TEST", 0, "<p> %1.6f\n", avr_plaquette());

        full_plaquette();

        double p1 = calc_plaq(u_gauge);
        lprintf("TEST", 0, "u_gauge plaq %1.6f\n", p1);

        copy_suNg_field(fixed_gauge, u_gauge);

        double p2 = calc_plaq(fixed_gauge);
        lprintf("TEST", 0, "fixed_gauge plaq %1.6f\n", p2);
        if (strcmp(gaugefix_var.make, "true") == 0) {
            random_gauge_transform(fixed_gauge);
            p2 = calc_plaq(fixed_gauge);
            lprintf("TEST", 0, "gt unit_gauge plaq %1.6f\n", p2);
        }
        double act = gaugefix(10, //= 0, 1, 2, 3 for Coulomb gauge else Landau
                              gaugefix_var.overrelax, //overrelax
                              gaugefix_var.maxit, //maxit
                              gaugefix_var.tol, //tolerance
                              fixed_gauge //gauge
        );
        lprintf("TEST", 0, "action  %1.6f\n", act);

        p2 = calc_plaq(fixed_gauge);
        lprintf("TEST", 0, "fixed_gauge plaq %1.6f\n", p2);

        copy_suNg_field(u_gauge, fixed_gauge); //u_gauge = fixed_gauge
        represent_gauge_field(); //u_gauge_f = represented fixed_gauge

        if (list == NULL) { break; }
    }

    if (list != NULL) { fclose(list); }

    finalize_process();

    return 0;
}
