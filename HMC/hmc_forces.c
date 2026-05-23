/****************************************************************************
* Copyright (c) 2014, Ari Hietanen
* All rights reserved.                                                      * 
\***************************************************************************/

/*******************************************************************************
*
* Main HMC program
*
*******************************************************************************/

#include "libhr.h"
#include "hmc_utils.h"

/* flow control variable */
hmc_flow flow = init_hmc_flow(flow);

typedef struct input_force_list {
    char configlist[256];
    input_record_t read[2];
} input_force_list;

#define init_configlist(varname)                                                                  \
    {                                                                                             \
        .read = {                                                                                 \
            { "Configuration list:", "forces:configlist = %s", STRING_T, &(varname).configlist }, \
            { NULL, NULL, INT_T, NULL }                                                           \
        }                                                                                         \
    }

input_force_list listinput = init_configlist(listinput);
char list_filename[256] = "";
char cnfg_filename[256] = "";

int main(int argc, char *argv[]) {
    FILE *list = NULL;
    int i = 0;

    setup_process(&argc, &argv);
    setup_gauge_fields();
    init_mc_ghmc(&flow, get_input_filename());

    read_input(listinput.read, get_input_filename());
    strcpy(list_filename, listinput.configlist);

    lprintf("MAIN", 0, "list_filename = %s %s\n", list_filename, listinput.configlist);
    if (strcmp(list_filename, "") != 0) {
        error((list = fopen(list_filename, "r")) == NULL, 1, "main [measure_spectrum.c]", "Failed to open list file\n");
    }

    while (++i) {
        if (list != NULL) {
            if (fscanf(list, "%s", cnfg_filename) == 0 || feof(list)) { break; }
        }

        lprintf("MAIN", 0, "Configuration from %s\n", cnfg_filename);
        read_gauge_field(cnfg_filename);
        represent_gauge_field();
        lprintf("TEST", 0, "<p> %1.6f\n", avr_plaquette());

        Timer clock;
        timer_set(&clock);

        // Measure forces
        print_force_summary();
        double elapsed = timer_lap(&clock) * 1.e-6;
        lprintf("MAIN", 0, "Configuration #%d: analysed in [%lf sec]\n", i, elapsed);
        if (list == NULL) { break; }
    }

    finalize_process();
    return 0;
}
