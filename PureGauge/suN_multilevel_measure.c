/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *   
* All rights reserved.                                                      * 
\***************************************************************************/

/*******************************************************************************
*
* Main pure gauge program
*
*******************************************************************************/

#include "suN_utils_multilevel.h"
#include "libhr.h"
#include <string.h>
#include <ctype.h>

pg_flow_ml_measure flow = init_pg_flow_ml_measure(flow);

int main(int argc, char *argv[]) {
    int i;
    FILE *list = NULL;
    char cnfg_filename[256];
    Timer clock;

    setup_process(&argc, &argv);

    setup_gauge_fields();

    /* Init Monte Carlo */
    init_mc_ml_measure(&flow, get_input_filename());

    /* Measures */
    lprintf("MAIN", 0, "Configurations list from %s\n", flow.configlist);

    error((list = fopen(flow.configlist, "r")) == NULL, 1, "main [suN_multilevel_measure.c]",
          "Failed to open config list file\n");

    while (1) {
        if (fscanf(list, "%s", cnfg_filename) == 0 || feof(list)) { break; }

        char *str;
        str = strrchr(cnfg_filename, 'n');
        error(sscanf(str, "n%d", &i) != 1, 1, "main [suN_multilevel_measure.c]",
              "Malformed configuration name (not ending by ...n<number>) \n");

        lprintf("MAIN", 0, "\n\nConfiguration %d from %s\n", i, cnfg_filename);

        read_gauge_field(cnfg_filename);

        apply_BCs_on_fundamental_gauge_field();

        timer_set(&clock);

        update_hb_multilevel_gb_measure(0, &(flow.pg_v->beta), flow.pg_v->nhb, flow.pg_v->nor, flow.pg_v->ml_niteration,
                                        flow.pg_v->ml_nskip, flow.pg_v->nblkstart, flow.pg_v->nblkend, &(flow.pg_v->APEsmear),
                                        &(flow.pg_v->corrs));

        double elapsed_sec = timer_lap(&clock) * 1.e-6; //time in seconds
        lprintf("MAIN", 0, "ML Measure & update#%d: generated in [%lf sec]\n", i, elapsed_sec);
        lprintf("MAIN", 0, "Plaquette %1.18e\n", avr_plaquette());

        if (strcmp(flow.wf->make, "true") == 0) {
            static suNg_field *Vwf = NULL;
            if (Vwf == NULL) { Vwf = alloc_gfield(&glattice); }
            elapsed_sec = timer_lap(&clock) * 1.e-6; //time in seconds
            suNg_field_copy(Vwf, u_gauge);
            WF_update_and_measure(RK3_ADAPTIVE, Vwf, &(flow.wf->tmax), &(flow.wf->eps), &(flow.wf->delta), flow.wf->nmeas,
                                  DONTSTORE);
            elapsed_sec = timer_lap(&clock) * 1.e-6; //time in seconds
            lprintf("MAIN", 0, "WF Measure #%d: generated in [%lf sec]\n", i, elapsed_sec);
        }

        if (strcmp(flow.poly->make, "true") == 0) {
            elapsed_sec = timer_lap(&clock) * 1.e-6; //time in seconds
            polyakov();
            elapsed_sec = timer_lap(&clock) * 1.e-6; //time in seconds
            lprintf("MAIN", 0, "Polyakov Measure #%d: generated in [%lf sec]\n", i, elapsed_sec);
        }
    }

    /* close communications */
    finalize_process();

    return 0;
}
