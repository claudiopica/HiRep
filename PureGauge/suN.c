/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *   
* All rights reserved.                                                      * 
\***************************************************************************/

/*******************************************************************************
*
* Main pure gauge program
*
*******************************************************************************/

#include "suN_utils.h"
#include "libhr.h"
#include <string.h>

pg_flow flow = init_pg_flow(flow);

int main(int argc, char *argv[]) {
    int i, n;
    Timer clock;

    setup_process(&argc, &argv);

    setup_gauge_fields();

    /* Init Monte Carlo */
    init_mc(&flow, get_input_filename());

    timer_set(&clock);
    /* Thermalization */
    for (i = 0; i < flow.therm; ++i) {
        update(&(flow.pg_v->beta), flow.pg_v->nhb, flow.pg_v->nor);
        if (flow.therm > 20) {
            if (i % (flow.therm / 5) == 0) {
                lprintf("MAIN", 0, "%d", ((i * 100) / flow.therm));
            } else if (i % (flow.therm / 20) == 0) {
                lprintf("MAIN", 0, ".");
            }
        }
    }
    if (i) {
        double elapsed_sec = timer_lap(&clock) * 1.e-6; //time in seconds
        lprintf("MAIN", 0, "100\nThermalized %d Trajectories: [%lf sec]\n", flow.therm, elapsed_sec);
        save_conf(&flow, GENERIC_MAX(0, flow.start - 1));
    }
    /* updates and measure  */
    i = flow.start - 1;

    while (i < flow.end) {
        /* nit updates */
        for (n = 0; n < flow.nit; n++) {
            i++;
            lprintf("MAIN", 0, "Trajectory #%d...\n", i);

            double elapsed_sec = timer_lap(&clock) * 1.e-6; //time in seconds

            update(&(flow.pg_v->beta), flow.pg_v->nhb, flow.pg_v->nor);

            elapsed_sec = timer_lap(&clock) * 1.e-6; //time in seconds
            lprintf("MAIN", 0, "Trajectory #%d: generated in [%lf sec]\n", i, elapsed_sec);
            lprintf("MAIN", 0, "Plaquette %1.18e\n", avr_plaquette());

            if ((i % flow.save_freq) == 0) {
                save_conf(&flow, i);
                if (rlx_var.rlxd_state[0] != '\0') {
                    lprintf("MAIN", 0, "Saving rlxd state to file %s\n", rlx_var.rlxd_state);
                    write_ranlxd_state(rlx_var.rlxd_state);
                }
            }
        }

        if (strcmp(flow.wf->make, "true") == 0) {
            static suNg_field *Vwf = NULL;
            if (Vwf == NULL) { Vwf = alloc_suNg_field(&glattice); }
            double elapsed_sec = timer_lap(&clock) * 1.e-6; //time in seconds
            copy_suNg_field(Vwf, u_gauge);
            WF_update_and_measure(flow.wf->ittype, Vwf, &(flow.wf->tmax), &(flow.wf->eps), &(flow.wf->delta), flow.wf->nmeas,
                                  DONTSTORE);
            elapsed_sec = timer_lap(&clock) * 1.e-6; //time in seconds
            lprintf("MAIN", 0, "WF Measure #%d: generated in [%lf sec]\n", i, elapsed_sec);
        }

        if (strcmp(flow.poly->make, "true") == 0) {
            double elapsed_sec = timer_lap(&clock) * 1.e-6; //time in seconds
            polyakov();
            elapsed_sec = timer_lap(&clock) * 1.e-6; //time in seconds
            lprintf("MAIN", 0, "Polyakov Measure #%d: generated in [%lf sec]\n", i, elapsed_sec);
        }
    }

    /* save final configuration */
    if ((i % flow.save_freq) != 0) {
        save_conf(&flow, i);
        /* Only save state if we have a file to save to */
        if (rlx_var.rlxd_state[0] != '\0') {
            lprintf("MAIN", 0, "Saving rlxd state to file %s\n", rlx_var.rlxd_state);
            write_ranlxd_state(rlx_var.rlxd_state);
        }
    }

    /* close communications */
    finalize_process();

    return 0;
}
