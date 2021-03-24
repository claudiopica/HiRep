/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *   
* All rights reserved.                                                      * 
\***************************************************************************/

/*******************************************************************************
*
* Main pure gauge program
*
*******************************************************************************/

#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include "io.h"
#include "ranlux.h"
#include "geometry.h"
#include "update.h"
#include "global.h"
#include "observables.h"
#include "dirac.h"
#include "logger.h"
#include "memory.h"
#include "communications.h"
#include "observables.h"
#include "utils.h"
#include "suN_utils.h"
#include "setup.h"
#include "wilsonflow.h"

pg_flow flow = init_pg_flow(flow);

int main(int argc, char *argv[])
{
    int i, n;
    struct timeval start, end, etime; /* //for trajectory timing */

    setup_process(&argc, &argv);

    setup_gauge_fields();

    /* Init Monte Carlo */
    init_mc(&flow, get_input_filename());

    /* Thermalization */
    gettimeofday(&start, 0);
    for (i = 0; i < flow.therm; ++i)
    {

        update(&(flow.pg_v->beta), flow.pg_v->nhb, flow.pg_v->nor);
        if (flow.therm > 20)
        {
            if (i % (flow.therm / 5) == 0)
                lprintf("MAIN", 0, "%d", ((i * 100) / flow.therm));
            else if (i % (flow.therm / 20) == 0)
                lprintf("MAIN", 0, ".");
        }
    }
    gettimeofday(&end, 0);
    timeval_subtract(&etime, &end, &start);
    if (i)
    {
        lprintf("MAIN", 0, "100\nThermalized %d Trajectories: [%ld sec %ld usec]\n", flow.therm, etime.tv_sec, etime.tv_usec);
        save_conf(&flow, GENERIC_MAX(0, flow.start - 1));
    }
    /* updates and measure  */
    i = flow.start - 1;

    while (i < flow.end)
    {
        for (n = 0; n < flow.nit; n++) /* nit updates */
        {
            i++;
            lprintf("MAIN", 0, "Trajectory #%d...\n", i);

            gettimeofday(&start, 0);

            update(&(flow.pg_v->beta), flow.pg_v->nhb, flow.pg_v->nor);

            gettimeofday(&end, 0);
            timeval_subtract(&etime, &end, &start);
            lprintf("MAIN", 0, "Trajectory #%d: generated in [%ld sec %ld usec]\n", i, etime.tv_sec, etime.tv_usec);
            lprintf("MAIN", 0, "Plaquette %1.18e\n", avr_plaquette());

            if ((i % flow.save_freq) == 0)
            {
                save_conf(&flow, i);
                if (rlx_var.rlxd_state[0] != '\0')
                {
                    lprintf("MAIN", 0, "Saving rlxd state to file %s\n", rlx_var.rlxd_state);
                    write_ranlxd_state(rlx_var.rlxd_state);
                }
            }
        }

        if (strcmp(flow.wf->make, "true") == 0)
        {
            static suNg_field *Vwf = NULL;
            if (Vwf == NULL)
                Vwf = alloc_gfield(&glattice);
            gettimeofday(&start, 0);
            suNg_field_copy(Vwf,u_gauge);
            WF_update_and_measure(flow.wf->ittype, Vwf, &(flow.wf->tmax), &(flow.wf->eps), &(flow.wf->delta), flow.wf->nmeas, DONTSTORE);
            gettimeofday(&end, 0);
            timeval_subtract(&etime, &end, &start);
            lprintf("MAIN", 0, "WF Measure #%d: generated in [%ld sec %ld usec]\n", i, etime.tv_sec, etime.tv_usec);
        }

        if (strcmp(flow.poly->make, "true") == 0)
        {
            gettimeofday(&start, 0);
            polyakov();
            gettimeofday(&end, 0);
            timeval_subtract(&etime, &end, &start);
            lprintf("MAIN", 0, "Polyakov Measure #%d: generated in [%ld sec %ld usec]\n", i, etime.tv_sec, etime.tv_usec);
        }
    }

    /* save final configuration */
    if ((i % flow.save_freq) != 0)
    {
        save_conf(&flow, i);
        /* Only save state if we have a file to save to */
        if (rlx_var.rlxd_state[0] != '\0')
        {
            lprintf("MAIN", 0, "Saving rlxd state to file %s\n", rlx_var.rlxd_state);
            write_ranlxd_state(rlx_var.rlxd_state);
        }
    }

    /* finalize Monte Carlo */
    end_mc();

    /* close communications */
    finalize_process();

    return 0;
}
