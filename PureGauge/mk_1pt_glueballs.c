/***************************************************************************\
* Copyright (c) 2025, Antonio Rago
* All rights reserved.                                                      * 
\***************************************************************************/

/*******************************************************************************
*
* Main program to measure glueballs with pure gauge interpolators 
*
*******************************************************************************/

#include "libhr.h"
#include "suN_utils.h"
#include <unistd.h>

pg_flow_glueballs_measure flow = init_pg_flow_glueballs_measure(flow);

int main(int argc, char *argv[]) {
    FILE *list = NULL;
    char cnfg_filename[256];

    struct timeval start, end, etime; /* //for measurment timing */

    setup_process(&argc, &argv);

    setup_gauge_fields();

    /* Init Monte Carlo */
    init_mk_glueballs(&flow, get_input_filename());

    int nblocking = flow.pg_v->nblkend - flow.pg_v->nblkstart + 1;

    hr_complex *one_point_gb = (hr_complex *)amalloc(sizeof(hr_complex) * total_n_glue_op * nblocking * n_active_slices, ALIGN);
    error(one_point_gb == NULL, 1, __func__, "Could not allocate memory space for field structure");
    hr_complex *one_point_tor = (hr_complex *)amalloc(sizeof(hr_complex) * total_n_tor_op * n_active_slices, ALIGN);
    error(one_point_tor == NULL, 1, __func__, "Could not allocate memory space for field structure");
    hr_complex **polyf = (hr_complex **)amalloc(sizeof(hr_complex *) * 3, ALIGN);
    error(polyf == NULL, 1, __func__, "Could not allocate memory space for field structure");
    polyf[0] = (hr_complex *)amalloc(sizeof(hr_complex) * Y * Z * T, ALIGN);
    error(polyf[0] == NULL, 1, __func__, "Could not allocate memory space for field structure");
    polyf[1] = (hr_complex *)amalloc(sizeof(hr_complex) * X * Z * T, ALIGN);
    error(polyf[1] == NULL, 1, __func__, "Could not allocate memory space for field structure");
    polyf[2] = (hr_complex *)amalloc(sizeof(hr_complex) * X * Y * T, ALIGN);
    error(polyf[2] == NULL, 1, __func__, "Could not allocate memory space for field structure");

    /* Measures */
    lprintf("MAIN", 0, "Configurations list from %s\n", flow.configlist);

    error((list = fopen(flow.configlist, "r")) == NULL, 1, "main [suN_multilevel_measure.c]",
          "Failed to open config list file\n");

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

#if total_n_glue_op > 0
        for (int i = 0; i < total_n_glue_op * nblocking * n_active_slices; i++) {
            one_point_gb[i] = 0.;
        }
        measure_1pt_glueballs(flow.pg_v->nblkstart, flow.pg_v->nblkend, &(flow.pg_v->APEsmear), one_point_gb);
        for (int j = 0; j < n_active_slices; j++) {
            lprintf("step", 0, "%d ", zerocoord[0] + j);
            for (int i = 0; i < total_n_glue_op * nblocking; i++) {
                lprintf("step", 0, " (%.10e +I*(%.10e)) ", creal(one_point_gb[i]), cimag(one_point_gb[i]));
            }
            lprintf("step", 0, "\n");
        }
#endif

#if total_n_tor_op > 0
        for (int i = 0; i < total_n_tor_op * n_active_slices; i++) {
            one_point_tor[i] = 0.;
        }
        for (int i = 0; i < Y * Z * T; i++) {
            polyf[0][i] = 0.;
        }
        for (int i = 0; i < X * Z * T; i++) {
            polyf[1][i] = 0.;
        }
        for (int i = 0; i < X * Y * T; i++) {
            polyf[2][i] = 0.;
        }
        measure_1pt_torellons(&(flow.pg_v->APEsmear), one_point_tor, polyf);
#endif

        gettimeofday(&end, 0);
        timeval_subtract(&etime, &end, &start);
        lprintf("MAIN", 0, "Glueballs & Torellons 1pt #%d: generated in [%ld sec %ld usec]\n", confid, etime.tv_sec,
                etime.tv_usec);
        lprintf("MAIN", 0, "Plaquette %1.18e\n", avr_plaquette());

        if (strcmp(flow.wf->make, "true") == 0) {
            static suNg_field *Vwf = NULL;
            if (Vwf == NULL) { Vwf = alloc_suNg_field(&glattice); }
            gettimeofday(&start, 0);
            copy_suNg_field(Vwf, u_gauge);
            WF_update_and_measure(RK3_ADAPTIVE, Vwf, &(flow.wf->tmax), &(flow.wf->eps), &(flow.wf->delta), flow.wf->nmeas,
                                  DONTSTORE);
            gettimeofday(&end, 0);
            timeval_subtract(&etime, &end, &start);
            lprintf("MAIN", 0, "WF Measure #%d: generated in [%ld sec %ld usec]\n", confid, etime.tv_sec, etime.tv_usec);
        }

        if (strcmp(flow.poly->make, "true") == 0) {
            gettimeofday(&start, 0);
            polyakov();
            gettimeofday(&end, 0);
            timeval_subtract(&etime, &end, &start);
            lprintf("MAIN", 0, "Polyakov Measure #%d: generated in [%ld sec %ld usec]\n", confid, etime.tv_sec, etime.tv_usec);
        }
    }
    afree(one_point_gb);
    afree(one_point_tor);

    afree(polyf[0]);
    afree(polyf[1]);
    afree(polyf[2]);
    afree(polyf);

    /* close communications */
    finalize_process();

    return 0;
}
