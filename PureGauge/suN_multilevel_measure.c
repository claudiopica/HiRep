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
#include "suN_utils_multilevel.h"
#include "setup.h"
#include "glueballs.h"
#include "wilsonflow.h"

pg_flow_ml_measure flow = init_pg_flow_ml_measure(flow);

int main(int argc, char *argv[])
{
  int i;
  FILE *list = NULL;
  char cnfg_filename[256];

  struct timeval start, end, etime; /* //for trajectory timing */

  setup_process(&argc, &argv);

  setup_gauge_fields();

  /* Init Monte Carlo */
  init_mc_ml_measure(&flow, get_input_filename());

  /* Measures */
  lprintf("MAIN", 0, "Configurations list from %s\n", flow.configlist);

  error((list = fopen(flow.configlist, "r")) == NULL, 1, "main [suN_multilevel_measure.c]",
        "Failed to open config list file\n");

  while (1)
  {
    if (fscanf(list, "%s", cnfg_filename) == 0 || feof(list))
      break;

    char *str;
    str = strrchr(cnfg_filename, 'n');
    error(sscanf(str, "n%d", &i) != 1, 1, "main [suN_multilevel_measure.c]",
          "Malformed configuration name (not ending by ...n<number>) \n");

    lprintf("MAIN", 0, "\n\nConfiguration %d from %s\n", i, cnfg_filename);

    read_gauge_field(cnfg_filename);

    apply_BCs_on_fundamental_gauge_field();

    gettimeofday(&start, 0);

    update_hb_multilevel_gb_measure(0, &(flow.pg_v->beta), flow.pg_v->nhb, flow.pg_v->nor, flow.pg_v->ml_niteration, flow.pg_v->ml_nskip, flow.pg_v->nblkstart, flow.pg_v->nblkend, &(flow.pg_v->APEsmear), &(flow.pg_v->corrs));

    gettimeofday(&end, 0);
    timeval_subtract(&etime, &end, &start);
    lprintf("MAIN", 0, "ML Measure & update#%d: generated in [%ld sec %ld usec]\n", i, etime.tv_sec, etime.tv_usec);
    lprintf("MAIN", 0, "Plaquette %1.18e\n", avr_plaquette());

    if (strcmp(flow.wf->make, "true") == 0)
    {
      static suNg_field *Vwf = NULL;
      if (Vwf == NULL)
        Vwf = alloc_gfield(&glattice);
      gettimeofday(&start, 0);
      suNg_field_copy(Vwf, u_gauge);
      WF_update_and_measure(RK3_ADAPTIVE, Vwf, &(flow.wf->tmax), &(flow.wf->eps), &(flow.wf->delta), flow.wf->nmeas, DONTSTORE);
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

  /* close communications */
  finalize_process();

  return 0;
}
