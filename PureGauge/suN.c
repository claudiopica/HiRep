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
#include <string.h>
#include <math.h>
#include "io.h"
#include "random.h"
#include "error.h"
#include "geometry.h"
#include "memory.h"
#include "statistics.h"
#include "update.h"
#include "global.h"
#include "logger.h"
#include "observables.h"
#include "representation.h"
#include "suN_utils.h"

pg_flow flow=init_pg_flow(flow);

int main(int argc,char *argv[])
{
  char tmp[256];
  int i,n;

  /* setup process id and communications */
  setup_process(&argc,&argv);

  /* logger setup */
  logger_setlevel(0,40);
  sprintf(tmp,">out_%d",PID); logger_stdout(tmp);
  sprintf(tmp,"err_%d",PID); freopen(tmp,"w",stderr);

  lprintf("MAIN",0,"PId =  %d [world_size: %d]\n\n",PID,WORLD_SIZE); 

  /* read input file */
  read_input(glb_var.read,"input_file");
  lprintf("MAIN",0,"RLXD [%d,%d]\n",glb_var.rlxd_level,glb_var.rlxd_seed+PID);
  rlxd_init(glb_var.rlxd_level,glb_var.rlxd_seed+PID);

  /* setup communication geometry */
  if (geometry_init() == 1) {
    finalize_process();
    return 0;
  }

  lprintf("MAIN",0,"Gauge group: SU(%d)\n",NG);
  lprintf("MAIN",1,"Fermion representation: " REPR_NAME " [dim=%d]\n",NF);
  lprintf("MAIN",0,"global size is %dx%dx%dx%d\n",GLB_T,GLB_X,GLB_Y,GLB_Z);
  lprintf("MAIN",0,"proc grid is %dx%dx%dx%d\n",NP_T,NP_X,NP_Y,NP_Z);

  /* setup lattice geometry */
  geometry_mpi_eo();
  /* test_geometry_mpi_eo(); */

  /* Init Monte Carlo */
  init_mc(&flow, "input_file");
  lprintf("MAIN",0,"Initial plaquette: %1.8e\n",avr_plaquette());


  /* Termalizzazione */
  for (i=0;i<flow.therm;++i) {
    update(flow.pg_v->beta,flow.pg_v->nhb,flow.pg_v->nor);
    if ((i%10)==0)
      lprintf("MAIN",0,"%d",i);
    else
      lprintf("MAIN",0,".");
  }
  if(i) lprintf("MAIN",0,"%d\nThemalization done.\n",i);


  /* Misure */
  for(i=flow.start;i<flow.end;++i) {
    lprintf("MAIN",0,"Trajectory #%d...\n",i);

    for (n=0;n<flow.pg_v->nit;n++) /* nit updates */
      update(flow.pg_v->beta,flow.pg_v->nhb,flow.pg_v->nor);

    if((i%flow.save_freq)==0) {
      save_conf(&flow, i);
    }

    if((i%flow.meas_freq)==0) {
      lprintf("MAIN",0,"Plaquette: %1.8e\n",avr_plaquette());
      /* do something */
    }
  }


  /* finalize Monte Carlo */
  end_mc();

  /* close communications */
  finalize_process();

  return 0;
}
