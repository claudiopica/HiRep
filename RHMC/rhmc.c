/***************************************************************************\
 * Copyright (c) 2008, Claudio Pica                                          *   
 * All rights reserved.                                                      * 
 \***************************************************************************/

/*******************************************************************************
 *
 * Main RHMC program
 *
 *******************************************************************************/

#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
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
#include "rhmc_utils.h"

/* flow control variable */
rhmc_flow flow=init_rhmc_flow(flow);

int main(int argc,char *argv[])
{
  int i, acc, rc;
  char sbuf[128];

  /* setup process id and communications */
  setup_process(&argc,&argv);

  /* logger setup */
  logger_setlevel(0,10);
  /* disable logger for MPI processes != 0 */
  if (PID!=0) { logger_disable(); }

  if (PID==0) {
    sprintf(sbuf,">>out_%d",PID); logger_stdout(sbuf);
    sprintf(sbuf,"err_%d",PID); freopen(sbuf,"w",stderr);
  }

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
  lprintf("MAIN",0,"Fermion representation: " REPR_NAME " [dim=%d]\n",NF);
  lprintf("MAIN",0,"global size is %dx%dx%dx%d\n",GLB_T,GLB_X,GLB_Y,GLB_Z);
  lprintf("MAIN",0,"proc grid is %dx%dx%dx%d\n",NP_T,NP_X,NP_Y,NP_Z);

  /* setup lattice geometry */
  geometry_mpi_eo();
  /* test_geometry_mpi_eo(); */

  /* Init Monte Carlo */
  init_mc(&flow, "input_file");
  lprintf("MAIN",0,"MVM during RHMC initialzation: %ld\n",getMVM());
  lprintf("MAIN",0,"Initial plaquette: %1.8e\n",avr_plaquette());

  rc=acc=0;
  for(i=flow.start;i<flow.end;++i) {
    int rr;
    double perc;
    lprintf("MAIN",0,"Trajectory #%d...\n",i);
    rr=update_rhmc();
    if(rr<0) {
      lprintf("MAIN",0,"Error in updating the gauge field!!\n");
      return 1;
    } else if(rr!=0) {
      acc++;
    }
    rc++;
    perc=(acc==0)?0.:(float)(100*acc)/(float)(rc);

    lprintf("MAIN",0,"Trajectory #%d: %d/%d (%3.4f%%) MVM = %ld\n",i,acc,rc,perc,getMVM());

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
