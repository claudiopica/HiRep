/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *   
* All rights reserved.                                                      * 
\***************************************************************************/

/*******************************************************************************
*
* Computation of the average plaquette
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


/* Input file variables */
typedef struct _input_pg {

  double beta;
  int nth, nms, nit, nhb, nor;

  /* for the reading function */
  input_record_t read[7];
  
} input_pg;

#define init_input_pg(varname) \
{ \
  .read={\
    {"beta", "beta = %lf", DOUBLE_T, &(varname).beta},\
    {"nth", "nth = %d", INT_T, &(varname).nth},\
    {"nms", "nms = %d", INT_T, &(varname).nms},\
    {"nit", "nit = %d", INT_T, &(varname).nit},\
    {"nhb", "nhb = %d", INT_T, &(varname).nhb},\
    {"nor", "nor = %d", INT_T, &(varname).nor},\
    {NULL, NULL, 0, NULL}\
  }\
}

input_pg pg_var = init_input_pg(pg_var);

int main(int argc,char *argv[])
{
  char tmp[256];
  int i,n,flag;
  double *p,pbar,sig,tau;

  /* setup process id and communications */
  setup_process(&argc,&argv);

  /* logger setup */
  logger_setlevel(0,40);
  sprintf(tmp,">out_%d",PID); logger_stdout(tmp);
  sprintf(tmp,"err_%d",PID); freopen(tmp,"w",stderr);

  lprintf("MAIN",0,"PId =  %d [world_size: %d]\n\n",PID,WORLD_SIZE); 

  /* read input file */
  read_input(glb_var.read,"common_input");
  lprintf("MAIN",0,"RLXD [%d,%d]\n",glb_var.rlxd_level,glb_var.rlxd_seed+PID);
  rlxd_init(glb_var.rlxd_level,glb_var.rlxd_seed+PID);
  read_input(pg_var.read,"input_file");

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

  /* disable logger for MPI processes != 0 */
  //if (PID!=0) { logger_disable(); }

  /* alloc global gauge fields */
  u_gauge=alloc_gfield(&glattice);
#ifndef REPR_FUNDAMENTAL
  u_gauge_f=alloc_gfield_f(&glattice);
#endif
  random_u(u_gauge); 
  represent_gauge_field();

  lprintf("MAIN",0,"beta = %2.4f\n",pg_var.beta);
  lprintf("MAIN",0,"nth  = %d\tNumber of thermalization cycles\n",pg_var.nth);
  lprintf("MAIN",0,"nms  = %d\tNumber of measure cycles\n",pg_var.nms);
  lprintf("MAIN",0,"nit  = %d\tNumber of hb-or iterations per cycle\n",pg_var.nit);
  lprintf("MAIN",0,"nhb  = %d\tNumber of heatbaths per iteration\n",pg_var.nhb);   
  lprintf("MAIN",0,"nor  = %d\tNumber or overrelaxations per iteration\n",pg_var.nor);

  /* random_u(); */

  p=malloc(pg_var.nms*sizeof(double)); /* array per le misure */

  /* Termalizzazione */
  for (i=0;i<pg_var.nth;++i) {
    update(pg_var.beta,pg_var.nhb,pg_var.nor);
    if ((i%10)==0)
      lprintf("MAIN",0,"%d",i);
    else
      lprintf("MAIN",0,".");
  }
  if(i) lprintf("MAIN",0,"%d\nThemalization done.\n",i);


  /* Misure */
  for (i=0;i<pg_var.nms;i++){ /* nms misure */
    p[i]=avr_plaquette();
    lprintf("MAIN",0,"[%d] <p> = %1.6f\n",i,p[i]);

    for (n=0;n<pg_var.nit;n++) /* nit updates */
      update(pg_var.beta,pg_var.nhb,pg_var.nor);

    /* Plaquette media */
  }

  lprintf("MAIN",0,"[%d] <p> = %1.6f\n",i,avr_plaquette());

  pbar=average(pg_var.nms,p);
  sig=sigma(pg_var.nms,p,&tau,&flag);

  if (flag!=0)
  {
    lprintf("MAIN",0,"Warning: unreliable error estimation ");
    lprintf("MAIN",0,"(data series is too short)\n\n");
  }

  lprintf("MAIN",0,"<p>   = %1.6f\n",pbar);
  lprintf("MAIN",0,"sigma = %1.6f\n",sig);
  lprintf("MAIN",0,"tau   = %1.2f [iterations]\n\n",tau);


  write_gauge_field("Config");

  /* free memory */
  free_gfield(u_gauge);
#ifndef REPR_FUNDAMENTAL
  free_gfield_f(u_gauge_f);
#endif
  /* close communications */
  finalize_process();

  return 0;
}
