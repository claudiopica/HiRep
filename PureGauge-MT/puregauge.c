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
#include "puregauge_utils.h"
#include "memory.h"
#include "communications.h"
#include "observables.h"
#include "utils.h"

#include "cinfo.c"


/* Polyakov-loop parameters */
typedef struct _input_polyakov {
  char make[256];

  /* for the reading function */
  input_record_t read[2];

} input_polyakov;

#define init_input_polyakov(varname)					\
  {									\
    .read={								\
      {"make polyakov loops", "poly:make = %s", STRING_T, (varname).make}, \
      {NULL, NULL, INT_T, NULL}						\
    }									\
  }


input_polyakov poly_var = init_input_polyakov(poly_var);



/* flow control variable */
puregauge_flow flow=init_puregauge_flow(flow);

char input_filename[256] = "input_file";
char output_filename[256] = "out_0";
char error_filename[256] = "err_0";
static void read_cmdline(int argc, char* argv[]) {
  int i, ai=0, ao=0, am=0, requested=1;
  
  for (i=1;i<argc;i++) {
    if (strcmp(argv[i],"-i")==0) {ai=i+1;requested+=2;}
    else if (strcmp(argv[i],"-o")==0) {ao=i+1;requested+=2;}
    else if (strcmp(argv[i],"-m")==0) {am=i;requested+=1;}
  }
  
  if (am != 0) {
    print_compiling_info();
    exit(0);
  }
  
  error(argc!=requested,1,"read_cmdline [puregauge.c]",
        "Arguments: [-i <input file>] [-o <output file>] [-m]");
  
  if (ao!=0) strcpy(output_filename,argv[ao]);
  if (ai!=0) strcpy(input_filename,argv[ai]);
}



int main(int argc,char *argv[])
{
  int i, acc, rc;
  char sbuf[128];
  
  read_cmdline(argc,argv);
  
  /* setup process communications */
  setup_process(&argc,&argv);

  /* read global variables file */
  read_input(glb_var.read,input_filename);

  setup_replicas();

  /* logger setup */
  read_input(logger_var.read,input_filename);
  logger_set_input(&logger_var);
  if (PID!=0) { logger_disable(); }   /* disable logger for MPI processes != 0 */
  else {
    sprintf(sbuf,">>%s",output_filename);  logger_stdout(sbuf);
    freopen(error_filename,"w",stderr);
  }

  lprintf("MAIN",0,"Compiled with macros: %s\n",MACROS); 
  lprintf("MAIN",0,"[RepID: %d][world_size: %d]\n[MPI_ID: %d][MPI_size: %d]\n\n",RID,WORLD_SIZE,MPI_PID,MPI_WORLD_SIZE);

  /* setup lattice geometry */
  if (geometry_init() == 1) { finalize_process(); return 0; }
  geometry_mpi_eo();
  /* test_geometry_mpi_eo(); */
    /* setup random numbers */
    read_input(rlx_var.read,input_filename);
    //slower(rlx_var.rlxd_start); //convert start variable to lowercase
    if(strcmp(rlx_var.rlxd_start,"continue")==0 && rlx_var.rlxd_state[0]!='\0') {
        /*load saved state*/
        lprintf("MAIN",0,"Loading rlxd state from file [%s]\n",rlx_var.rlxd_state);
        read_ranlxd_state(rlx_var.rlxd_state);
    } else {
        lprintf("MAIN",0,"RLXD [%d,%d]\n",rlx_var.rlxd_level,rlx_var.rlxd_seed+MPI_PID);
        rlxd_init(rlx_var.rlxd_level,rlx_var.rlxd_seed+MPI_PID); /* use unique MPI_PID to shift seeds */
    }
  
  lprintf("MAIN",0,"Gauge group: SU(%d)\n",NG);
  lprintf("MAIN",0,"Fermion representation: " REPR_NAME " [dim=%d]\n",NF);
  

  /* read input for measures */
  read_input(poly_var.read,input_filename);
  
  /* Init Monte Carlo */
  init_mc(&flow, input_filename);
  lprintf("MAIN",0,"Initial plaquette: %1.8e\n",avr_plaquette());


  
  rc=acc=0;
  for(i=flow.start;i<flow.end;++i) {
    int rr;
    double perc;
    struct timeval start, end, etime; /* //for trajectory timing */
    lprintf("BLOCK",0," Start %d\n",i);
    lprintf("MAIN",0,"Trajectory #%d...\n",i);
    
    gettimeofday(&start,0);
    
    rr=update_puregauge();

    gettimeofday(&end,0);
    timeval_subtract(&etime,&end,&start);
    lprintf("MAIN",0,"Trajectory #%d: generated in [%ld sec %ld usec]\n",i,etime.tv_sec,etime.tv_usec);
      
    if(rr<0) {
      lprintf("MAIN",0,"Error in updating the gauge field!!\n");
      return 1;
    } else if(rr!=0) {
      acc++;
    }
    rc++;
    perc=(acc==0)?0.:(float)(100*acc)/(float)(rc);

    lprintf("MAIN",0,"Trajectory #%d: %d/%d (%3.4f%%)\n",i,acc,rc,perc);

    if((i%flow.save_freq)==0) {
      save_conf(&flow, i);
      /* Only save state if we have a file to save to */
      if(rlx_var.rlxd_state[0]!='\0') {
        lprintf("MAIN",0,"Saving rlxd state to file %s\n",rlx_var.rlxd_state);
        write_ranlxd_state(rlx_var.rlxd_state);
      }
    }

    if((i%flow.meas_freq)==0) {
      /* plaquette */
      lprintf("MAIN",0,"Plaquette: %1.8e\n",avr_plaquette());

      /* Polyakov loops */
      if(strcmp(poly_var.make,"true")==0) {
        polyakov();
      }
    }

    lprintf("BLOCK",0," End %d\n",i);

  }
  /* save final configuration */
  if(((--i)%flow.save_freq)!=0) {
    save_conf(&flow, i);
    /* Only save state if we have a file to save to */
    if(rlx_var.rlxd_state[0]!='\0') {
      lprintf("MAIN",0,"Saving rlxd state to file %s\n",rlx_var.rlxd_state);
      write_ranlxd_state(rlx_var.rlxd_state);
    }
  }

  
  /* finalize Monte Carlo */
  end_mc();

  /* close communications */
  finalize_process();

  return 0;

}
