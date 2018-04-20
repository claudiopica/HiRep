/***************************************************************************\
 * Copyright (c) 2008, Claudio Pica                                          *
 * All rights reserved.                                                      *
 \***************************************************************************/

/*******************************************************************************
 *
 * File process_init.c
 *
 * Inizialization of geometry structures
 *
 *******************************************************************************/

#include "global.h"
#include "error.h"
#include "logger.h"
#include "hr_omp.h"
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include "utils.h"
#include "io.h"

/* setup_process
 * Assign a unique RID, PID to each process and setup
 * communications as necessary
 *
 * return codes:
 * 0 => success
 *
 * OUTPUT:
 * MPI: GLB_COMM, RID, PID, WORLD_SIZE
 */

static char input_filename[256]="input_file";
static char output_filename[256]="out_0";
static char error_filename[256]="err_0";

char* get_input_filename(){return input_filename;}
char* get_output_filename(){return output_filename;}
char* get_error_filename(){return error_filename;}


static void read_cmdline(int argc, char** argv) {
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

  error(argc!=requested,1,"read_cmdline [hmc.c]",
      "Arguments: [-i <input file>] [-o <output file>] [-m]");

  if (ao!=0) strcpy(output_filename,argv[ao]);
  if (ai!=0) strcpy(input_filename,argv[ai]);

}



int setup_process(int *argc, char ***argv) {
  
  read_cmdline(*argc,*argv) ;
  
  
#ifdef WITH_MPI
  /* INIT MPI*/
  int mpiret;
  int required=MPI_THREAD_SINGLE;
  int provided;
  mpiret=MPI_Init_thread(argc,argv,required,&provided);
  if (mpiret != MPI_SUCCESS) {
    char mesg[MPI_MAX_ERROR_STRING];
    int mesglen;
    MPI_Error_string(mpiret,mesg,&mesglen);
    lprintf("MPI",0,"ERROR: %s\n",mesg);
    error(1,1,"setup_process " __FILE__,"MPI inizialization failed");
  }
  MPI_Comm_rank(MPI_COMM_WORLD, &MPI_PID);
  MPI_Comm_size(MPI_COMM_WORLD, &MPI_WORLD_SIZE);
  PID = MPI_PID;
  WORLD_SIZE = MPI_WORLD_SIZE;
  RID=0;
  
#else
  RID=0;
  PID=0;
  WORLD_SIZE=1;
#endif
  

 


  /* read global variables file */
  read_input(glb_var.read,input_filename);
  
  setup_replicas();
  
  /* logger setup */
  read_input(logger_var.read,input_filename);
  logger_set_input(&logger_var);
  if (PID!=0) { logger_disable(); }   /* disable logger for MPI processes != 0 */
  else {
    FILE* stderrp;
    char sbuf[256];
    sprintf(sbuf,">>%s",output_filename);  
    logger_stdout(sbuf);
    stderrp=freopen(error_filename,"w",stderr);
    error(stderrp==NULL,1,"main [hmc.c]",
	  "Cannot redirect the stderr");
  }
  

#ifdef _OPENMP
#pragma omp parallel 
  {
#pragma omp master
    {
      lprintf ("OMP", 0, "Number of Threads requested = %i\n", omp_get_num_threads());
    }
  }
#endif
  
  
  lprintf("SYSTEM",0,"Gauge group: SU(%d)\n",NG);                           
  lprintf("SYSTEM",0,"Fermion representation: dim = %d\n",NF);
  lprintf("SYSTEM",0,"[RepID: %d][world_size: %d]\n[MPI_ID: %d][MPI_size: %d]\n",RID,WORLD_SIZE,MPI_PID,MPI_WORLD_SIZE);
  print_compiling_info_short();

  //  lprintf("MAIN",0,"Logger lelvel: %d\n",logger_getlevel(0));
  
  /* setup lattice geometry */
  if (geometry_init() == 1) { finalize_process(); return 0; }


  return 0;
  
}

/* this function is intended to clean up before process ending
 *
 * return codes:
 * 0 => success
 */
int finalize_process() {
#ifdef WITH_MPI
  /* MPI variables */
  int init;
#endif
  
#ifdef WITH_MPI
  MPI_Initialized(&init);
  if(init) MPI_Finalize();
#endif
  
  return 0;
  
}

/* setup_replicas
 * Split MPI_COMM_WORLD into replica communicators GLB_COMM
 *
 * return codes:
 * 0 => success
 *
 * AFFECTS THE GLOBAL VARIABLES: GLB_COMM, RID, PID, WORLD_SIZE
 */
int setup_replicas() {
#ifdef WITH_MPI
  
  if (N_REP>1) {
    int mpiret;
    
    MPI_Initialized(&mpiret);
    if(!mpiret) {
      lprintf("MPI",0,"ERROR: MPI has not been initialized!!!\n");
      error(1,1,"setup_replicas " __FILE__,"Cannot create replicas");
    }
    
    /* set up replicas */
    char sbuf[64];
    if ((MPI_WORLD_SIZE%N_REP)!=0) {
      error(1,1,"setup_replicas " __FILE__,"MPI_WORLD_SIZE is not a multiple of the number of replicas!");
    }
    
    RID = MPI_PID/(MPI_WORLD_SIZE/N_REP); /* Replica ID */
    mpiret=MPI_Comm_split(MPI_COMM_WORLD, RID, 0, &GLB_COMM);
    if (mpiret != MPI_SUCCESS) {
      char mesg[MPI_MAX_ERROR_STRING];
      int mesglen;
      MPI_Error_string(mpiret,mesg,&mesglen);
      lprintf("MPI",0,"ERROR: %s\n",mesg);
      error(1,1,"setup_replicas " __FILE__,"Inizialization of replicas failed");
    }
    MPI_Comm_rank(GLB_COMM, &PID);
    MPI_Comm_size(GLB_COMM, &WORLD_SIZE);
    
    /* chdir to replica dir */
    sprintf(sbuf,"Rep_%d",RID);
    mpiret = chdir(sbuf);
  }
  
#endif //ifdef WITH_MPI
  
  return 0;
}
