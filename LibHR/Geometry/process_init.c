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
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

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
int setup_process(int *argc, char ***argv) {
  
#ifdef WITH_MPI
  /* MPI variables */
  int mpiret;
#endif
  
#ifdef WITH_MPI
  /* init MPI */
  mpiret=MPI_Init(argc,argv);
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
