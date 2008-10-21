/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *   
* All rights reserved.                                                      * 
\***************************************************************************/

#include "io.h"
#include "error.h"
#include "global.h"
#include "logger.h"
#include <stdio.h>

static void mpi_broadcast_parameters() {
#ifdef WITH_MPI
  MPI_Bcast(&input_p.GLB_T,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&input_p.GLB_X,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&input_p.GLB_Y,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&input_p.GLB_Z,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&input_p.NP_T,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&input_p.NP_X,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&input_p.NP_Y,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&input_p.NP_Z,1,MPI_INT,0,MPI_COMM_WORLD);

  MPI_Bcast(&input_p.rlxd_level,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&input_p.rlxd_seed,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&input_p.rhmc_p.beta,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
  MPI_Bcast(&input_p.rhmc_p.nf,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&input_p.rhmc_p.mass,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
  MPI_Bcast(&input_p.rhmc_p.MT_prec,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
  MPI_Bcast(&input_p.rhmc_p.MD_prec,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
  MPI_Bcast(&input_p.rhmc_p.HB_prec,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
  MPI_Bcast(&input_p.rhmc_p.force_prec,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
  MPI_Bcast(&input_p.rhmc_p.n_pf,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&input_p.int_p.tlen,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
  MPI_Bcast(&input_p.int_p.nsteps,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&input_p.int_p.gsteps,1,MPI_INT,0,MPI_COMM_WORLD);

  /*
    MPI_Bcast(&nhb,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&nor,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&nit,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&nth,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&nms,1,MPI_INT,0,MPI_COMM_WORLD);
  */
#endif

  GLB_T=input_p.GLB_T;
  GLB_X=input_p.GLB_X;
  GLB_Y=input_p.GLB_Y;
  GLB_Z=input_p.GLB_Z;
  NP_T=input_p.NP_T;
  NP_X=input_p.NP_X;
  NP_Y=input_p.NP_Y;
  NP_Z=input_p.NP_Z;

}

void read_input(char *filename) {
  FILE *inputfile;
  fpos_t pos;
  char buf[256];
  int count;

  /* when using mpi only PID=0 reads the input file */
  if (PID!=0) {
    mpi_broadcast_parameters();
    return;
  }

  error((inputfile=fopen(filename,"r"))==NULL,1,"read_input " __FILE__ ,
        "Failed to open input file\n");

  do {
    /*skip white space*/
    fscanf(inputfile," ");

    fgetpos(inputfile,&pos);

    /* skip comments */
    count=0;
    fscanf(inputfile,"//%n",&count);
    if(count==2) { goto nextline; } if(feof(inputfile)) break; fsetpos(inputfile,&pos);

    /* read global variables */
#define READVAR(NAME,VAR)\
    if(fscanf(inputfile, NAME, &VAR)==1) { continue;} if(feof(inputfile)) break; fsetpos(inputfile,&pos)

    /* geometry related */
    READVAR("GLB_T = %d",input_p.GLB_T);
    READVAR("GLB_X = %d",input_p.GLB_X);
    READVAR("GLB_Y = %d",input_p.GLB_Y);
    READVAR("GLB_Z = %d",input_p.GLB_Z);
    READVAR("NP_T = %d",input_p.NP_T);
    READVAR("NP_X = %d",input_p.NP_X);
    READVAR("NP_Y = %d",input_p.NP_Y);
    READVAR("NP_Z = %d",input_p.NP_Z);

    /* random numbers */
    READVAR("level = %d",input_p.rlxd_level);
    READVAR("seed = %d",input_p.rlxd_seed);

    /* simulation parameters */
    READVAR("beta = %lf",input_p.rhmc_p.beta);
    READVAR("nf = %d",input_p.rhmc_p.nf);
    READVAR("mass = %lf",input_p.rhmc_p.mass);
    READVAR("MT_prec = %lf",input_p.rhmc_p.MT_prec);
    READVAR("MD_prec = %lf",input_p.rhmc_p.MD_prec);
    READVAR("HB_prec = %lf",input_p.rhmc_p.HB_prec);
    READVAR("force_prec = %lf",input_p.rhmc_p.force_prec);
    READVAR("n_pf = %u",input_p.rhmc_p.n_pf);
    READVAR("tlen = %lf",input_p.int_p.tlen);
    READVAR("nsteps = %u",input_p.int_p.nsteps);
    READVAR("gsteps = %u",input_p.int_p.gsteps);
   

#undef READVAR
    
    if(ferror(inputfile)) {
      lprintf("READINPUT",0,"Cannot read input file.\n");
      perror(0);
    }

    fscanf(inputfile,"%s",buf);
    lprintf("READINPUT",0,"Unknown token: [%s]\n",buf);
	
    error(1,1,"read_input " __FILE__, "Unknown directive in input file\n");

  nextline:
    fscanf(inputfile,"%*[^\n]");
  } while (!feof(inputfile) && !ferror(inputfile));

  fclose(inputfile);

  mpi_broadcast_parameters();

}




