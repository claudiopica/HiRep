/*******************************************************************************
*
* Testing geometry
*
*******************************************************************************/

#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "io.h"
#include "geometry.h"
#include "global.h"
#include "logger.h"


int main(int argc,char *argv[])
{
  char tmp[256];

/*  GLB_T=16;*/
/*  GLB_X=4;*/
/*  GLB_Y=4;*/
/*  GLB_Z=4;*/

/*  NP_T=4;*/
/*  NP_X=1;*/
/*  NP_Y=1;*/
/*  NP_Z=1;*/

/*  T_BORDER=1;*/
/*  X_BORDER=0;*/
/*  Y_BORDER=0;*/
/*  Z_BORDER=0;*/

  setup_process(&argc,&argv);
  
  logger_setlevel(0,10000); /* log all */
  logger_map("DEBUG","debug");
#ifdef WITH_MPI
  sprintf(tmp,">out_%d",PID); logger_stdout(tmp);
  sprintf(tmp,"err_%d",PID); freopen(tmp,"w",stderr);
#endif

  lprintf("MAIN",0,"PId =  %d [world_size: %d]\n\n",PID,WORLD_SIZE); 

  read_input(glb_var.read,"test_input");

  /* setup communication geometry */
  if (geometry_init() == 1) {
    finalize_process();
    return 0;
  }

	geometry_mpi_eo();

	test_geometry_mpi_eo();
  
  finalize_process();
	exit(0);
}
