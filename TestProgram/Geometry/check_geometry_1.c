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
#include "random.h"
#include <mpi.h>


int main(int argc,char *argv[])
{
  char tmp[512];
 setup_process(&argc,&argv);
  
  logger_setlevel(0,100); /* log all */
  if (PID!=0) { logger_disable(); }
  else{
    sprintf(tmp,"out_%d",PID); logger_stdout(tmp);
    sprintf(tmp,"err_%d",PID); freopen(tmp,"w",stderr);
  }
  logger_map("DEBUG","debug");
  
  lprintf("MAIN",0,"PId =  %d [world_size: %d]\n\n",PID,WORLD_SIZE); 
  
  read_input(glb_var.read,"test_input");
  rlxd_init(glb_var.rlxd_level,glb_var.rlxd_seed);
  
  /* setup communication geometry */
  if (geometry_init() == 1) {
    finalize_process();
    return 0;
  }

  geometry_mpi_eo();

  fflush(stdout);
  MPI_Barrier( GLB_COMM ) ; 

  test_geometry_mpi_eo();
  
  finalize_process();

  return 0;
}
