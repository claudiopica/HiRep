/*******************************************************************************
*
* Testing geometry
*
*******************************************************************************/

#define MAIN_PROGRAM

#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "global.h"
#include "io.h"
#include "geometry.h"
#include "logger.h"
#include "random.h"


int main(int argc,char *argv[])
{
  char tmp[256];
  
  setup_process(&argc,&argv);
  
  logger_setlevel(0,100); /* log all */
  if (PID!=0) { 
    logger_disable();}
  else{
    sprintf(tmp,">out_%d",PID); logger_stdout(tmp);
    sprintf(tmp,"err_%d",PID); freopen(tmp,"w",stderr);
  }
  
  logger_map("DEBUG","debug");
  
  lprintf("MAIN",0,"PId =  %d [world_size: %d]\n\n",PID,WORLD_SIZE); 
  
  read_input(glb_var.read,"test_input");
  
  /* setup communication geometry */
  if (geometry_init() == 1) {
    finalize_process();
    return 0;
  }
  
  geometry_mpi_eo();
  /* setup random numbers */
  read_input(rlx_var.read,"test_input");
  lprintf("MAIN",0,"RLXD [%d,%d]\n",rlx_var.rlxd_level,rlx_var.rlxd_seed+MPI_PID);
  rlxd_init(rlx_var.rlxd_level,rlx_var.rlxd_seed+MPI_PID); /* use unique MPI_PID to shift seeds */
  
  
  
  lprintf("MAIN",0,"Gauge group: SU(%d)\n",NG);
  lprintf("MAIN",0,"Fermion representation: dim = %d\n",NF);
  lprintf("MAIN",0,"The lattice size is %dx%dx%dx%d\n",T,X,Y,Z);
  lprintf("MAIN",0,"The lattice global size is %dx%dx%dx%d\n",GLB_T,GLB_X,GLB_Y,GLB_Z);
  lprintf("MAIN",0,"The lattice borders are (%d,%d,%d,%d)\n",T_BORDER,X_BORDER,Y_BORDER,Z_BORDER);
  lprintf("MAIN",0,"\n");
  fflush(stdout);

#ifdef WITH_MPI
  MPI_Barrier(GLB_COMM) ; 
#endif

  test_geometry_mpi_eo();
  
  finalize_process();

  return 0;
}
