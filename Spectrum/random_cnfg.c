/*******************************************************************************
*
* Generate random cnfg
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
#include "observables.h"
#include "suN.h"
#include "suN_types.h"
#include "dirac.h"
#include "linear_algebra.h"
#include "inverters.h"
#include "representation.h"
#include "utils.h"
#include "logger.h"


char input_filename[256] = "input_file";
char output_filename[256] = "random_cnfg.out";

input_glb glb_ip = init_input_glb(glb_ip);


void read_cmdline(int argc, char* argv[]) {
  int i, ai=0, ao=0;
  
  for (i=1;i<argc;i++) {
    if (strcmp(argv[i],"-i")==0) ai=i+1;
    else if (strcmp(argv[i],"-o")==0) ao=i+1;
  }

  if (ao!=0) strcpy(output_filename,argv[ao]);
  if (ai!=0) strcpy(input_filename,argv[ai]);
}


int main(int argc,char *argv[]) {
  char cnfg_filename[256];
  char tmp[256];

  /* setup process id and communications */
  setup_process(&argc,&argv);


  /* logger setup */
  /* disable logger for MPI processes != 0 */
  if (PID!=0) { logger_disable(); }
  logger_setlevel(0,40);
  sprintf(tmp,">%s",output_filename); logger_stdout(tmp);
  sprintf(tmp,"err_%d",PID); freopen(tmp,"w",stderr);

  lprintf("MAIN",0,"PId =  %d [world_size: %d]\n\n",PID,WORLD_SIZE); 


  /* read & broadcast parameters */
  read_cmdline(argc, argv);
  read_input(glb_ip.read,input_filename);


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

  lprintf("MAIN",0,"local size is %dx%dx%dx%d\n",T,X,Y,Z);
  lprintf("MAIN",0,"extended local size is %dx%dx%dx%d\n",T_EXT,X_EXT,Y_EXT,Z_EXT);

  lprintf("MAIN",0,"RLXD [%d,%d]\n",glb_ip.rlxd_level,glb_ip.rlxd_seed);
  rlxd_init(glb_ip.rlxd_level,glb_ip.rlxd_seed+PID);
 
  /* alloc global gauge fields */
  u_gauge=alloc_gfield(&glattice);
#ifndef REPR_FUNDAMENTAL
  u_gauge_f=alloc_gfield_f(&glattice);
#endif
  random_u(u_gauge);

  sprintf(cnfg_filename,"%dx%dx%dx%dNc%d",GLB_T,GLB_X,GLB_Y,GLB_Z,NG);
/*  write_gauge_field_eo_lexi(cnfg_filename);*/
  write_gauge_field(cnfg_filename);

	free_gfield(u_gauge);
#ifndef REPR_FUNDAMENTAL
	free_gfield_f(u_gauge_f);
#endif

	return 0;
}

