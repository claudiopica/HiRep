/*******************************************************************************
*
* Test of modules
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
#include "global.h"
#include "suN.h"
#include "suN_types.h"
#include "linear_algebra.h"
#include "inverters.h"
#include "representation.h"
#include "utils.h"
#include "communications.h"
#include "logger.h"
#include "memory.h"


void test_wilson_action_and_force(double beta);
void test_lw_force(double beta, double c0, double c1);
void test_ginv_lw_action(double beta, double c0, double c1);
void test_gcov_lw_force(double beta, double c0, double c1);


int main(int argc,char *argv[])
{
   char pame[256];

   setup_process(&argc,&argv);
   
   logger_setlevel(0,100); /* log all */
	if (PID!=0) { logger_disable(); }   /* disable logger for MPI processes != 0 */
	else {
   //sprintf(pame,">out_%d",PID); logger_stdout(pame);
   sprintf(pame,"err_%d",PID); freopen(pame,"w",stderr);
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

   
   BCs_pars_t BCs_pars = {
     .fermion_twisting_theta = {0.,0.,0.,0.},
     .gauge_boundary_improvement_cs = 1.,
     .gauge_boundary_improvement_ct = 1.,
     .chiSF_boundary_improvement_ds = 1.,
     .SF_BCs = 0
   };
   init_BCs(&BCs_pars);

   lprintf("MAIN",0,"Gauge group: SU(%d)\n",NG);
   lprintf("MAIN",0,"global size is %dx%dx%dx%d\n",GLB_T,GLB_X,GLB_Y,GLB_Z);
   lprintf("MAIN",0,"proc grid is %dx%dx%dx%d\n",NP_T,NP_X,NP_Y,NP_Z);

   lprintf("MAIN",0,"local size is %dx%dx%dx%d\n",T,X,Y,Z);
   lprintf("MAIN",0,"extended local size is %dx%dx%dx%d\n",T_EXT,X_EXT,Y_EXT,Z_EXT);

   /* alloc global gauge fields */
   u_gauge=alloc_gfield(&glattice);
  
   
   test_wilson_action_and_force(1.);
   test_wilson_action_and_force(3.);
   
   test_ginv_lw_action(1.,1.,0.);
   test_ginv_lw_action(1.,0.,1.);
   
   test_gcov_lw_force(1.,1.,0.);
   test_gcov_lw_force(1.,0.,1.);

   test_lw_force(1.,1.,0.);
   test_lw_force(1.,0.,1.);
   

   free_gfield(u_gauge);

   finalize_process();

   exit(0);
}
