/******************************************************************************* 
* Check that the molecular dynamics evolution is reversible
* 
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
#include "communications.h"
#include "../RHMC/rhmc_utils.h"

int nhb,nor,nit,nth,nms,level,seed;
double beta;

extern suNg_av_field *momenta;

rhmc_flow flow=init_rhmc_flow(flow);

static void flip_mom()
{

  suNg_algebra_vector  *dptr;
 
  _DECLARE_INT_ITERATOR(ix);

  geometry_descriptor *gd=momenta->type;
  
  int dx;

   _MASTER_FOR(gd,ix) {
     for(dx=0;dx <4 ; dx++)
       {
	 dptr=(suNg_algebra_vector*)(&momenta->ptr[4*ix+dx]);

	_algebra_vector_mul_g(*dptr,-1.0,*dptr);
	
       }
   }
}




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
  lprintf("MAIN",0,"RLXD [%d,%d]\n",glb_var.rlxd_level,glb_var.rlxd_seed);


  rlxd_init(glb_var.rlxd_level,glb_var.rlxd_seed);
  
  
  
  
  /* setup communication geometry */
  if (geometry_init() == 1) {
    finalize_process();
    return 0;
  }
  
  geometry_mpi_eo();

  
  if(glb_var.rlxd_state[0]!='\0')
    {
      /*load saved state*/
      lprintf("MAIN",0,"Loading rlxd state from file %s\n",glb_var.rlxd_state);
      read_ranlxd_state(glb_var.rlxd_state);
    }
  else
    {
      lprintf("MAIN",0,"RLXD [%d,%d]\n",glb_var.rlxd_level,glb_var.rlxd_seed+PID);
      rlxd_init(glb_var.rlxd_level,glb_var.rlxd_seed+PID);
    }
  lprintf("MAIN",0,"Gauge group: SU(%d)\n",NG);
  lprintf("MAIN",0,"Fermion representation: " REPR_NAME " [dim=%d]\n",NF);
  

  /* test_geometry_mpi_eo(); */

  /* Init Monte Carlo */
  init_mc(&flow, "test_input");
  lprintf("MAIN",0,"MVM during RHMC initialzation: %ld\n",getMVM());
  lprintf("MAIN",0,"Initial plaquette: %1.8e\n",avr_plaquette());



  
  int rr=update_rhmc_o();

  if(rr<0) {
    lprintf("REV TEST",0,"Error in updating the gauge field!!\n");
    return 1;
  }
  lprintf("REV TEST",0,"Plaquette: %1.8e\n",avr_plaquette());
  
   flip_mom();

   rr=update_rhmc_o();
   if(rr<0) {
     lprintf("REV TEST",0,"Error in updating the gauge field!!\n");
     return 1;
   }
   lprintf("REV TEST",0,"Plaquette: %1.8e\n",avr_plaquette());
   
   free_rhmc();
   
   free_gfield(u_gauge);
#ifndef REPR_FUNDAMENTAL
   free_gfield_f(u_gauge_f);
#endif
   
   finalize_process();
   return 0;
}
