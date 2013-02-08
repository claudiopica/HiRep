/*******************************************************************************
*
*
*******************************************************************************/

#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "io.h"
#include "update.h"
#include "observables.h"
#include "geometry.h"
#include "global.h"
#include "logger.h"
#include "random.h"
#include "memory.h"
#include "dirac.h"
#include "linear_algebra.h"
#include "representation.h"
#include "communications.h"
#include "random.h"

rhmc_par _update_par;
/* double M_PI=3.141592653589793238462643383279502884197; */


int main(int argc,char *argv[])
{
  setup_process(&argc,&argv);
  char tmp[256];

  BCs_pars_t BCs_pars = {
    .fermion_twisting_theta = {0.,M_PI/5.,M_PI/5.,M_PI/5.},

/*      .fermion_twisting_theta={0.,0.,0.,0.},  */

    .gauge_boundary_improvement_cs = 1.,
    .gauge_boundary_improvement_ct = 1.,
    .chiSF_boundary_improvement_ds = 0.5,
    .SF_BCs = 0
  };
  

  double acc=1.e-20;


   double mass=-0.3099358667; 
/*  double mass=0.0; */

  _update_par.mass=-0.3099358667;
  _update_par.SF_ds=0.5;
  _update_par.SF_sign=1;
  _update_par.SF_ct=1;
  _update_par.SF_zf=1.3;

  
  logger_setlevel(0,100); /* log all */
  if (PID!=0) { 
    logger_disable();}
  else{
    sprintf(tmp,">out_%d",PID); logger_stdout(tmp);
    sprintf(tmp,"err_%d",PID); freopen(tmp,"w",stderr);
  }
  
  lprintf("MAIN",0,"PId =  %d [world_size: %d]\n\n",PID,WORLD_SIZE); 
  
  read_input(glb_var.read,"test_input");
  lprintf("MAIN",0,"RLXD [%d,%d]\n",glb_var.rlxd_level,glb_var.rlxd_seed);


  rlxd_init(glb_var.rlxd_level,glb_var.rlxd_seed);
  
#if NG!=3 || NF!=3
#error "Can work only with NC=3 and Nf==3"
#endif
  
  
  /* setup communication geometry */
  if (geometry_init() == 1) {
    finalize_process();
    return 0;
  }
  
  geometry_mpi_eo();

  init_BCs(&BCs_pars);




  lprintf("MAIN",0,"This test implements a comparison with a working code of Stefan Sint\n");
  
  u_gauge=alloc_gfield(&glattice);
  u_gauge_f=alloc_gfield_f(&glattice); 
  
  /* unit_u(u_gauge); */
  read_gauge_field_nocheck("gfield_sint.dat"); 
  apply_BCs_on_fundamental_gauge_field();
   
  lprintf("MAIN",0,"mass = %f\n",mass);
/*   lprintf("MAIN",0,"ds = %f\n",_update_par.SF_ds); */
/*   lprintf("MAIN",0,"zf = %f\n",_update_par.SF_zf); */
/*   lprintf("MAIN",0,"theta = %f\n",_update_par.SF_theta); */
/*   lprintf("MAIN",0,"sign = %d\n",_update_par.SF_sign); */
  lprintf("MAIN",0,"To be compared with results in file gxly_testconfig\n");

  represent_gauge_field(); 
  full_plaquette();
 
  SF_PCAC_wall_mass(mass,acc); 
  
  free_gfield(u_gauge);
  free_gfield_f(u_gauge_f);
 
  free_BCs(); 
  finalize_process();
  exit(0);
}
