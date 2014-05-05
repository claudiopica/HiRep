/*******************************************************************************
*
* Does the SF_quark_propagator compute H^{-1}?
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

rhmc_par _update_par={0};
/* double M_PI=3.141592653589793238462643383279502884197; */


int main(int argc,char *argv[])
{
  double acc=1.e-20,tau,sig;
  spinor_field *s0,*s1,*s2,*s3;
  setup_process(&argc,&argv);
 
  BCs_pars_t BCs_pars = {
    .fermion_twisting_theta = {0.,M_PI/5.,M_PI/5.,M_PI/5.},
    .gauge_boundary_improvement_cs = 1.,
    .gauge_boundary_improvement_ct = 1.,
    .chiSF_boundary_improvement_ds = 0.5,
    .SF_BCs = 0
  };

 _update_par.mass=-0.3099358667;
  _update_par.SF_ds=0.5;
  _update_par.SF_sign=1;
  _update_par.SF_ct=1;
  _update_par.SF_zf=1.3;

  char tmp[256];

  logger_setlevel(0,100); /* log all */
  if (PID!=0) { 
    logger_disable();}
  else{
    sprintf(tmp,">out_%d",PID); logger_stdout(tmp);
    sprintf(tmp,"err_%d",PID); freopen(tmp,"w",stderr);
  }
  
  lprintf("MAIN",0,"PId =  %d [world_size: %d]\n\n",PID,WORLD_SIZE); 
  
  read_input(glb_var.read,"test_input");
  read_input(rlx_var.read,"test_input");

  rlxd_init(rlx_var.rlxd_level,rlx_var.rlxd_seed+MPI_PID);
  
  
  
  
  /* setup communication geometry */
  if (geometry_init() == 1) {
    finalize_process();
    return 0;
  }
  
  geometry_mpi_eo();

  init_BCs(&BCs_pars);

  
  u_gauge=alloc_gfield(&glattice);
  u_gauge_f=alloc_gfield_f(&glattice);
  s0=alloc_spinor_field_f(4,&glattice);
  s1=s0+1;
  s2=s1+1;
  s3=s2+1;
   
  spinor_field_zero_f(s0);
  gaussian_spinor_field(s0);

  tau = 1./sqrt(spinor_field_sqnorm_f(s0));
  spinor_field_mul_f(s0,tau,s0);

  double mass=0.0;
  _update_par.SF_zf=10.;
  _update_par.SF_ds=-3.;
  _update_par.SF_sign=1;
    

  random_u(u_gauge);

  apply_BCs_on_fundamental_gauge_field();

  represent_gauge_field();
 

  
  SF_quark_propagator(s0, mass, s1 , acc);
  g5Dphi(mass,s2,s1);

  spinor_field_mul_add_assign_f(s0,-1.0,s2);
  sig=spinor_field_sqnorm_f(s0);
  
  lprintf("MAIN",0,"Maximal normalized difference = %.2e\n",sqrt(sig));
  lprintf("MAIN",0,"(should be around 1*10^(-10) or so)\n\n");


  free_gfield(u_gauge);
  free_gfield_f(u_gauge_f);
  free_spinor_field_f(s0);
  
  finalize_process();
  exit(0);
}
