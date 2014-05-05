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


int main(int argc,char *argv[])
{
  setup_process(&argc,&argv);
  char sbuf[350];


  BCs_pars_t BCs_pars = {
    .fermion_twisting_theta = {0.,0.5,0.5,0.5},
    .gauge_boundary_improvement_cs = 1.,
    .gauge_boundary_improvement_ct = 1.,
    .chiSF_boundary_improvement_ds = 0.5,
    .SF_BCs = 0
  };
  


  double sig;
  spinor_field *s0,*s1,*s2,*s3,*s4,*s5;


  double kappa=0.1354990000;  
  double mass=1.0/(2*kappa)-4.0;
  double acc=1.e-20;

  _update_par.mass=1.0/(2*kappa)-4.0;
  _update_par.SF_ds=0.5;
  _update_par.SF_sign=1;
  _update_par.SF_ct=1;
  _update_par.SF_zf=1.3;


  logger_setlevel(0,99); /* log all */
  logger_map("DEBUG","debug");
  
  if (PID!=0) { logger_disable(); }

  if (PID==0) {
    sprintf(sbuf,"out_%d",PID);  logger_stdout(sbuf); 
    sprintf(sbuf,"err_%d",PID); freopen(sbuf,"w",stderr);
  }

  lprintf("MAIN",0,"PId =  %d [world_size: %d]\n\n",PID,WORLD_SIZE); 
  
  read_input(glb_var.read,"test_input");
  read_input(rlx_var.read,"test_input");

  rlxd_init(rlx_var.rlxd_level,rlx_var.rlxd_seed+MPI_PID);
  
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
  
  lprintf("MAIN",0,"Reading gauge configuration from file :gfield_sint.dat\n");
  read_gauge_field_nocheck("gfield_sint.dat");
  
  lprintf("MAIN",0,"Value of the plaquette in Sint's normalization %f\n",(avr_plaquette()*6*GLB_T*GLB_X*GLB_Y*GLB_Z*NG-3*2*GLB_X*GLB_Y*GLB_Z*NG)/(3*NG*GLB_X*GLB_Y*GLB_Z*(2*(GLB_T-2)-1)));
  
  lprintf("MAIN",0,"mass = %f\n",mass);
/*   lprintf("MAIN",0,"ds = %f\n",_update_par.SF_ds); */
/*   lprintf("MAIN",0,"zf = %f\n",_update_par.SF_zf); */
/*   lprintf("MAIN",0,"theta = %f\n",_update_par.SF_theta); */
/*   lprintf("MAIN",0,"sign = %d\n",_update_par.SF_sign); */

  represent_gauge_field(); 
  
  
  s0=alloc_spinor_field_f(6,&glattice);
  s1=s0+1;
  s2=s1+1;
  s3=s2+1;
  s4=s3+1;
  s5=s4+1;
 
  lprintf("MAIN",0,"Reading input spinor field configuration from file :test_volsource_zf1.3_ds_0.5\n");
  read_spinor_field_ascii("test_volsource_zf1.3_ds_0.5",s0);
  apply_BCs_on_spinor_field(s0);
 
  lprintf("MAIN",0,"Reading input spinor field configuration from file :test_Qvolsource_zf1.3_ds_0.5\n");
  read_spinor_field_ascii("test_Qvolsource_zf1.3_ds_0.5",s1);
  apply_BCs_on_spinor_field(s1);
 
  lprintf("MAIN",0,"Reading input spinor field configuration from file :test_Qinversevolsource_zf1.3_ds_0.5\n");
  read_spinor_field_ascii("test_Qinversevolsource_zf1.3_ds_0.5",s2);
  apply_BCs_on_spinor_field(s2);
  
  lprintf("MAIN",0,"Testing difference beetween the our implementation of the dirac operator on the delta source (test_volsource_zf1.3_ds_0.5) and the out spinor\n");
  g5Dphi(mass,s4,s0);
  apply_BCs_on_spinor_field(s4);

  spinor_field_mul_f(s4,2*(kappa)/(1+8*kappa),s4);

  spinor_field_mul_add_assign_f(s4,-1,s1);
  sig=spinor_field_sqnorm_f(s4);
  
  lprintf("MAIN",0,"Maximal normalized difference = %.2e\n",sqrt(sig/(GLB_T*GLB_X*GLB_Y*GLB_Z*24)));
  lprintf("MAIN",0,"(should be around 1*10^(-10) or so)\n\n");
  
  lprintf("MAIN",0,"Testing difference beetween the our implementation of the dirac on the out spinor of Stefan (test_Qinversevolsource_zf1.3_ds_0.5) and the in spinor\n");
  g5Dphi(mass,s3,s2);
  apply_BCs_on_spinor_field(s3);
  spinor_field_mul_f(s3,(2*kappa)/(1+8*kappa),s3);

  spinor_field_mul_add_assign_f(s3,-1,s0);
  
  sig=spinor_field_sqnorm_f(s3);
  
  lprintf("MAIN",0,"Maximal normalized difference = %.2e\n",sqrt(sig/(GLB_T*GLB_X*GLB_Y*GLB_Z*24)));
  lprintf("MAIN",0,"(should be around 1*10^(-10) or so)\n\n");

//  for(ix0=1;ix0<T;ix0++) for(ix1=0;ix1<X;ix1++) for(ix2=0;ix2<Y;ix2++) for(ix3=0;ix3<Z;ix3++) {
//    lprintf("MAIN",0,"%d %d %d %d\n",ix0,ix1,ix2,ix3);
//    i=ipt(ix0,ix1,ix2,ix3);
//    for(j=0;j<12;j++)
//      lprintf("SPINOR",0,"%.10f %.10f\n",((complex*)(_FIELD_AT(s3,i)))[j].re,((complex*)(_FIELD_AT(s3,i)))[j].im);
//  }

  
  lprintf("MAIN",0,"Testing the outspinor generated with Hirep Inverter\n");
  spinor_field_zero_f(s4);
  SF_quark_propagator(s0, mass, s4 , acc);

  sig=0.;

  spinor_field_mul_add_assign_f(s4,-(2*kappa)/(1+8*kappa),s2);
  sig=spinor_field_sqnorm_f(s4);
  
  lprintf("MAIN",0,"Maximal normalized difference = %.2e\n",sqrt(sig/(GLB_T*GLB_X*GLB_Y*GLB_Z*24)));
  lprintf("MAIN",0,"(should be around 1*10^(-10) or so)\n\n");


/*
  _update_par.SF_zf=2.3;
  _update_par.SF_ds=2.3;


  lprintf("MAIN",0,"Reading input spinor field configuration from file :test_volsource_zf1.3_ds_1.0\n");
  read_spinor_field_ascii("test_volsource_zf1.3_ds_1.0",s0);
  apply_BCs_on_spinor_field(s0);
 
  lprintf("MAIN",0,"Reading input spinor field configuration from file :test_Qvolsource_zf1.3_ds_1.0\n");
  read_spinor_field_ascii("test_Qvolsource_zf1.3_ds_1.0",s1);
  apply_BCs_on_spinor_field(s1);
 
  lprintf("MAIN",0,"Reading input spinor field configuration from file :test_Qinversevolsource_zf1.3_ds_1.0\n");
  read_spinor_field_ascii("test_Qinversevolsource_zf1.3_ds_1.0",s2);
  apply_BCs_on_spinor_field(s2);
  
  lprintf("MAIN",0,"Testing difference beetween the our implementation of the dirac operator on the delta source and the out spinor\n");
  g5Dphi(mass,s4,s0);
  apply_BCs_on_spinor_field(s4);

  spinor_field_mul_f(s4,2*(kappa)/(1+8*kappa),s4);

  spinor_field_mul_add_assign_f(s4,-1,s1);
  sig=spinor_field_sqnorm_f(s4);
  
  lprintf("MAIN",0,"Maximal normalized difference = %.2e\n",sqrt(sig/(GLB_T*GLB_X*GLB_Y*GLB_Z*24)));
  lprintf("MAIN",0,"(should be around 1*10^(-10) or so)\n\n");
  
  lprintf("MAIN",0,"Testing difference beetween the our implementation of the dirac on the out spinor of Stefan and the in spinor\n");
  g5Dphi(mass,s3,s2);
  apply_BCs_on_spinor_field(s3);
  spinor_field_mul_f(s3,(2*kappa)/(1+8*kappa),s3);

  spinor_field_mul_add_assign_f(s3,-1,s0);
  
  sig=spinor_field_sqnorm_f(s3);
  
  lprintf("MAIN",0,"Maximal normalized difference = %.2e\n",sqrt(sig/(GLB_T*GLB_X*GLB_Y*GLB_Z*24)));
  lprintf("MAIN",0,"(should be around 1*10^(-10) or so)\n\n");

  for(ix0=1;ix0<2;ix0++) for(ix1=0;ix1<1;ix1++) for(ix2=0;ix2<1;ix2++) for(ix3=0;ix3<1;ix3++) {
    lprintf("MAIN",0,"%d %d %d %d\n",ix0,ix1,ix2,ix3);
    i=ipt(ix0,ix1,ix2,ix3);
    for(j=0;j<12;j++)
      lprintf("SPINOR",0,"%.10f %.10f\n",((complex*)(_FIELD_AT(s3,i)))[j].re,((complex*)(_FIELD_AT(s3,i)))[j].im);
  }
*/
  free_gfield(u_gauge);
  free_gfield_f(u_gauge_f);
  free_spinor_field_f(s0);
  
  finalize_process();
  exit(0);
}
