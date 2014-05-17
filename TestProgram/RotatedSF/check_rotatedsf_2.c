/*******************************************************************************
*
* Gauge covariance of the SF observables
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
static suNg_field *g;
/* double M_PI=3.141592653589793238462643383279502884197; */

static void random_g(void)
{
  
  _MASTER_FOR(&glattice,ix) {
    random_suNg(_FIELD_AT(g,ix));
  }
  
  if(COORD[0] == 0) {
    for (int ix1=0;ix1<X;++ix1) for (int iy1=0;iy1<Y;++iy1) for (int iz1=0;iz1<Z;++iz1){
/*      ix=ipt(2,ix1,iy1,iz1);
      _suNg_unit(*_FIELD_AT(g,ix));*/
      int ix=ipt(1,ix1,iy1,iz1);
      _suNg_unit(*_FIELD_AT(g,ix));
      ix=ipt(0,ix1,iy1,iz1);
      _suNg_unit(*_FIELD_AT(g,ix));
    }
  }

  if(COORD[0] == NP_T -1) {
    for (int ix1=0;ix1<X;++ix1) for (int iy1=0;iy1<Y;++iy1) for (int iz1=0;iz1<Z;++iz1){
      int ix=ipt(T-1,ix1,iy1,iz1);
      _suNg_unit(*_FIELD_AT(g,ix));
/*      ix=ipt(T-2,ix1,iy1,iz1);
      _suNg_unit(*_FIELD_AT(g,ix));
      ix=ipt(T-3,ix1,iy1,iz1);
      _suNg_unit(*_FIELD_AT(g,ix));*/
    }
  }
}

static void transform_u(void)
{
  _MASTER_FOR(&glattice,ix) {
    for (int mu=0;mu<4;mu++) {
      int iy=iup(ix,mu);
      suNg *u=pu_gauge(ix,mu);
      suNg v;
      
      _suNg_times_suNg_dagger(v,*u,*_FIELD_AT(g,iy));
      _suNg_times_suNg(*u,*_FIELD_AT(g,ix),v);
    }
  }
  
  start_gf_sendrecv(u_gauge);
  represent_gauge_field();
}

int main(int argc,char *argv[])
{
  
  setup_process(&argc,&argv);

  double acc=1.e-20;
  
  BCs_pars_t BCs_pars = {
    .fermion_twisting_theta = {0.,M_PI/5.,M_PI/5.,M_PI/5.},
    .gauge_boundary_improvement_cs = 1.,
    .gauge_boundary_improvement_ct = 1.,
    .chiSF_boundary_improvement_ds = 0.5,
    .SF_BCs = 0
  };
  
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
  
  
  lprintf("MAIN",0,"Gauge group: SU(%d)\n",NG);
  lprintf("MAIN",0,"Fermion representation: dim = %d\n",NF);
  lprintf("MAIN",0,"The lattice size is %dx%dx%dx%d\n",T,X,Y,Z);
  lprintf("MAIN",0,"The lattice global size is %dx%dx%dx%d\n",GLB_T,GLB_X,GLB_Y,GLB_Z);
  lprintf("MAIN",0,"The lattice borders are (%d,%d,%d,%d)\n",T_BORDER,X_BORDER,Y_BORDER,Z_BORDER);
  lprintf("MAIN",0,"\n");
  fflush(stdout);
  
  double mass=0.0;
  _update_par.SF_zf=6.;
  _update_par.SF_ds=3.;
  _update_par.SF_sign=1;

  
  u_gauge=alloc_gfield(&glattice);
  u_gauge_f=alloc_gfield_f(&glattice);
  random_u(u_gauge);
  apply_BCs_on_fundamental_gauge_field();

  represent_gauge_field();
 

  g=alloc_gtransf(&glattice);
  random_g();
  start_gt_sendrecv(g);
  complete_gt_sendrecv(g);

  lprintf("MAIN",0,"Plaquette before the random gauge transf %f\n",avr_plaquette());
  SF_PCAC_wall_mass(mass,acc);

  transform_u();

  lprintf("MAIN",0,"Plaquette after the random gauge transf %f\n",avr_plaquette());
  SF_PCAC_wall_mass(mass,acc);

  
   


  free_gfield(u_gauge);
  free_gfield_f(u_gauge_f);
  
  free_gtransf(g);
  
  finalize_process();
  exit(0);
}
