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

static suNg_field *g;
rhmc_par _update_par={0};
/* double M_PI=3.141592653589793238462643383279502884197; */

static void random_g(void)
{
  _MASTER_FOR(&glattice,ix) {
    random_suNg(_FIELD_AT(g,ix));
  }
/*  if(COORD[0] == 0) {
    for (ix1=0;ix1<X;++ix1) for (iy1=0;iy1<Y;++iy1) for (iz1=0;iz1<Z;++iz1){
      ix=ipt(2,ix1,iy1,iz1);
      _suNg_unit(*_FIELD_AT(g,ix));
      ix=ipt(1,ix1,iy1,iz1);
      _suNg_unit(*_FIELD_AT(g,ix));
      ix=ipt(0,ix1,iy1,iz1);
      _suNg_unit(*_FIELD_AT(g,ix));
    }
  }
  
  if(COORD[0] == NP_T -1) {
    for (ix1=0;ix1<X;++ix1) for (iy1=0;iy1<Y;++iy1) for (iz1=0;iz1<Z;++iz1){
      ix=ipt(T-1,ix1,iy1,iz1);
      _suNg_unit(*_FIELD_AT(g,ix));
      ix=ipt(T-2,ix1,iy1,iz1);
      _suNg_unit(*_FIELD_AT(g,ix));
      ix=ipt(T-3,ix1,iy1,iz1);
      _suNg_unit(*_FIELD_AT(g,ix));
    }
  }*/
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

static void transform_s(spinor_field *in,spinor_field *out)
{
  _MASTER_FOR(&glattice,ix) {
    suNf_spinor *s = _FIELD_AT(in,ix);
    suNf_spinor *r = _FIELD_AT(out,ix);
    suNf gfx;
    
    _group_represent2(&gfx,_FIELD_AT(g,ix));
    
    _suNf_multiply(r->c[0],gfx,s->c[0]);
    _suNf_multiply(r->c[1],gfx,s->c[1]);
    _suNf_multiply(r->c[2],gfx,s->c[2]);
    _suNf_multiply(r->c[3],gfx,s->c[3]);
  }
}

int main(int argc,char *argv[])
{
  char sbuf[128];

  double acc=1.e-26,tau,sig;
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


  logger_setlevel(0,10000); /* log all */
  logger_map("DEBUG","debug");
  
  if (PID!=0) { logger_disable(); }

  if (PID==0) {
    sprintf(sbuf,">>out_%d",PID);  logger_stdout(sbuf); 
    sprintf(sbuf,"err_%d",PID); freopen(sbuf,"w",stderr);
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
  gaussian_spinor_field(&(s0[0]));
  tau = 1./sqrt(spinor_field_sqnorm_f(s0));
  spinor_field_mul_f(s0,tau,s0);
    
  double mass=0.0;
  _update_par.SF_zf=6.;
  _update_par.SF_ds=3.;
  _update_par.SF_sign=1;

  random_u(u_gauge);

  apply_BCs_on_fundamental_gauge_field();


  represent_gauge_field();
 
  g=alloc_gtransf(&glattice);
  random_g();
  start_gt_sendrecv(g);
  complete_gt_sendrecv(g);

 
  
  SF_quark_propagator(s0, mass, s1 , acc);
 
  transform_u();
  transform_s(s0, s2);
  transform_s(s1, s0);


  SF_quark_propagator(s2, mass, s1 , acc);

  spinor_field_mul_add_assign_f(s1,-1.0,s0);
  sig=spinor_field_sqnorm_f(s1);
  
  lprintf("MAIN",0,"Maximal normalized difference = %.2e\n",sqrt(sig));
  lprintf("MAIN",0,"(should be around 1*10^(-10) or so)\n\n");


  free_gfield(u_gauge);
  free_gfield_f(u_gauge_f);
  free_spinor_field_f(s0);

  free_gtransf(g);
  
  finalize_process();
  exit(0);
}
