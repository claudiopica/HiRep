/*******************************************************************************
*
* Gauge covariance of the Dirac operator
*
*******************************************************************************/

#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "io.h"
#include "update.h"
#include "geometry.h"
#include "global.h"
#include "logger.h"
#include "random.h"
#include "memory.h"
#include "dirac.h"
#include "linear_algebra.h"
#include "representation.h"
#include "communications.h"

static double hmass=0.1;
static suNg_field *g;


static void D(spinor_field *out, spinor_field *in){
   Dphi(hmass,out,in);
}

static void random_g(void)
{
  _MASTER_FOR(&glattice,ix) {
      random_suNg(_FIELD_AT(g,ix));
  }
}

static void transform_u(void)
{
   _MASTER_FOR(&glattice,ix) {
     suNg v;
     for (int mu=0;mu<4;mu++) {
         int iy=iup(ix,mu);
         suNg *u=pu_gauge(ix,mu);
         _suNg_times_suNg_dagger(v,*u,*_FIELD_AT(g,iy));
         _suNg_times_suNg(*u,*_FIELD_AT(g,ix),v);
      }
   }
   
   start_gf_sendrecv(u_gauge);
   represent_gauge_field();
}

static void transform_s(spinor_field *out, spinor_field *in)
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
  char tmp[256];
  double sig,tau;
  spinor_field *s0,*s1,*s2,*s3;
  
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
  
  u_gauge=alloc_gfield(&glattice);
#ifndef REPR_FUNDAMENTAL
  u_gauge_f=alloc_gfield_f(&glattice);
#endif
  represent_gauge_field();
  
  /* allocate memory */
  g=alloc_gtransf(&glattice);
  s0=alloc_spinor_field_f(4,&glattice);
  s1=s0+1;
  s2=s1+1;
  s3=s2+1;
  
  lprintf("MAIN",0,"Generating a random gauge field... ");
  fflush(stdout);
  random_u(u_gauge);
  start_gf_sendrecv(u_gauge);
  represent_gauge_field();
  lprintf("MAIN",0,"done.\n");
  
  spinor_field_zero_f(s0);
  gaussian_spinor_field(&(s0[0]));
  tau = 1./sqrt(spinor_field_sqnorm_f(s0));
  spinor_field_mul_f(s0,tau,s0);
  
  lprintf("MAIN",0,"Generating a random gauge transf... ");
  fflush(stdout);
  random_g();
  start_gt_sendrecv(g);
  complete_gt_sendrecv(g);
  lprintf("MAIN",0,"done.\n");
  
  lprintf("MAIN",0,"Gauge covariance of the Dirac operator:\n");
  
  D(s1,s0);
  
  transform_s(s2,s1);
  
  
  transform_s(s3,s0);
  
  transform_u();
  
  spinor_field_zero_f(s1);
  D(s1,s3);
  
  
  spinor_field_mul_add_assign_f(s1,-1.0,s2);
  sig=spinor_field_sqnorm_f(s1);
  
  lprintf("MAIN",0,"Maximal normalized difference = %.2e\n",sqrt(sig));
  lprintf("MAIN",0,"(should be around 1*10^(-15) or so)\n\n");
  
  free_gfield(u_gauge);
#ifndef REPR_FUNDAMENTAL
  free_gfield_f(u_gauge_f);
#endif
  free_spinor_field_f(s0);
  
  free_gtransf(g);
  
  finalize_process();
  exit(0);
}
