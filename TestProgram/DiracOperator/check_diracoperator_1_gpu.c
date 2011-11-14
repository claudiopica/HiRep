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

double sfdiff (spinor_field* sf){
  spinor_field *tmp;
  double res;
  tmp=alloc_spinor_field_f(1, &glattice);
  alloc_spinor_field_f_gpu(1, tmp);
  spinor_field_copy_f_cpu(tmp,sf);
  spinor_field_copy_to_gpu_f(tmp);
  spinor_field_sub_f(tmp,tmp,sf);
  res= spinor_field_sqnorm_f(tmp);
  free_spinor_field_gpu(tmp);
  free_spinor_field(tmp);
  return res;
}



static void D(spinor_field *out, spinor_field *in){
   Dphi(hmass,out,in);
}

static void random_g(void)
{
   _DECLARE_INT_ITERATOR(ix);

   _MASTER_FOR(&glattice,ix)
      random_suNg(_FIELD_AT(g,ix));
}

static void transform_u(void)
{
   _DECLARE_INT_ITERATOR(ix);
   int iy,mu;
   suNg *u,v;

   _MASTER_FOR(&glattice,ix) {
      for (mu=0;mu<4;mu++) {
         iy=iup(ix,mu);
         u=pu_gauge(ix,mu);
         _suNg_times_suNg_dagger(v,*u,*_FIELD_AT(g,iy));
         _suNg_times_suNg(*u,*_FIELD_AT(g,ix),v);
      }
   }
   
   start_gf_sendrecv(u_gauge);
   represent_gauge_field();
}

static void transform_s(spinor_field *out, spinor_field *in)
{
   _DECLARE_INT_ITERATOR(ix);
   suNf gfx;
   suNf_spinor *r,*s;

   _MASTER_FOR(&glattice,ix) {
      s = _FIELD_AT(in,ix);
      r = _FIELD_AT(out,ix);
      
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
  double res1,res2,res_cpu,res_gpu;
  spinor_field *s0,*s1,*s2,*s3;
  
  setup_process(&argc,&argv);
  
  logger_setlevel(0,10000); /* log all */
  logger_map("DEBUG","debug");
#ifdef WITH_MPI
  sprintf(tmp,">out_%d",PID); logger_stdout(tmp);
  sprintf(tmp,"err_%d",PID); freopen(tmp,"w",stderr);
#endif
  
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
  
  
  lprintf("MAIN",0,"Gauge group: SU(%d)\n",NG);
  lprintf("MAIN",0,"Fermion representation: dim = %d\n",NF);
  lprintf("MAIN",0,"The lattice size is %dx%dx%dx%d\n",T,X,Y,Z);
  lprintf("MAIN",0,"The lattice global size is %dx%dx%dx%d\n",GLB_T,GLB_X,GLB_Y,GLB_Z);
  lprintf("MAIN",0,"The lattice borders are (%d,%d,%d,%d)\n",T_BORDER,X_BORDER,Y_BORDER,Z_BORDER);
  lprintf("MAIN",0,"\n");
  fflush(stdout);
  
  u_gauge=alloc_gfield(&glattice);
  alloc_gfield_gpu(u_gauge);
#if (!defined(REPR_FUNDAMENTAL) && !defined(WITH_QUATERNIONS)) || defined(ROTATED_SF) 
  u_gauge_f=alloc_gfield_f(&glattice);
  alloc_gfield_f_gpu(u_gauge_f);
#endif
  /* allocate memory */
  s0=alloc_spinor_field_f(4,&glattice);
  alloc_spinor_field_f_gpu(4,s0);
  s1=s0+1;
  s2=s1+1;
  s3=s2+1;
  
  gaussian_spinor_field(s0);
  spinor_field_copy_to_gpu_f(s0);

  gaussian_spinor_field(s1);
  spinor_field_copy_to_gpu_f(s1);

  lprintf("MAIN",0,"Generating a random gauge field... ");
  fflush(stdout);
  random_u(u_gauge);
  start_gf_sendrecv(u_gauge);
  represent_gauge_field();
  gfield_copy_to_gpu_f(u_gauge_f);

  lprintf("MAIN",0,"done.\n");

  
  Dphi_(s1,s0);
  Dphi__cpu(s1,s0);
  
  res1 = sfdiff(s0);
  res2 = sfdiff(s1);

  res_gpu = spinor_field_sqnorm_f(s1);
  res_cpu = spinor_field_sqnorm_f_cpu(s1);

  lprintf("LA TEST",0,"Check diracoperator, mass=0, \nsqnorm(qpu)=%1.10g, sqnorm(cpu)=%1.10g,\nsqnorm(gpu-cpu)= %1.10g (check %1.10g=0?) \n\n",res_gpu,res_cpu,res2,res1);


  Dphi(hmass,s1,s0);
  Dphi_cpu(hmass,s1,s0);

  res1 = sfdiff(s0);
  res2 = sfdiff(s1);

  res_gpu = spinor_field_sqnorm_f(s1);
  res_cpu = spinor_field_sqnorm_f_cpu(s1);

  lprintf("LA TEST",0,"Check diracoperator, mass=%1.2g, \nsqnorm(qpu)=%1.10g, sqnorm(cpu)=%1.10g,\nsqnorm(gpu-cpu)= %1.10g (check %1.10g=0?) \n\n",hmass,res_gpu,res_cpu,res2,res1);



  free_spinor_field(s0);
  free_spinor_field_gpu(s0);
  
  
  finalize_process();
  exit(0);
}
