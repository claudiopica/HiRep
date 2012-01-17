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
#include "communications_flt.h"

static double hmass=0.1;
static suNg_field_flt *g;

double sfdiff (spinor_field_flt* sf){
  spinor_field *tmp;
  double res;
  tmp=alloc_spinor_field_f(2, sf->type);
  alloc_spinor_field_f_gpu(2, tmp);
  assign_s2sd(&tmp[0], sf);
  spinor_field_copy_from_gpu_f_flt(sf);
  assign_s2sd(&tmp[1], sf);
  spinor_field_sub_f_cpu(&tmp[1],&tmp[0],&tmp[1]);
  res= spinor_field_sqnorm_f_cpu(&tmp[1]);
  assign_sd2s(sf,&tmp[0]);
  free_spinor_field_gpu(tmp);
  free_spinor_field(tmp);
}



static void D(spinor_field_flt *out, spinor_field_flt *in){
   Dphi_flt(hmass,out,in);
}

static void random_g(void)
{
   _DECLARE_INT_ITERATOR(ix);

   _MASTER_FOR(&glattice,ix)
      random_suNg_flt(_FIELD_AT(g,ix));				// This one might not exist
}

static void transform_u(void)
{
   _DECLARE_INT_ITERATOR(ix);
   int iy,mu;
   suNg_flt *u,v;

   _MASTER_FOR(&glattice,ix) {
      for (mu=0;mu<4;mu++) {
         iy=iup(ix,mu);
         u=pu_gauge_flt(ix,mu);
         _suNg_times_suNg_dagger(v,*u,*_FIELD_AT(g,iy));
         _suNg_times_suNg(*u,*_FIELD_AT(g,ix),v);
      }
   }
   
   start_gf_sendrecv_flt(u_gauge_flt);
   represent_gauge_field_flt();
}

static void transform_s(spinor_field_flt *out, spinor_field_flt *in)
{
   _DECLARE_INT_ITERATOR(ix);
   suNf_flt gfx;
   suNf_spinor_flt *r,*s;

   _MASTER_FOR(&glattice,ix) {
      s = _FIELD_AT(in,ix);
      r = _FIELD_AT(out,ix);
      
      _group_represent2_flt(&gfx,_FIELD_AT(g,ix));

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
  spinor_field_flt *s0,*s1,*s2,*s3;
  gpu_timer t1;
  float elapsed, gflops;
  int i;
  int n_times=50;
  
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
  
  u_gauge_flt=alloc_gfield_flt(&glattice);
  alloc_gfield_flt_gpu(u_gauge_flt);
#if (!defined(REPR_FUNDAMENTAL) && !defined(WITH_QUATERNIONS)) || defined(ROTATED_SF) 
  u_gauge_f_flt=alloc_gfield_f_flt(&glattice);
  alloc_gfield_f_flt_gpu(u_gauge_f_flt);
#endif
  /* allocate memory */
  s0=alloc_spinor_field_f_flt(1,&glat_even);
  alloc_spinor_field_f_flt_gpu(1,s0);
  s1=alloc_spinor_field_f_flt(1,&glat_even);
  alloc_spinor_field_f_flt_gpu(1,s1);

  s2=alloc_spinor_field_f_flt(1,&glat_odd);
  alloc_spinor_field_f_flt_gpu(1,s2);
  s3=alloc_spinor_field_f_flt(1,&glat_odd);
  alloc_spinor_field_f_flt_gpu(1,s3);
  
  gaussian_spinor_field_flt(s0);
  spinor_field_copy_to_gpu_f_flt(s0);

  gaussian_spinor_field_flt(s1);
  spinor_field_copy_to_gpu_f_flt(s1);

  gaussian_spinor_field_flt(s2);
  spinor_field_copy_to_gpu_f_flt(s2);

  gaussian_spinor_field_flt(s3);
  spinor_field_copy_to_gpu_f_flt(s3);

  res1 = sfdiff(s0);
  res2 = sfdiff(s2);  
  lprintf("LA TEST",0,"Copy works if 0=%1.10g and 0=%1.10 .\n",res1,res2);

	  
  lprintf("MAIN",0,"Generating a random gauge field... ");
  fflush(stdout);
  random_u_flt(u_gauge_flt);
  gfield_copy_to_gpu_flt(u_gauge_flt);
  start_gf_sendrecv_flt(u_gauge_flt);
  represent_gauge_field_flt();
  gfield_copy_to_gpu_f_flt(u_gauge_f_flt);

  //Diracoperator with zero mass
  
  lprintf("LA TEST",0,"Testing dirac operator with mass=0\n");
  lprintf("LA TEST",0,"Even lattice\n");

  // Norm before operation
  res_gpu = spinor_field_sqnorm_f_flt(s0);
  res_cpu = spinor_field_sqnorm_f_flt_cpu(s0);
  lprintf("LA TEST",0,"Before operator even, \nsqnorm(qpu)=%1.10g, sqnorm(cpu)=%1.10g\nsqnorm(gpu-cpu)=%1.10g\n",res_gpu,res_cpu,sfdiff(s0));

  Dphi_eopre_flt(0,s1,s0);
  Dphi_eopre_flt_cpu(0,s1,s0);

  res_gpu = spinor_field_sqnorm_f_flt(s1);
  res_cpu = spinor_field_sqnorm_f_flt_cpu(s1);
  
  res1 = sfdiff(s0); 
  res2 = sfdiff(s1);

  lprintf("LA TEST",0,"Result, mass=0, \nsqnorm(qpu)=%1.10g, sqnorm(cpu)=%1.10g,\nsqnorm(gpu-cpu)= %1.10g (check %1.10g=0?),\n",res_gpu,res_cpu,res2,res1);

  // Norm before operation
  res_gpu = spinor_field_sqnorm_f_flt(s2);
  res_cpu = spinor_field_sqnorm_f_flt_cpu(s2);
  lprintf("LA TEST",0,"Before operator odd, \nsqnorm(qpu)=%1.10g, sqnorm(cpu)=%1.10g\nsqnorm(gpu-cpu)=%1.10g\n",res_gpu,res_cpu,sfdiff(s2));
	
  lprintf("MAIN",0,"done even.\n");
  lprintf("LA TEST",0,"Odd lattice\n");

  Dphi_oepre_flt(0,s2,s3);
  Dphi_oepre_flt_cpu(0,s2,s3);

  res_gpu = spinor_field_sqnorm_f_flt(s3);
  res_cpu = spinor_field_sqnorm_f_flt_cpu(s3);
  
  res1 = sfdiff(s2); 
  res2 = sfdiff(s3);

  lprintf("LA TEST",0,"Result, mass=0, \nsqnorm(qpu)=%1.10g, sqnorm(cpu)=%1.10g,\nsqnorm(gpu-cpu)= %1.10g (check %1.10g=0?),\n",res_gpu,res_cpu,res2,res1);
  lprintf("MAIN",0,"done odd.\n");

  //Diracoperator with zero m=0.13
  
  // Norm before operation
  lprintf("LA TEST",0,"Testing dirac operator with mass=0.13\n");
  lprintf("LA TEST",0,"Even lattice\n");

  Dphi_eopre_flt(0.13,s1,s0);
  Dphi_eopre_flt_cpu(0.13,s1,s0);

  res_gpu = spinor_field_sqnorm_f_flt(s1);
  res_cpu = spinor_field_sqnorm_f_flt_cpu(s1);
  
  res1 = sfdiff(s0); 
  res2 = sfdiff(s1);

  lprintf("LA TEST",0,"Result, mass=0, \nsqnorm(qpu)=%1.10g, sqnorm(cpu)=%1.10g,\nsqnorm(gpu-cpu)= %1.10g (check %1.10g=0?),\n",res_gpu,res_cpu,res2,res1);

  // Norm before operation
  lprintf("MAIN",0,"done even.\n");

  lprintf("LA TEST",0,"Odd lattice\n");
  Dphi_oepre_flt(0.13,s2,s3);
  Dphi_oepre_flt_cpu(0.13,s2,s3);

  res_gpu = spinor_field_sqnorm_f_flt(s3);
  res_cpu = spinor_field_sqnorm_f_flt_cpu(s3);
  
  res1 = sfdiff(s2); 
  res2 = sfdiff(s3);

  lprintf("LA TEST",0,"Result, mass=0, \nsqnorm(qpu)=%1.10g, sqnorm(cpu)=%1.10g,\nsqnorm(gpu-cpu)= %1.10g (check %1.10g=0?),\n",res_gpu,res_cpu,res2,res1);
  lprintf("MAIN",0,"done odd.\n");


  //Gamma_5xDiracoperator with m=0.13
  
  // Norm before operation
  lprintf("LA TEST",0,"Testing gamma_5 x dirac operator with mass=0.13");
  lprintf("LA TEST",0,"Even lattice\n");

  g5Dphi_eopre_flt(0.13,s1,s0);
  g5Dphi_eopre_flt_cpu(0.13,s1,s0);

  res_gpu = spinor_field_sqnorm_f_flt(s1);
  res_cpu = spinor_field_sqnorm_f_flt_cpu(s1);
  
  res1 = sfdiff(s0); 
  res2 = sfdiff(s1);

  lprintf("LA TEST",0,"Result, mass=0, \nsqnorm(qpu)=%1.10g, sqnorm(cpu)=%1.10g,\nsqnorm(gpu-cpu)= %1.10g (check %1.10g=0?),\n",res_gpu,res_cpu,res2,res1);

  // Norm before operation
  lprintf("MAIN",0,"done even.\n");


  //Gamma_5xDiracoperator^2 with m=0.13
  
  // Norm before operation
  lprintf("LA TEST",0,"Testing gamma_5 x (dirac operator)^2 with mass=0.13");
  lprintf("LA TEST",0,"Even lattice\n");

  g5Dphi_eopre_sq_flt(0.13,s1,s0);
  g5Dphi_eopre_sq_flt_cpu(0.13,s1,s0);

  res_gpu = spinor_field_sqnorm_f_flt(s1);
  res_cpu = spinor_field_sqnorm_f_flt_cpu(s1);
  
  res1 = sfdiff(s0); 
  res2 = sfdiff(s1);

  lprintf("LA TEST",0,"Result, mass=0, \nsqnorm(qpu)=%1.10g, sqnorm(cpu)=%1.10g,\nsqnorm(gpu-cpu)= %1.10g (check %1.10g=0?),\n",res_gpu,res_cpu,res2,res1);

  // Norm before operation
  lprintf("MAIN",0,"done even.\n");



  lprintf("LA TEST",0,"Calculating Diracoperator %d times on even lattice.\n",n_times);

  t1 = gpuTimerStart();

  for (i=0;i<n_times;++i){
    Dphi_eopre_flt(0,s1,s0);
  }

  elapsed = gpuTimerStop(t1);
  lprintf("LA TEST",0,"Time: %1.10gms\n",elapsed);

  gflops=n_times*GLB_T*GLB_X*GLB_Y*GLB_Z*744./elapsed/1.e6;   // 536 //240
  lprintf("LA TEST",0,"GFLOPS: %1.6g\n\n",gflops);

  gflops=8.;
  gflops=n_times*GLB_T*GLB_X*GLB_Y*GLB_Z*((24.+28*gflops)*4.+gflops*4.)/elapsed/1.e6; 
  //gflops=n_times*GLB_T*GLB_X*GLB_Y*GLB_Z*328./elapsed/1.e6; 
  //gflops=n_times*GLB_T*GLB_X*GLB_Y*GLB_Z*212./elapsed/1.e6; 
  lprintf("LA TEST",0,"BAND: %1.6g\n\n",gflops);

  lprintf("LA TEST",0,"Calculating Diracoperator %d times on odd lattice.\n",n_times);

  t1 = gpuTimerStart();

  for (i=0;i<n_times;++i){
    Dphi_oepre_flt(0,s3,s2);
  }

  elapsed = gpuTimerStop(t1);
  lprintf("LA TEST",0,"Time: %1.10gms\n",elapsed);

  gflops=n_times*GLB_T*GLB_X*GLB_Y*GLB_Z*744./elapsed/1.e6;   // 536 //240
  lprintf("LA TEST",0,"GFLOPS: %1.6g\n\n",gflops);

  gflops=8.;
  gflops=n_times*GLB_T*GLB_X*GLB_Y*GLB_Z*((24.+28*gflops)*4.+gflops*4.)/elapsed/1.e6; 
  lprintf("LA TEST",0,"BAND: %1.6g\n\n",gflops);

  lprintf("LA TEST",0,"Calculating gamma_5 x Diracoperator %d times on even lattice.\n",n_times);

  t1 = gpuTimerStart();

  for (i=0;i<n_times;++i){
    g5Dphi_eopre_flt(0,s1,s0);
  }

  elapsed = gpuTimerStop(t1);
  lprintf("LA TEST",0,"Time: %1.10gms\n",elapsed);

  gflops=n_times*GLB_T*GLB_X*GLB_Y*GLB_Z*(744.+24.)/elapsed/1.e6;   // 536 //240
  lprintf("LA TEST",0,"GFLOPS: %1.6g\n\n",gflops);

  gflops=8.;
  gflops=n_times*GLB_T*GLB_X*GLB_Y*GLB_Z*((24.+28*gflops)*4.+gflops*4.+24.)/elapsed/1.e6; 
  lprintf("LA TEST",0,"BAND: %1.6g\n\n",gflops);

  lprintf("LA TEST",0,"Calculating (gamma_5 x Diracoperator)^2 %d times on even lattice.\n",n_times);

  t1 = gpuTimerStart();

  for (i=0;i<n_times;++i){
    g5Dphi_eopre_sq_flt(0,s1,s0);
  }

  elapsed = gpuTimerStop(t1);
  lprintf("LA TEST",0,"Time: %1.10gms\n",elapsed);

  gflops=2*n_times*GLB_T*GLB_X*GLB_Y*GLB_Z*(744.+24.)/elapsed/1.e6;   // 536 //240
  lprintf("LA TEST",0,"GFLOPS: %1.6g\n\n",gflops);

  gflops=8.;
  gflops=2*n_times*GLB_T*GLB_X*GLB_Y*GLB_Z*((24.+28*gflops)*4.+gflops*4.+24.)/elapsed/1.e6; 
  lprintf("LA TEST",0,"BAND: %1.6g\n\n",gflops);




/*  t1 = gpuTimerStart();
  for (i=0;i<n_times;++i){
    spinor_field_sub_f(s2,s0,s1);
  }
  elapsed = gpuTimerStop(t1);
  gflops=n_times*GLB_T*GLB_X*GLB_Y*GLB_Z*24.*2./elapsed/1.e6; 
    lprintf("LA TEST",0,"SQnorm GFLOPS: %1.4g\n\n",gflops);
 

*/


  lprintf("LA TEST",0,"DONE!\n\n");
  free_spinor_field_flt(s0);
  free_spinor_field_flt_gpu(s0);

  free_spinor_field_flt(s1);
  free_spinor_field_flt_gpu(s1);
  
  
  finalize_process();
  exit(0);
}
