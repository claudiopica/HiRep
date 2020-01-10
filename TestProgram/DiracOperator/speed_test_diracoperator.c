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
#include "setup.h"
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
#include "utils.h"
#include "communications.h"

static double hmass=0.1;
static suNg_field *g;



//static void loc_D(spinor_field *out, spinor_field *in){
//   Dphi(hmass,out,in);
//}



int main(int argc,char *argv[])
{

  double res1,res2,res3,res_cpu,res_gpu;
  spinor_field *s0,*s1,*s2,*s3, *tmps;
  float elapsed, gflops;
  int flopsite, bytesite;
  int n_times=5000;
  int n_warmup=1000;
  struct timeval start, end, etime;
  
  setup_process(&argc,&argv);
  
  setup_gauge_fields();

  //geometry_mpi_eo();

  /* setup random numbers */
 // read_input(rlx_var.read,get_input_filename());
  //lprintf("MAIN",0,"RLXD [%d,%d]\n",rlx_var.rlxd_level,rlx_var.rlxd_seed+MPI_PID);
  lprintf("MAIN",0,"n_warmup = %d, n_times = %d \n",n_warmup,n_times);
  
  //rlxd_init(rlx_var.rlxd_level,rlx_var.rlxd_seed+MPI_PID); /* use unique MPI_PID to shift seeds */
  
  
//  lprintf("MAIN",0,"Allocating gauge field\n");  
//  u_gauge=alloc_gfield(&glattice);
//#if (!defined(REPR_FUNDAMENTAL) && !defined(WITH_QUATERNIONS)) || defined(ROTATED_SF) 
//  u_gauge_f=alloc_gfield_f(&glattice);
//#endif

  /* allocate memory */
  lprintf("MAIN",0,"Allocating spinor field\n");  
  s0=alloc_spinor_field_f(4,&glattice);
  s1=s0+1;
  s2=s1+1;
  s3=s2+1;
  

  lprintf("MAIN",0,"Randomizing spinor field...\n");  
  gaussian_spinor_field(s0);
  gaussian_spinor_field(s1);
  lprintf("LA_TEST",0,"un sito %lf\n",creal(s0->ptr[2].c[0].c[0]));

  //
  //#pragma omp parallel num_threads(1) default(shared)
#pragma omp parallel
  {
    res1=spinor_field_sqnorm_f(s0);

_OMP_BARRIER
    
    res2=spinor_field_sqnorm_f(s1);
  }
  lprintf("LA_TEST",0,"Square norm for check after gaussian spinor field %lf and %lf\n",res1,res2);
  
  
 
  lprintf("MAIN",0,"Generating a random gauge field... ");
  fflush(stdout);
  random_u(u_gauge);
  start_gf_sendrecv(u_gauge);
  represent_gauge_field();
 
  lprintf("MAIN",0,"done.\n");

  //Check speed diracoperator

  
#if defined(REPR_ADJOINT)
  flopsite=8*NF*(7+8*NF);
#else
  flopsite=8*NF*(7+16*NF);
#endif
  bytesite=36*sizeof(suNf_vector)+8*sizeof(suNf); //add integers for geometry indexes?
  
  lprintf("LA TEST",0,"Flop per site = %d\n",flopsite);
  lprintf("LA TEST",0,"Byte per site = %d\n",bytesite);


  //speed test Dirac operator
  lprintf("LA TEST",0,"Warmup application of the Diracoperator %d times.\n",n_warmup);
  
  
#pragma omp parallel default(shared)
  {
    for (int i=0;i<n_warmup;++i){ Dphi_(s1,s0); }
  }

  lprintf("LA TEST",0,"Calculating massless Diracoperator %d times.\n",n_times);
  
  gettimeofday(&start,0);
#pragma omp parallel default(shared)
  {
    for (int i=0;i<n_times;++i){ Dphi_(s1,s0); }
  }
  gettimeofday(&end,0);
  timeval_subtract(&etime,&end,&start);
  elapsed=etime.tv_sec*1000.+etime.tv_usec*0.001;
  lprintf("LA TEST",0,"Time: [%ld sec %ld usec]\n",etime.tv_sec,etime.tv_usec);
  
  
#pragma omp parallel num_threads(1) default(shared)
  {
    Dphi_(s2,s0);
    
    spinor_field_sub_assign_f(s2,s1);
    res1=spinor_field_sqnorm_f(s2);
#pragma omp barrier 
    
    res2=spinor_field_sqnorm_f(s1);
#pragma omp barrier 
    
    res3=spinor_field_sqnorm_f(s0);
    
  }
  lprintf("LA_TEST",0,"Square norm for check %lf originally %lf and from the begining %lf \n",res1,res2,res3);
  

  //lprintf("LA TEST",0,"Time: %1.10g ms\n",elapsed);
  gflops=((double)n_times*GLB_T*GLB_X*GLB_Y*GLB_Z*flopsite)/elapsed/1.e6;
  lprintf("LA TEST",0,"GFLOPS: %1.6g\n\n",gflops);
  gflops=(double)n_times*GLB_T*GLB_X*GLB_Y*GLB_Z*bytesite/elapsed/1.e6;
  lprintf("LA TEST",0,"BAND: %1.6g GB/s\n\n",gflops);
  lprintf("LA TEST",0,"DONE!");
  
  
  free_spinor_field_f(s0);
  
  finalize_process();
  exit(0);
}


