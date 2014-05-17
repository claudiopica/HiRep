/*******************************************************************************
*
* Test of modules
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
#include "global.h"
#include "suN.h"
#include "suN_types.h"
#include "linear_algebra.h"
#include "inverters.h"
#include "representation.h"
#include "utils.h"
#include "communications.h"
#include "logger.h"
#include "memory.h"


static void random_g(suNg_field* g)
{
  _MASTER_FOR(&glattice,ix) {
    random_suNg(_FIELD_AT(g,ix));
  }
}

static void transform(suNg_field* gtransf, suNg_field* gfield)
{
  _MASTER_FOR(&glattice,ix) {
    for (int mu=0;mu<4;mu++) {
      int iy=iup(ix,mu);
      suNg * u=_4FIELD_AT(gfield,ix,mu);
      suNg v;
      _suNg_times_suNg_dagger(v,*u,*_FIELD_AT(gtransf,iy));
      _suNg_times_suNg(*u,*_FIELD_AT(gtransf,ix),v);
    }
  }
  
  start_gf_sendrecv(gfield);
  complete_gf_sendrecv(gfield);
}


int main(int argc,char *argv[])
{


  char pame[256];
  
  setup_process(&argc,&argv);
  
  logger_setlevel(0,100); /* log all */
  if (PID!=0) {
    logger_disable();}
  else{
    sprintf(pame,">out_%d",PID); logger_stdout(pame);
    sprintf(pame,"err_%d",PID); freopen(pame,"w",stderr);
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

  
  BCs_pars_t BCs_pars = {
    .fermion_twisting_theta = {0.,0.,0.,0.},
    .gauge_boundary_improvement_cs = 1.,
    .gauge_boundary_improvement_ct = 1.,
    .chiSF_boundary_improvement_ds = 1.,
    .SF_BCs = 0
  };
  init_BCs(&BCs_pars);

 


  suNg_field* u_gauge_1;
  suNg_field* u_gauge_2;
  suNg_field* g;

  double norm, tmp;
  double weight[3]={0.6,0.4,0.3};
  
  
 
  lprintf("MAIN",0,"Gauge group: SU(%d)\n",NG);
  lprintf("MAIN",0,"global size is %dx%dx%dx%d\n",GLB_T,GLB_X,GLB_Y,GLB_Z);
  lprintf("MAIN",0,"proc grid is %dx%dx%dx%d\n",NP_T,NP_X,NP_Y,NP_Z);

  lprintf("MAIN",0,"local size is %dx%dx%dx%d\n",T,X,Y,Z);
  lprintf("MAIN",0,"extended local size is %dx%dx%dx%d\n",T_EXT,X_EXT,Y_EXT,Z_EXT);
   
  lprintf("MAIN",0,"HYP smearing weights: %f %f %f\n",weight[0],weight[1],weight[2]);

  /* alloc global gauge fields */
  u_gauge=alloc_gfield(&glattice);
  u_gauge_1=alloc_gfield(&glattice);
  u_gauge_2=alloc_gfield(&glattice);
  g=alloc_gtransf(&glattice);

 
  lprintf("MAIN",0,"Generating a random gauge field... ");
  fflush(stdout);
  random_u(u_gauge);
  start_gf_sendrecv(u_gauge);
  complete_gf_sendrecv(u_gauge);
  lprintf("MAIN",0,"done.\n");
  
  lprintf("MAIN",0,"Generating a random gauge transf... ");
  fflush(stdout);
  random_g(g);
  start_gt_sendrecv(g);
  complete_gt_sendrecv(g);
  lprintf("MAIN",0,"done.\n");

  HYP_smearing(u_gauge_1,u_gauge,weight);
  transform(g,u_gauge_1);
  
  transform(g,u_gauge);
  HYP_smearing(u_gauge_2,u_gauge,weight);
  
  lprintf("MAIN",0,"Gauge covariance of the HYP smearing operator...\n");
  
  norm=0.;
  _MASTER_FOR_SUM(&glattice,ix,norm) {
    for (int mu=0;mu<4;mu++) {
      _suNg_sub_assign(*_4FIELD_AT(u_gauge_1,ix,mu),*_4FIELD_AT(u_gauge_2,ix,mu));
      _suNg_sqnorm(tmp,*_4FIELD_AT(u_gauge_1,ix,mu));
      norm+=tmp;
    }
  }
  global_sum(&norm,1);
  norm=sqrt(norm/(GLB_T*GLB_X*GLB_Y*GLB_Z*4*(NG*NG-1)*2));

  lprintf("MAIN",0,"Mean difference (must be small): %e\n",norm);
  free_gfield(u_gauge);
  free_gfield(u_gauge_1);
  free_gfield(u_gauge_2);
  free_gfield(g);
  finalize_process();

  exit(0);
}
