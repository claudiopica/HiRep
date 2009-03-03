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
  _DECLARE_INT_ITERATOR(ix);

  _MASTER_FOR(&glattice,ix)
    random_suNg(_FIELD_AT(g,ix));
}


static void transform(suNg_field* gtransf, suNg_field* gfield)
{
  _DECLARE_INT_ITERATOR(ix);
  int iy,mu;
  suNg *u,v;

  _MASTER_FOR(&glattice,ix) {
    for (mu=0;mu<4;mu++) {
      iy=iup(ix,mu);
      u=_4FIELD_AT(gfield,ix,mu);
      _suNg_times_suNg_dagger(v,*u,*_FIELD_AT(gtransf,iy));
      _suNg_times_suNg(*u,*_FIELD_AT(gtransf,ix),v);
    }
  }

  start_gf_sendrecv(gfield);
  complete_gf_sendrecv(gfield);
}


int main(int argc,char *argv[])
{
  suNg_field* u_gauge_1;
  suNg_field* u_gauge_2;
  suNg_field* g;
  int level,seed;
  double norm, tmp;
  double weight[3]={0.6,0.4,0.3};
  _DECLARE_INT_ITERATOR(ix);
  int mu;

  
  GLB_T=GLB_X=GLB_Y=GLB_Z=8;
  NP_T=NP_X=NP_Y=NP_Z=1;
  
  /* setup communication geometry */
  if (geometry_init() == 1) {
    finalize_process();
    return 0;
  }

  lprintf("MAIN",0,"Gauge group: SU(%d)\n",NG);
  lprintf("MAIN",0,"global size is %dx%dx%dx%d\n",GLB_T,GLB_X,GLB_Y,GLB_Z);
  lprintf("MAIN",0,"proc grid is %dx%dx%dx%d\n",NP_T,NP_X,NP_Y,NP_Z);
  lprintf("MAIN",0,"Fermion boundary conditions: %.2f,%.2f,%.2f,%.2f\n",bc[0],bc[1],bc[2],bc[3]);

  /* setup lattice geometry */
  geometry_mpi_eo();
  /* test_geometry_mpi_eo(); */

  lprintf("MAIN",0,"local size is %dx%dx%dx%d\n",T,X,Y,Z);
  lprintf("MAIN",0,"extended local size is %dx%dx%dx%d\n",T_EXT,X_EXT,Y_EXT,Z_EXT);
   
  level=0;
  seed=123;
  rlxs_init(level,seed);
  lprintf("MAIN",0,"ranlux: level = %d, seed = %d\n\n",level,seed); 
  fflush(stdout);

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
  _MASTER_FOR(&glattice,ix) {
    for (mu=0;mu<4;mu++) {
      _suNg_sub_assign(*_4FIELD_AT(u_gauge_1,ix,mu),*_4FIELD_AT(u_gauge_2,ix,mu));
      _suNg_sqnorm(tmp,*_4FIELD_AT(u_gauge_1,ix,mu));
      norm+=tmp;
    }
  }
  norm=sqrt(norm/(GLB_T*GLB_X*GLB_Y*GLB_Z*4*(NG*NG-1)*2));

  lprintf("MAIN",0,"Mean difference (must be small): %e\n",norm);
  
  exit(0);
}
