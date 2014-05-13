/*******************************************************************************
*
* Check gauge covariance of the HYP smearing
*
*******************************************************************************/

#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "io.h"
#include "geometry.h"
#include "global.h"
#include "logger.h"
#include "random.h"
#include "memory.h"
#include "linear_algebra.h"
#include "communications.h"
#include "observables.h"
#include "error.h"
#include "utils.h"


static void random_g(suNg_field* g) {
   _DECLARE_INT_ITERATOR(ix);

   _MASTER_FOR(g->type,ix)
      random_suNg(_FIELD_AT(g,ix));
}

static void transform_u(suNg_field* out, suNg_field* in, suNg_field* g) {
   _DECLARE_INT_ITERATOR(ix);
   int iy,mu;
   suNg v;

   _MASTER_FOR(&glattice,ix) {
      for (mu=0;mu<4;mu++) {
         iy=iup(ix,mu);
         _suNg_times_suNg_dagger(v,*_4FIELD_AT(in,ix,mu),*_FIELD_AT(g,iy));
         _suNg_times_suNg(*_4FIELD_AT(out,ix,mu),*_FIELD_AT(g,ix),v);
      }
   }
}


extern struct {
  int c[3];
  int* path;
  int length;
  int** perm;
  int nperms;
  int nsteps;
} WL_path[256];
extern int WL_npaths;
extern int WL_max_nsteps;



int main(int argc,char *argv[])
{
  setup_process(&argc,&argv);
  
  logger_setlevel(0,10000); /* log all */
  if (PID!=0) { logger_disable(); }
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

  lprintf("MAIN",0,"Gauge group: SU(%d)\n",NG);
  lprintf("MAIN",0,"The lattice size is %dx%dx%dx%d\n",T,X,Y,Z);
  lprintf("MAIN",0,"The lattice global size is %dx%dx%dx%d\n",GLB_T,GLB_X,GLB_Y,GLB_Z);
  lprintf("MAIN",0,"The lattice borders are (%d,%d,%d,%d)\n",T_BORDER,X_BORDER,Y_BORDER,Z_BORDER);
  lprintf("MAIN",0,"RLXD [%d,%d]\n",rlx_var.rlxd_level,rlx_var.rlxd_seed+MPI_PID);
  rlxd_init(rlx_var.rlxd_level,rlx_var.rlxd_seed+MPI_PID); /* use unique MPI_PID to shift seeds */
  
  suNg_field* u[5];
  u[0]=alloc_gfield(&glattice);
  u[1]=alloc_gfield(&glattice);
  u[2]=alloc_gfield(&glattice);
  u[3]=alloc_gfield(&glattice);
  u[4]=alloc_gfield(&glattice);

  suNg_field* g=alloc_gtransf(&glattice);

  suNg* poly[2];
  poly[0]=amalloc(sizeof(suNg)*X*Y*Z,ALIGN);
  poly[1]=amalloc(sizeof(suNg)*X*Y*Z,ALIGN);
  
  random_u(u[0]);
  start_gf_sendrecv(u[0]);
  complete_gf_sendrecv(u[0]);
  
  random_g(g);
  start_gt_sendrecv(g);
  complete_gt_sendrecv(g);
  
  double HYP_weight[3]={1.,1.,1.};

  transform_u(u[1],u[0],g);
  start_gf_sendrecv(u[1]);
  complete_gf_sendrecv(u[1]);
  HYP_smearing(u[2],u[1],HYP_weight);
  
  HYP_smearing(u[3],u[0],HYP_weight);
  start_gf_sendrecv(u[3]);
  complete_gf_sendrecv(u[3]);
  transform_u(u[4],u[3],g);


  double err=0.;
  double dtmp;
  int t,x,y,z,i;
  for(t=0;t<T;t++) for(x=0;x<X;x++) for(y=0;y<Y;y++) for(z=0;z<Z;z++) {
    i=ipt(t,x,y,z);
    for(int mu=0; mu<4; mu++) {
      _suNg_sub_assign(*_4FIELD_AT(u[4],i,mu),*_4FIELD_AT(u[2],i,mu));
      _suNg_sqnorm(dtmp,*_4FIELD_AT(u[4],i,mu));
      if(dtmp>err) err=dtmp;
    }
  }


  
  lprintf("MAIN",0,"Checking gauge covariance of the HYP smearing\n");
  lprintf("MAIN",0,"\n");


#ifdef WITH_MPI
  int mpiret;

  dtmp=err/(2*NG*NG);
  mpiret=MPI_Allreduce(&dtmp,&err,1,MPI_DOUBLE,MPI_MAX,GLB_COMM);

  if (mpiret != MPI_SUCCESS) {
    char mesg[MPI_MAX_ERROR_STRING];
    int mesglen;
    MPI_Error_string(mpiret,mesg,&mesglen);
    lprintf("MPI",0,"ERROR: %s\n",mesg);
    error(1,1,"main [check_wilsonloops_3.c]","Cannot compute global maximum");
  }
#endif

  lprintf("MAIN",0,"Maximal normalized difference = %.2e\n",err);
  lprintf("MAIN",0,"(should be around 1*10^(-14) or so)\n");
  lprintf("MAIN",0,"\n");
  if(err<1.e-14)
    lprintf("MAIN",0,"check_wilsonloops_6 ... OK\n");
  else
    lprintf("MAIN",0,"check_wilsonloops_6 ... FAILED\n");

  WL_free();
  

  free_gfield(u[0]);
  free_gfield(u[1]);
  free_gfield(u[2]);
  free_gfield(u[3]);
  free_gfield(u[4]);
  free_gtransf(g);

  finalize_process();
  exit(0);
}
