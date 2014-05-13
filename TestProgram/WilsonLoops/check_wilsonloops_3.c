/*******************************************************************************
*
* Check the Polyakov loop under Hamiltonian gauge
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



void polyakov_test(suNg* poly, suNg_field* gf) {
  suNg_field* sh=alloc_gtransf(&glattice);
  suNg tmp;
  
  for(int x=0;x<X;x++) for(int y=0;y<Y;y++) for(int z=0;z<Z;z++) {
    _suNg_unit(poly[_WL_3VOL_INDEX(x,y,z)]);
  }
  
  for(int t=0;t<T;t++) for(int x=0;x<X;x++) for(int y=0;y<Y;y++) for(int z=0;z<Z;z++) {
    *_FIELD_AT(sh,ipt(t,x,y,z)) = *_4FIELD_AT(gf,ipt(t,x,y,z),0);
  }

  start_gt_sendrecv(sh);
  complete_gt_sendrecv(sh);

  for(int s=0;s<GLB_T;s++) {
    for(int x=0;x<X;x++) for(int y=0;y<Y;y++) for(int z=0;z<Z;z++) {
      _suNg_times_suNg(tmp,poly[_WL_3VOL_INDEX(x,y,z)],*_FIELD_AT(sh,ipt(0,x,y,z)));
      poly[_WL_3VOL_INDEX(x,y,z)]=tmp;
    }
    
    for(int t=0;t<T;t++) for(int x=0;x<X;x++) for(int y=0;y<Y;y++) for(int z=0;z<Z;z++) {
      int i=ipt(t,x,y,z);
      *_FIELD_AT(sh,i) = *_FIELD_AT(sh,iup(i,0));
    }
    start_gt_sendrecv(sh);
    complete_gt_sendrecv(sh);
  }
  
  free_gtransf(sh);
}




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
  
  u_gauge=alloc_gfield(&glattice);
  random_u(u_gauge);
  start_gf_sendrecv(u_gauge);
  complete_gf_sendrecv(u_gauge);
  
  suNg* poly[2];
  poly[0]=amalloc(sizeof(suNg)*X*Y*Z,ALIGN);
  poly[1]=amalloc(sizeof(suNg)*X*Y*Z,ALIGN);
  
  polyakov_test(poly[0],u_gauge);
  
  WL_initialize();
  
  WL_Hamiltonian_gauge(u_gauge,u_gauge);

  WL_broadcast_polyakov(poly[1],u_gauge);
  
  lprintf("MAIN",0,"Checking that Polyakov loop do not change in the Hamiltonian gauge\n");
  lprintf("MAIN",0,"\n");

 
  double dtmp;
  double err=0.;
  int x,y,z;

  for(x=0;x<X;x++) for(y=0;y<Y;y++) for(z=0;z<Z;z++) {
    _suNg_sub_assign(poly[0][_WL_3VOL_INDEX(x,y,z)],poly[1][_WL_3VOL_INDEX(x,y,z)]);
    _suNg_sqnorm(dtmp,poly[0][_WL_3VOL_INDEX(x,y,z)]);
    if(dtmp>err) err=dtmp;
  }
  
  err=sqrt(err/(2*NG*NG*GLB_T));

#ifdef WITH_MPI
  int mpiret;

  dtmp=err;
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
    lprintf("MAIN",0,"check_wilsonloops_3 ... OK\n");
  else
    lprintf("MAIN",0,"check_wilsonloops_3 ... FAILED\n");
  
  WL_free();
  
  free_gfield(u_gauge);
  afree(poly[0]);
  afree(poly[1]);
  
  finalize_process();
  exit(0);
}
