/*******************************************************************************
*
* Check that the Hamiltonian gauge sets the temporal links to the identity
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
  
  WL_initialize();
  
  WL_Hamiltonian_gauge(u_gauge,u_gauge);
  start_gf_sendrecv(u_gauge);
  complete_gf_sendrecv(u_gauge);
  
  double dtmp;

  lprintf("MAIN",0,"Checking that the Hamiltonian gauge sets the temporal links to the identity\n");
  lprintf("MAIN",0,"\n");

  
  double err=0.;
  int t,x,y,z,i;
  if(COORD[0]!=NP_T-1) {
    for(t=0;t<T;t++) for(x=0;x<X;x++) for(y=0;y<Y;y++) for(z=0;z<Z;z++) {
      i=ipt(t,x,y,z);
      _suNg_sqnorm_m1(dtmp,*_4FIELD_AT(u_gauge,i,0));
      if(dtmp>err) err=dtmp;
    }
  } else {
    for(t=0;t<T-1;t++) for(x=0;x<X;x++) for(y=0;y<Y;y++) for(z=0;z<Z;z++) {
      i=ipt(t,x,y,z);
      _suNg_sqnorm_m1(dtmp,*_4FIELD_AT(u_gauge,i,0));
      if(dtmp>err) err=dtmp;
    }
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
    error(1,1,"main [check_wilsonloops_1.c]","Cannot compute global maximum");
  }
#endif

  lprintf("MAIN",0,"CID=%d Maximal normalized difference = %.2e\n",CID,err);
  lprintf("MAIN",0,"(should be around 1*10^(-15) or so)\n");
  lprintf("MAIN",0,"\n");
  if(err<1.e-15)
    lprintf("MAIN",0,"CID=%d check_wilsonloops_1 ... OK\n",CID);
  else
    lprintf("MAIN",0,"CID=%d check_wilsonloops_1 ... FAILED\n",CID);
  
  WL_free();
  
  free_gfield(u_gauge);
  
  finalize_process();
  exit(0);
}
