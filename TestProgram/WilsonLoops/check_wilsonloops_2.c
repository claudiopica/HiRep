/*******************************************************************************
*
* Check plaquettes do not change in the Hamiltonian gauge
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
  
  complex plaq[T*X*Y*Z][6];
  
  int t,x,y,z,i,j,k,mu,nu;
  
  j=0;
  for(t=0;t<T;t++) for(x=0;x<X;x++) for(y=0;y<Y;y++) for(z=0;z<Z;z++) {
    i=ipt(t,x,y,z);
    k=0;
    for(mu=0;mu<4;mu++) for(nu=0;nu<mu;nu++) {
      cplaq(&(plaq[j][k]),i,mu,nu);
      k++;
    }
    j++;
  }
  
  
  WL_initialize();
  
  WL_Hamiltonian_gauge(u_gauge,u_gauge);
 
  lprintf("MAIN",0,"Checking that plaquettes do not change in the Hamiltonian gauge\n");
  lprintf("MAIN",0,"\n");

 
  complex ctmp;
  double dtmp;
  double err=0.;

  j=0;
  for(t=0;t<T;t++) for(x=0;x<X;x++) for(y=0;y<Y;y++) for(z=0;z<Z;z++) {
    i=ipt(t,x,y,z);
    k=0;
    for(mu=0;mu<4;mu++) for(nu=0;nu<mu;nu++) {
      cplaq(&ctmp,i,mu,nu);
      dtmp=(plaq[j][k].re-ctmp.re)*(plaq[j][k].re-ctmp.re)+(plaq[j][k].im-ctmp.im)*(plaq[j][k].im-ctmp.im);
      if(dtmp>err) err=dtmp;
/*      printf("%d\t%d\t%d\t%d\t##\t%d\t%d\t##\t%.2e\n",t,x,y,z,mu,nu,dtmp);*/
      k++;
    }
    j++;
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
    error(1,1,"main [check_wilsonloops_2.c]","Cannot compute global maximum");
  }
#endif

  lprintf("MAIN",0,"Maximal normalized difference = %.2e\n",err);
  lprintf("MAIN",0,"(should be around 1*10^(-14) or so)\n");
  lprintf("MAIN",0,"\n");
  if(err<1.e-14)
    lprintf("MAIN",0,"check_wilsonloops_2 ... OK\n");
  else
    lprintf("MAIN",0,"check_wilsonloops_2 ... FAILED\n");
  
  WL_free();
  
  free_gfield(u_gauge);
  
  finalize_process();
  exit(0);
}
