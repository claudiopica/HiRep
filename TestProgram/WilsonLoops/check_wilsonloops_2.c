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
#include "setup.h"



int main(int argc,char *argv[])
{
  int return_value = 0;
  logger_map("DEBUG","debug");
  setup_process(&argc,&argv);
  setup_gauge_fields();


  random_u(u_gauge);
  start_gf_sendrecv(u_gauge);
  complete_gf_sendrecv(u_gauge);

  double complex plaq[T*X*Y*Z][6];
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


  double complex ctmp;
  double dtmp;
  double err=0.;

  j=0;
  for(t=0;t<T;t++) for(x=0;x<X;x++) for(y=0;y<Y;y++) for(z=0;z<Z;z++) {
    i=ipt(t,x,y,z);
    k=0;
    for(mu=0;mu<4;mu++) for(nu=0;nu<mu;nu++) {
      cplaq(&ctmp,i,mu,nu);
      dtmp=(plaq[j][k]-ctmp)*conj(plaq[j][k]-ctmp);
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
  {
    return_value +=1 ;
    lprintf("MAIN",0,"check_wilsonloops_2 ... FAILED\n");
  }
  WL_free();

  free_gfield(u_gauge);

  finalize_process();
  return return_value;
}
