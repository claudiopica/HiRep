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

  WL_initialize();

  WL_Hamiltonian_gauge(u_gauge,u_gauge);
  start_gf_sendrecv(u_gauge);
  complete_gf_sendrecv(u_gauge);

  double dtmp;

  lprintf("MAIN",0,"Checking that the Hamiltonian gauge sets the temporal links to the identity\n");
  lprintf("MAIN",0,"\n");


  double err=0.;
  int t,x,y,z,i;
  if(COORD[0]==NP_T-1){
    for(t=0;t<T-1;t++) for(x=0;x<X;x++) for(y=0;y<Y;y++) for(z=0;z<Z;z++) {
            i=ipt(t,x,y,z);
            _suNg_sqnorm_m1(dtmp,*_4FIELD_AT(u_gauge,i,0));
            if(dtmp>err) err=dtmp;
          }
  } else {
    for(t=0;t<T;t++) for(x=0;x<X;x++) for(y=0;y<Y;y++) for(z=0;z<Z;z++) {
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
    lprintf("MAIN",0,"(should be around 1*10^(-8) or so)\n");
    lprintf("MAIN",0,"\n");
    if(err<1.e-8)
      lprintf("MAIN",0,"CID=%d check_wilsonloops_1 ... OK\n",CID);
    else
    {
      return_value +=1 ;
      lprintf("MAIN",0,"CID=%d check_wilsonloops_1 ... FAILED\n",CID);
    }

    WL_free();

    free_gfield(u_gauge);

    finalize_process();
    return return_value;
  }
