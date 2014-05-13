/******************************************************************************
*
* Test of modules
*
******************************************************************************/

#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "io.h"
#include "random.h"
#include "error.h"
#include "geometry.h"
#include "memory.h"
#include "statistics.h"
#include "update.h"
#include "global.h"
#include "observables.h"
#include "suN.h"
#include "suN_types.h"
#include "dirac.h"
#include "linear_algebra.h"
#include "inverters.h"
#include "representation.h"
#include "utils.h"
#include "logger.h"
#include "communications.h"

int nhb,nor,nit,nth,nms,level,seed;
double beta;

static double hmass=0.1;


void D(spinor_field *out, spinor_field *in){
   Dphi(hmass,out,in);
}

void H(spinor_field *out, spinor_field *in){
   g5Dphi(hmass,out,in);
}

static spinor_field *tmp;
void M(spinor_field *out, spinor_field *in){
   g5Dphi(-hmass,tmp,in); 
   g5Dphi(-hmass,out,tmp);
}


int main(int argc,char *argv[])
{
  char pame[256];
  int i;
  double tau;
  spinor_field *s1, *s2;
  spinor_field *res;

  mshift_par par;
  MINRES_par MINRESpar2;
  
  int cgiters;
  /* setup process id and communications */
  setup_process(&argc,&argv);
  
  /* logger setup */

  logger_setlevel(0,100); /* log all */
  if (PID!=0) { 
    logger_disable();}
  else{
    sprintf(pame,">out_%d",PID); logger_stdout(pame);
    sprintf(pame,"err_%d",PID); freopen(pame,"w",stderr);
  }
   
  lprintf("MAIN",0,"PId =  %d [world_size: %d]\n\n",PID,WORLD_SIZE); 
  
  /* read input file */
  read_input(glb_var.read,"test_input");
  
  /* setup communication geometry */
  if (geometry_init() == 1) {
    finalize_process();
    return 0;
  }
  
   
   lprintf("MAIN",0,"Gauge group: SU(%d)\n",NG);
   lprintf("MAIN",0,"Fermion representation: " REPR_NAME " [dim=%d]\n",NF);
   lprintf("MAIN",0,"global size is %dx%dx%dx%d\n",GLB_T,GLB_X,GLB_Y,GLB_Z);
   lprintf("MAIN",0,"proc grid is %dx%dx%dx%d\n",NP_T,NP_X,NP_Y,NP_Z);
   
   /* setup lattice geometry */
   geometry_mpi_eo();
   /* test_geometry_mpi_eo(); */
    
    
    /* setup random numbers */
    read_input(rlx_var.read,"test_input");
    lprintf("MAIN",0,"RLXD [%d,%d]\n",rlx_var.rlxd_level,rlx_var.rlxd_seed+MPI_PID);
    rlxd_init(rlx_var.rlxd_level,rlx_var.rlxd_seed+MPI_PID); /* use unique MPI_PID to shift seeds */

   
   u_gauge=alloc_gfield(&glattice);
#ifndef REPR_FUNDAMENTAL
   u_gauge_f=alloc_gfield_f(&glattice);
#endif
  represent_gauge_field();

  lprintf("MAIN",0,"Generating a random gauge field... ");
  fflush(stdout);
  random_u(u_gauge);
  
  start_gf_sendrecv(u_gauge);
  complete_gf_sendrecv(u_gauge);

  lprintf("MAIN",0,"done.\n");
  represent_gauge_field();

  par.n = 6;
  par.shift=(double*)malloc(sizeof(double)*(par.n));
  par.err2=1.e-28;
  par.max_iter=0;
  res=alloc_spinor_field_f(par.n+3,&glattice);
  s1=res+par.n;
  s2=s1+1;
  tmp=s2+1;

  par.shift[0]=+0.1;
  par.shift[1]=-0.21;
  par.shift[2]=+0.05;
  par.shift[3]=-0.01;
  par.shift[4]=-0.15;
  par.shift[5]=-0.05;

  gaussian_spinor_field(s1);

  /* TEST MINRES_M */

  lprintf("MR TEST",0,"\n");
  lprintf("MR TEST",0,"Testing MINRES multishift\n");
  lprintf("MR TEST",0,"-------------------------\n");

  cgiters=MINRES_mshift(&par, &H, s1, res);
  lprintf("MR TEST",0,"Converged in %d iterations\n",cgiters);

  for(i=0;i<par.n;++i){
    H(s2,&res[i]);
    spinor_field_mul_add_assign_f(s2,-par.shift[i],&res[i]);
    spinor_field_sub_assign_f(s2,s1);
    tau=spinor_field_sqnorm_f(s2)/spinor_field_sqnorm_f(s1);
    lprintf("MR TEST",0,"test MINRES[%d] = %e (req. %e)\n",i,tau,par.err2);
  }

  /* TEST MINRES_M */

  lprintf("MR TEST",0,"\n");
  lprintf("MR TEST",0,"Testing MINRES \n");
  lprintf("MR TEST",0,"-------------- \n");

  MINRESpar2.err2=par.err2;
  MINRESpar2.max_iter=0;

  cgiters=MINRES(&MINRESpar2, &H, s1, &res[0],0);
  for(i=1;i<par.n;++i){
    hmass=0.1-par.shift[i-1];
    cgiters+=MINRES(&MINRESpar2, &H, s1, &res[i],&res[i-1]);
  }
  lprintf("MR TEST",0,"Converged in %d iterations\n",cgiters);

  hmass=0.1;
  for(i=0;i<par.n;++i){
    if(i!=0)
      hmass=0.1-par.shift[i-1];
    H(s2,&res[i]);
    spinor_field_sub_f(s2,s2,s1);
    tau=spinor_field_sqnorm_f(s2)/spinor_field_sqnorm_f(s1);
    lprintf("MR TEST",0,"test MINRES[%d] = %e (req. %e)\n",i,tau,par.err2);
  }

  free_spinor_field_f(res);
  free(par.shift);

  finalize_process();

  return 0;
}
