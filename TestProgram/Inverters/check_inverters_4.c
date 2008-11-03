
/*******************************************************************************
*
* File check9.c
*
* Check of the program eva (random field)
*
* Author: Luigi Del Debbio
*
*******************************************************************************/

#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "suN.h"
#include "random.h"
#include "error.h"
#include "geometry.h"
#include "memory.h"
#include "linear_algebra.h"
#include "update.h"
#include "inverters.h"
#include "dirac.h"
#include "representation.h"
#include "global.h"
#include "logger.h"
#include "io.h"
#include "communications.h"

static int iw;
static double hmass=-7.94871867e-01f;
static spinor_field *ws,*ev;
static double EPSILON=1.e-12;




static double normalize(spinor_field *ps)
{
  double r;

  r=spinor_field_sqnorm_f(ps);
  r=sqrt(r);
  error(r<EPSILON,1,"normalize [check9.c]","vector has vanishing norm");

  r=1.0/r;
  spinor_field_mul_f(ps,r,ps);

  return (double)(1.0/r);
}

static void Op1(spinor_field *out,spinor_field *in)
{
  g5Dphi(hmass,&ev[iw],in);
  g5Dphi(hmass,out,&ev[iw]);
}


static double power(int nit,spinor_operator Op,spinor_field *ws)
{
  int i;
  double ubnd;

  gaussian_spinor_field(&ws[0]);
  normalize(&ws[0]);
  Op(&ws[1],&ws[0]);
  Op(&ws[0],&ws[1]);
  ubnd=normalize(&ws[0]);   

  for (i=1;i<nit;i++)
    {
      Op(&ws[1],&ws[0]);
      Op(&ws[0],&ws[1]);
      ubnd=normalize(&ws[0]);
    }
  return (double)sqrt((double)(ubnd));
}


int main(int argc,char *argv[])
{
  int i;
  int nev,nevt,ie,status;
  double omega1,omega2,d1[6],res,ubnd;
  complex z;
  char pame[256];

  /* setup process id and communications */
  setup_process(&argc,&argv);
  
  /* logger setup */
  logger_setlevel(0,10000); /* log all */
  logger_map("DEBUG","debug");
#ifdef WITH_MPI
  sprintf(pame,">out_%d",PID); logger_stdout(pame);
  sprintf(pame,"err_%d",PID); freopen(pame,"w",stderr);
#endif
  
  lprintf("MAIN",0,"PId =  %d [world_size: %d]\n\n",PID,WORLD_SIZE); 
  
  /* read input file */
  read_input("test_input");
  rlxd_init(1,12345+PID);
  
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
   
   u_gauge=alloc_gfield(&glattice);
#ifndef REPR_FUNDAMENTAL
   u_gauge_f=alloc_gfield_f(&glattice);
#endif
  lprintf("MAIN",0,"Generating a random gauge field... ");fflush(stdout);
  random_u(u_gauge);
  
  start_gf_sendrecv(u_gauge);
  complete_gf_sendrecv(u_gauge);
  lprintf("MAIN",0,"done.\n");
  represent_gauge_field();


  lprintf("MAIN",0,"Diagonalization of Q^2 (random fields)\n");
  lprintf("MAIN",0,"--------------------------------------\n\n");

 
  
  ws=alloc_spinor_field_f(2,&glattice);
  ev=alloc_spinor_field_f(7,&glattice);


   

  iw=4;
  nev=3;
  nevt=4;
  ubnd=1.05f*power(30,Op1,ws);
  lprintf("MAIN",0,"test-ubnd: %f\n",ubnd);
  omega1=1.0e-6f;
  omega2=1.0e-3f;

  lprintf("MAIN",0,"Accuracy parameters: omega1=%.1e, omega2=%.1e\n\n",
	 omega1,omega2);

  ie=eva(nev,nevt,0,100,20,ubnd,omega1,omega2,Op1,ws,ev,d1,&status);

  lprintf("MAIN",0,"\nEigenvalues of Q^2 (status = %d, ie = %d):\n\n",
	 status,ie);

  for (i=0;i<nevt;i++)
    {
      Op1(&ws[0],&ev[i]);
  
      z.re=-(double)d1[i];
      z.im=(double)0.0f;
      spinor_field_mulc_add_assign_f(&ws[0],z,&ev[i]);
      res=spinor_field_sqnorm_f(&ws[0]);
      res=(double)(sqrt((double)(res)));
      
      if (i==nev)
	lprintf("MAIN",0,"\n");
      
      lprintf("MAIN",0,"d[%d] = % .3e, acc = %.1e\n",i,d1[i],res);
    }
  finalize_process();


  exit(0);
}

