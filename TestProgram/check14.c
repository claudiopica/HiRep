
/******************************************************************************
*
* File check14.c
*
* Check of dirac_eva
*
* Author: Agostino Patella
*
******************************************************************************/

#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
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
#include "observables.h"


#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288419716939937510
#endif


static double mass = .1;


static void H(suNf_spinor *out, suNf_spinor *in){
	g5Dphi(mass,out,in);
}


int main(int argc,char *argv[]) {
  int i, j;
  FILE *log=NULL;   
  double norm;
  complex prod;
  suNf_spinor *test;
  suNf_spinor **ev;
  double *d;
  int nevt, nev;

  int status;
  /* On exit this variable reports the number of times the */
  /* operator was applied */

  int kmax;
  /* Maximal degree of the Chebyshev polynomials used for the  */
  /* acceleration of the algorithm */

  int imax;
  /* Maximal number of subspace iterations */

  float omega1, omega2;
  /* The first eigenvalues are obtained to an absolute precision omega1 or */
  /* a relative precision omega2 (whichever is reached first). */


  nev = 100;
  nevt = 100;
  kmax = 50;
  imax = nevt*5;
  omega1 = 1e-6;
  omega2 = 1e-6;


  log=freopen("check14.log","w",stdout);
  printf("\n");
  printf("Diagonalization of the g5D operator.\n");
  printf("---------------------------------------\n\n");
  printf("The lattice size is %dx%d^3\n",T,L);
  printf("size of the gluon rep: %d, size of the fermion rep: %d\n",NG,NF);
  printf("mass of the fermion: %f\n",mass);

  rlxd_init(1,12345);

  logger_setlevel(0,1000);

  geometry_eo_lexi();
  u_gauge=alloc_gfield();
  #ifndef REPR_FUNDAMENTAL
  u_gauge_f=alloc_gfield_f();
  #endif

  printf("Generating a random gauge field... ");
  fflush(stdout);
  random_u();
  printf("done.\n");
  represent_gauge_field();

  set_spinor_len(VOLUME);

  test = alloc_spinor_field_f(1);

  d = (double*)malloc(sizeof(double)*nevt);
  ev = (suNf_spinor**)malloc(sizeof(suNf_spinor*)*nevt);
  for(i = 0; i < nevt; i++)
    ev[i] = alloc_spinor_field_f(1);


  dirac_eva_onemass(nev,nevt,kmax,imax,omega1,omega2,mass,ev,d,&status);

  printf("\n");

  for(i = 0; i < nevt; i++) {
    H(test,ev[i]);
    spinor_field_mul_add_assign_f(test,-d[i],ev[i]);
    norm=spinor_field_sqnorm_f(test);
    printf("Eigenvector test [%d,%e,%e] = %e\n",i,d[i],d[i]*d[i],norm);
  }

  printf("\n");

  norm = 0.;
  for(i = 0; i < nevt; i++)
  for(j = 0; j < nevt; j++) {
    prod.re = spinor_field_prod_re_f(ev[i],ev[j]);
    prod.im = spinor_field_prod_im_f(ev[i],ev[j]);
    if(i == j) norm += (prod.re-1.)*(prod.re-1.);
    else norm += prod.re*prod.re;
    norm += prod.im*prod.im;
  }
  norm = sqrt(norm);
  printf("Norm test = %e\n",norm);

  printf("\n");
  fclose(log);

  exit(0);
}

