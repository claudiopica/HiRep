/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *   
* All rights reserved.                                                      * 
\***************************************************************************/

#include "update.h"
#include "dirac.h"
#include "inverters.h"
#include "linear_algebra.h"
#include "suN.h"
#include "random.h"
#include "memory.h"
#include "utils.h"
#include "global.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "logger.h"

static spinor_field *ev;

void find_spec_H2(double *max, double *min) {
  /* EVA parameters */
  const int nevt=5; /* use 5-dim space */
  const int nev=1; /* require only the smallest to be accurate */
  const int kmax=200; /* max degree of polynomial */
  const int maxiter=20; /* max number of subiterations */
  static double *d1;
  const double omega1=1.e-8; /* absolute precision */
  const double omega2=1.e-1; /* relative precision */
  int status,ie;
  /* END of EVA parameters */
  int MVM=0; /* counter for matrix-vector multiplications */


  d1=malloc(sizeof(*d1)*nevt);
#ifdef UPDATE_EO
  ev=alloc_spinor_field_f(nevt,&glat_even);
  MVM+=max_H(&H2, &glat_even, max);
#else
  ev=alloc_spinor_field_f(nevt,&glattice);
  MVM+=max_H(&H2, &glattice, max);
#endif

  ie=eva(nev,nevt,0,kmax,maxiter,*max,omega1,omega2,&H2,ev,d1,&status);
  MVM+=status;
  while (ie!=0) { /* if failed restart EVA */
    ie=eva(nev,nevt,2,kmax,maxiter,*max,omega1,omega2,&H2,ev,d1,&status);
    MVM+=status;
  }

  *min=d1[0]*(1-omega2)-omega1;

  lprintf("SPECLIMITS",0,"Range = [%1.8e,%1.8e] [MVM = %d]\n",*min,*max,MVM);

  free(d1);
  free_spinor_field_f(ev);

  return;
}


