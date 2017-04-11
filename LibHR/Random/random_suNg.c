/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *   
* All rights reserved.                                                      * 
\***************************************************************************/

/*******************************************************************************
*
* File random_suNg.c
*
* Generation of uniformly distributed SU(Ng) matrices
*
*******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "suN.h"
#include "random.h"
#include "utils.h"

extern void random_su2(double rho,double s[]);

void gaussian_suNg_vector(suNg_vector *v)
{
   gauss((double*)v, sizeof(suNg_vector)/sizeof(double));
}

/* generates a random SU(N) matrix via SU(2) rotations */
static void rotate(suNg_vector *pu1, suNg_vector *pu2, double s[4]) /* same as in cabmar */
{
	  complex z1,z2;
	  complex *cu1, *cu2;
  
	  cu1 = &((*pu1).c[0]);
	  cu2 = &((*pu2).c[0]);
						  
	  for (int i=0; i<NG; ++i) {
      z1.re=s[0]*(*cu1).re-s[1]*(*cu2).im+s[2]*(*cu2).re-s[3]*(*cu1).im;
      z1.im=s[0]*(*cu1).im+s[1]*(*cu2).re+s[2]*(*cu2).im+s[3]*(*cu1).re;
      z2.re=s[0]*(*cu2).re-s[1]*(*cu1).im-s[2]*(*cu1).re+s[3]*(*cu2).im;
      z2.im=s[0]*(*cu2).im+s[1]*(*cu1).re-s[2]*(*cu1).im-s[3]*(*cu2).re;
      (*cu1) = z1;
      (*cu2) = z2;
      ++cu1;
      ++cu2;
	  }
}

#ifndef GAUGE_SON
void random_suNg(suNg *u) {
#ifdef WITH_QUATERNIONS
  random_su2(0.,u->c);
#else
  double s[4];
  suNg_vector *pu1=(suNg_vector*)(u);
	
  _suNg_unit(*u);
  
  for (int i=0; i<NG; ++i) {
    suNg_vector *pu2 = pu1 + 1;
    for (int j=i+1; j<NG; ++j) {
		  random_su2(0.0,s);
      rotate(pu1, pu2, s);
      ++pu2; 
    } 
	  ++pu1; 
  }
#endif //WITH_QUATERNIONS
}
#else
void random_suNg(suNg *u) {
  suNg tmp;
  double gr[NG*NG];
  int i;
  do {
    gauss(gr,NG*NG);
    for (i=0;i<NG*NG;i++){
      tmp.c[i]=gr[i];
    }
  } while (!project_to_suNg_real(u,&tmp));
}

#endif

