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

#include "random.h"
#include "libhr_core.h"


void gaussian_suNg_vector(suNg_vector *v)
{
   gauss((double*)v, sizeof(suNg_vector)/sizeof(double));
}
#if !defined(WITH_QUATERNIONS) && !defined(GAUGE_SON)
/* generates a random SU(N) matrix via SU(2) rotations */
static void rotate(suNg_vector *pu1, suNg_vector *pu2, double s[4]) /* same as in cabmar */
{
	  hr_complex z1,z2;
	  hr_complex *cu1, *cu2;
  
	  cu1 = &((*pu1).c[0]);
	  cu2 = &((*pu2).c[0]);
						  
	  for (int i=0; i<NG; ++i) {
            z1=s[0]*creal(*cu1)-s[1]*cimag(*cu2)+s[2]*creal(*cu2)-s[3]*cimag(*cu1)+I*(s[0]*cimag(*cu1)+s[1]*creal(*cu2)+s[2]*cimag(*cu2)+s[3]*creal(*cu1));
            z2=s[0]*creal(*cu2)-s[1]*cimag(*cu1)-s[2]*creal(*cu1)+s[3]*cimag(*cu2)+I*(s[0]*cimag(*cu2)+s[1]*creal(*cu1)-s[2]*cimag(*cu1)-s[3]*creal(*cu2));
            (*cu1) = z1;
            (*cu2) = z2;
            ++cu1;
            ++cu2;
	  }
}

static void rotate_f(suNf_vector *pu1, suNf_vector *pu2, double s[4]) /* same as in cabmar */
{
	  double complex z1,z2;
	  double complex *cu1, *cu2;
  
	  cu1 = &((*pu1).c[0]);
	  cu2 = &((*pu2).c[0]);

    int dim = sizeof(suNf_vector)/sizeof(hr_complex);  
	  for (int i=0; i<dim; ++i) {
            z1=s[0]*creal(*cu1)-s[1]*cimag(*cu2)+s[2]*creal(*cu2)-s[3]*cimag(*cu1)+I*(s[0]*cimag(*cu1)+s[1]*creal(*cu2)+s[2]*cimag(*cu2)+s[3]*creal(*cu1));
            z2=s[0]*creal(*cu2)-s[1]*cimag(*cu1)-s[2]*creal(*cu1)+s[3]*cimag(*cu2)+I*(s[0]*cimag(*cu2)+s[1]*creal(*cu1)-s[2]*cimag(*cu1)-s[3]*creal(*cu2));
            (*cu1) = z1;
            (*cu2) = z2;
            ++cu1;
            ++cu2;
	  }
}
#endif 

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

void random_suNf(suNf *u) {
#ifdef WITH_QUATERNIONS
  random_su2(0.,u->c);
#else
  double s[4];
  suNf_vector *pu1=(suNf_vector*)(u);
	
  _suNf_unit(*u);
  
  for (int i=0; i<NG; ++i) {
    suNf_vector *pu2 = pu1 + 1;
    for (int j=i+1; j<NG; ++j) {
		  random_su2(0.0,s);
      rotate_f(pu1, pu2, s);
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

void random_suNf(suNf *u) {
  suNf tmp;
  double gr[NF*NF];
  int i;
  do {
    gauss(gr,NF*NF);
    for (i=0;i<NF*NF;i++){
      tmp.c[i]=gr[i];
    }
  } while (!project_to_suNg_real(u,&tmp));
}


#endif

