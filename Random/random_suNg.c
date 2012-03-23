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

/* static variables for random_suNg */
static double s[4];
static suNg_vector *pu1,*pu2;


void random_su2(double rho,double s[]);


void random_suNg(suNg *u) {
  suNg tmp;
  double gr[NG*NG];
  int i;
  gauss(gr,NG*NG);
  for (i=0;i<NG*NG;i++){
    tmp.c[i].re=gr[i];
    tmp.c[i].im=0;
  }
  project_to_suNg_real(u,&tmp);
}

void random_suNg_unit_vector(suNg_vector *v)
{
   double norm=0.f,fact;

   while ((1.0+norm)==1.0)   {
     gauss((double*)v,(sizeof(suNg_vector)/sizeof(double)));
     _vector_prod_re_g(norm,*v,*v);
     norm=sqrt(norm);
   }

   fact=1.0/norm;
   _vector_mul_g(*v,fact,*v);
}

void gaussian_suNg_vector(suNg_vector *v)
{
   gauss((double*)v,(sizeof(suNg_vector)/sizeof(double)));
}


/* generates a random SU(N) matrix via SU(2) rotations */
static void rotate(void) /* same as in cabmar */
{
	  int i;
	  complex z1,z2;
	  complex *cu1, *cu2;
				  
	  cu1 = &((*pu1).c[0]);
	  cu2 = &((*pu2).c[0]);
						  
	  for (i=0; i<NG; ++i) {
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

void random_suNg_suN(suNg *u) {
  int i,j;

	_suNg_unit(*u);
  pu1=(suNg_vector*)(u);
				  
  for (i=0; i<NG; ++i) {
    pu2 = pu1 + 1;
    for (j=i+1; j<NG; ++j) {
		  random_su2(0.0,s);
      rotate();
      ++pu2; 
    } 
	  ++pu1; 
  }

}

/* this generates a U(N) matrix but not necessarely in SU(N) */
void random_suNg_old(suNg *u)
{
  int i, j;
  double norm=0.f,fact;
  suNg_vector *v1,*v2;
  complex z;
  
  v1=(suNg_vector*)(u);
  v2=v1+1;
  
  random_suNg_unit_vector(v1);
  for (i=1; i<NG; ++i) {
    while ((1.0+norm)==1.0) {
      random_suNg_unit_vector(v2);
      for (j=i; j>0; --j) {
	_vector_prod_re_g(z.re,*v1, *v2);
	_vector_prod_im_g(z.im,*v1, *v2);
	_vector_project_g(*v2, z, *v1); 
	++v1;
      }
      _vector_prod_re_g(norm,*v2,*v2);
      norm=sqrt(norm);
    }        
    fact=1.0/norm;
    _vector_mul_g(*v2,fact,*v2); /* normalize v2 */
    norm=0.;
    ++v2;
    v1=(suNg_vector*)(u);
  }
  /* TEST */
  /*
  v1=(suNg_vector*)(u);
  v2=v1;
  for (i=0; i<NG; ++i) {
    norm = _vector_prod_re_g(*v2, *v2);
    printf("norm2 %d= %4.5f\n", i, norm);
    v2++;
  }
  for (i=0; i<NG-1; ++i) {
    for (j=1+i; j<NG; ++j){
      v1=(suNg_vector*)(u)+i;
      v2=(suNg_vector*)(u)+j;
      z.re = _vector_prod_re_g(*v1, *v2);
      z.im = _vector_prod_im_g(*v1, *v2);
      printf("prod %d %d= (%4.5f,%4.5f)\n", i+1,j+1, z.re,z.im);
    }
  }
  */
}

