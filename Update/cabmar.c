/*******************************************************************************
*
* File cabmar.c
*
* Cabibbo-Marinari rotations
*
*******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "suN.h"
#include "update.h"


static double s[4],w[4];
static suNg_vector *pu1,*pu2,*pv1,*pv2;

static void rotate(void)
{
  int i;
  complex z1,z2;
  complex *cu1, *cu2;
  
  cu1 = &((*pu1).c1);
  cu2 = &((*pu2).c1);
  
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

static void wmatrix(void) 
{
  w[0] = _vector_prod_re_g(*pu1, *pv1)+_vector_prod_re_g(*pu2,*pv2);
  w[1] = _vector_prod_im_g(*pu1, *pv2)+_vector_prod_im_g(*pu2,*pv1);
  w[2] = _vector_prod_re_g(*pu2, *pv1)-_vector_prod_re_g(*pu1,*pv2);
  w[3] = _vector_prod_im_g(*pu1, *pv1)-_vector_prod_im_g(*pu2,*pv2);
}

void cabmar(double beta,suNg *u,suNg *v,int type)
{
  int i,j;
  double b,bsq,wsq,rho,fact,r[4];
  
  const double invng = 1. / (double) NG;
  
  pu1=(suNg_vector*)(u);
  pv1=(suNg_vector*)(v);
  
  b=invng*beta; 
  bsq=b*b;
  
  for (i=0; i<NG; ++i) {
    pu2 = pu1 + 1;
    pv2 = pv1 + 1;
    for (j=i+1; j<NG; ++j) {
      wmatrix();
      wsq=w[0]*w[0]+w[1]*w[1]+w[2]*w[2]+w[3]*w[3];
      
      if ((bsq*wsq)>1.0e-16f) {
	if (type==1) {
	  fact=(w[0]+w[0])/wsq;
	  s[0]=fact*w[0]-1.0f;
	  s[1]=fact*w[1];
	  s[2]=fact*w[2];
	  s[3]=fact*w[3];         
	} else {
	  fact=(float)(sqrt((double)(wsq)));
	  rho=b*fact;
	  
	  random_su2(rho,r);
	  
	  fact=1.0f/fact;            
	  s[0]=fact*(r[0]*w[0]-r[1]*w[1]-r[2]*w[2]-r[3]*w[3]);
	  s[1]=fact*(r[1]*w[0]+r[0]*w[1]-r[2]*w[3]+r[3]*w[2]);
	  s[2]=fact*(r[2]*w[0]+r[0]*w[2]-r[3]*w[1]+r[1]*w[3]);
	  s[3]=fact*(r[3]*w[0]+r[0]*w[3]-r[1]*w[2]+r[2]*w[1]);         
	}      
      } else
	random_su2(0.0f,s);
      
      rotate();
      
      ++pu2; ++pv2;
    }
    ++pu1; ++pv1;
   }
}

