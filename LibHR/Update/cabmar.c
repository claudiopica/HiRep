/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *   
* All rights reserved.                                                      * 
\***************************************************************************/

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


#ifndef WITH_QUATERNIONS

static inline void rotate(suNg_vector *pu1, suNg_vector *pu2, double *s)
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

static inline void wmatrix(suNg_vector *pu1, suNg_vector *pu2, suNg_vector *pv1, suNg_vector *pv2, double *w)
{
	double prod1,prod2;
  _vector_prod_re_g(prod1,*pu1,*pv1);
	_vector_prod_re_g(prod2,*pu2,*pv2);
	w[0] = prod1+prod2;
  _vector_prod_im_g(prod1,*pu1,*pv2);
	_vector_prod_im_g(prod2,*pu2,*pv1);
	w[1] = prod1+prod2;
  _vector_prod_re_g(prod1,*pu2,*pv1);
	_vector_prod_re_g(prod2,*pu1,*pv2);
	w[2] = prod1-prod2;
  _vector_prod_im_g(prod1,*pu1,*pv1);
	_vector_prod_im_g(prod2,*pu2,*pv2);
	w[3] = prod1-prod2;
}

#endif

void cabmar(double beta,suNg *u,suNg *v,int type)
{
#ifdef WITH_QUATERNIONS
  double wsq,rho,fact;
	suNg s,w,r;
  
	/*w=u*v^+ */
	_suNg_times_suNg_dagger(w,*u,*v);
	
	/*wsq=w[0]*w[0]+w[1]*w[1]+w[2]*w[2]+w[3]*w[3];*/
	_suNg_sqnorm(wsq,w); //sqnorm gives twice what we need
  wsq*=0.5;
  
	if ((beta*beta*wsq)>1.0e-28) {
		if (type==1) {
			_suNg_trace_re(fact,w);
			fact/=wsq;
			_suNg_mul(s,fact,*v);
			/* _suNg_sub(*u,s,*u); */
      _suNg_minus(*u,*u);
      _suNg_add_assign(*u,s);
		} else {
			fact=sqrt(wsq);
			rho=beta*fact;
      
			random_su2(rho,&r.c[0]);
      
			fact=1./fact;
      
      _suNg_times_suNg(*u,r,*v);
			_suNg_mul(*u,fact,*u);
		}
	} else random_su2(0.0,&(u->c[0]));
	
  
#else
  int i,j;
  double b,bsq,wsq,rho,fact,r[4],w[4],s[4];
  
  const double invng = 1. / (double) NG;
  
  suNg_vector *pu1=(suNg_vector*)(u);
  suNg_vector *pv1=(suNg_vector*)(v);
  
  b=invng*beta;
  bsq=b*b;
  
  for (i=0; i<NG; ++i) {
    suNg_vector *pu2 = pu1 + 1;
    suNg_vector *pv2 = pv1 + 1;
    for (j=i+1; j<NG; ++j) {
      wmatrix(pu1, pu2, pv1, pv2, w);
      wsq=w[0]*w[0]+w[1]*w[1]+w[2]*w[2]+w[3]*w[3];
      
      if ((bsq*wsq)>1.0e-28) {
        if (type==1) {
          fact=(w[0]+w[0])/wsq;
          s[0]=fact*w[0]-1.0;
          s[1]=fact*w[1];
          s[2]=fact*w[2];
          s[3]=fact*w[3];
        } else {
          fact=sqrt(wsq);
          rho=b*fact;
          
          random_su2(rho,r);
          
          fact=1.0/fact;
          s[0]=fact*(r[0]*w[0]-r[1]*w[1]-r[2]*w[2]-r[3]*w[3]);
          s[1]=fact*(r[1]*w[0]+r[0]*w[1]-r[2]*w[3]+r[3]*w[2]);
          s[2]=fact*(r[2]*w[0]+r[0]*w[2]-r[3]*w[1]+r[1]*w[3]);
          s[3]=fact*(r[3]*w[0]+r[0]*w[3]-r[1]*w[2]+r[2]*w[1]);
        }
      } else random_su2(0.0,s);
      
      rotate(pu1, pu2, s);
      
      ++pu2; ++pv2;
    }
    ++pu1; ++pv1;
  }
#endif
}

