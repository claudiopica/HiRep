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


static double s[4],w[4];
static suNg_vector *pu1,*pu2,*pv1,*pv2;

static void rotate(void)
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

static void wmatrix(void) 
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
	  fact=(float)(sqrt((double)(wsq))); // k
	  rho=b*fact;						// total weight for a_0
	  
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


void cabmar_so(double beta,suNg *u,suNg *v,int type) 
// u is the link 
// v is the sum of staples
// if type = 0 -> heatbath
// if type = 1 -> Overrelaxation
{
  int i,j,k;
  double b,bsq,w,wsq,theta;
  suNg uv;	// Will contain link times staple
  double COSuv,SINuv,COSran,SINran,COSup,SINup;
  double Urow; 
  
  
  const double invng = 1. / (double) NG;
  
  b=invng*beta; 
  bsq=b*b;
 
  
  
      _suNg_times_suNg_dagger(uv,*u,*v);
  
  for (i=0; i<NG-1; ++i) { // Looping over subgroups
    for (j=i+1; j<NG; ++j) {
    
    // Reading components
    COSuv	= (uv.c[i*NG+i]+uv.c[j*NG+j])/2.;  // Assuming rows*NG+cols
    SINuv	= (uv.c[i*NG+j]-uv.c[j*NG+i])/2.;
    wsq		= COSuv*COSuv+SINuv*SINuv;
    w		= sqrt(wsq);
    COSuv	= COSuv/w;
    SINuv	= SINuv/w;
      
      
	if ((bsq*wsq)>1.0e-16f) {
		if (type==0) {
      	  
      		theta = random_so2(2*b*w);
      		COSran = cos(theta);
      		SINran = sin(theta);

      		// Rotating back the generated rotation matrix
      		COSup=COSran*COSuv+SINran*SINuv;
      		SINup=COSuv*SINran-COSran*SINuv;

      		// Rotate the subgroup using these angles
      		for (k=0;k<NG;k++){
      			Urow=u->c[NG*i+k];
      			u->c[NG*i+k]=COSup*Urow+SINup*u->c[NG*j+k];
      			u->c[NG*j+k]=COSup*u->c[NG*j+k]-SINup*Urow;
      			
      			// Update also the plaquettes for the next read
      			Urow=uv.c[NG*i+k];
      			uv.c[NG*i+k]=COSup*Urow+SINup*uv.c[NG*j+k];
      			uv.c[NG*j+k]=COSup*uv.c[NG*j+k]-SINup*Urow;
      		}
      		
 
      		
      
		} else {
			// Microcanonical Overrelaxation 

      		// Rotating back the generated rotation matrix
      		COSup=COSuv;
      		SINup=-SINuv;

      		for (int l=0;l<2;l++){
      			for (k=0;k<NG;k++){
      				Urow=u->c[NG*i+k];
      				u->c[NG*i+k]=COSup*Urow+SINup*u->c[NG*j+k];
      				u->c[NG*j+k]=COSup*u->c[NG*j+k]-SINup*Urow;
      			
      				// Update also the plaquettes for the next read
      				Urow=uv.c[NG*i+k];
      				uv.c[NG*i+k]=COSup*Urow+SINup*uv.c[NG*j+k];
      				uv.c[NG*j+k]=COSup*uv.c[NG*j+k]-SINup*Urow;
      			}
      		}
      		
			   
		}      
      } else {
			// If the element is approximately orthogonal to the subgroup
			}
    }// 
   }// Stop looping over subgroups
 
 
}

