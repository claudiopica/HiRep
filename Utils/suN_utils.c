/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *   
* All rights reserved.                                                      * 
\***************************************************************************/

/*******************************************************************************
*
* File su3_utils.c
*
* Functions to project to SU(3)
*
*******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "utils.h"
#include "suN.h"
#include "representation.h"
#include "error.h"

static void normalize(suNg_vector *v)
{
   double fact;
	 _vector_prod_re_g(fact,*v,*v);
   fact=1.0f/sqrt(fact);
   _vector_mul_g(*v, fact, *v);
}


static void normalize_flt(suNg_vector_flt *v)
{
   double fact;
	 _vector_prod_re_g(fact,*v,*v);
   fact=1.0/sqrt(fact);
   _vector_mul_g(*v, fact, *v);
}


void project_to_suNg(suNg *u)
{
#ifdef WITH_QUATERNIONS
	double norm;
    
    _suNg_sqnorm(norm,*u);
    norm=sqrt(0.5*norm);
    norm=1./norm;
    _suNg_mul(*u,norm,*u);
    
#else
  int i,j;
  suNg_vector *v1,*v2;
  complex z;

  v1=(suNg_vector*)(u);
  v2=v1+1;
   
  normalize(v1);
  for (i=1; i<NG; ++i ) {
    for (j=i; j>0; --j) {
      _vector_prod_re_g(z.re,*v1, *v2);
      _vector_prod_im_g(z.im,*v1, *v2);
      _vector_project_g(*v2, z, *v1);
      ++v1;
    }
    normalize(v2);
    ++v2;
    v1=(suNg_vector*)(u);
  }
#endif
}

void project_to_suNg_flt(suNg_flt *u)
{
#ifdef WITH_QUATERNIONS
	float norm;
    
    _suNg_sqnorm(norm,*u);
    norm=sqrtf(0.5*norm);
    norm=1.f/norm;
    _suNg_mul(*u,norm,*u);
    
#else
  int i,j;
  suNg_vector_flt *v1,*v2;
  complex_flt z;

  v1=(suNg_vector_flt*)(u);
  v2=v1+1;
   
  normalize_flt(v1);
  for (i=1; i<NG; ++i ) {
    for (j=i; j>0; --j) {
      _vector_prod_re_g(z.re,*v1, *v2);
      _vector_prod_im_g(z.im,*v1, *v2);
      _vector_project_g(*v2, z, *v1);
      ++v1;
    }
    normalize_flt(v2);
    ++v2;
    v1=(suNg_vector_flt*)(u);
  }
#endif
}



void project_cooling_to_suNg(suNg* g_out, suNg* g_in, int cooling)
{
#ifdef WITH_QUATERNIONS
    error(1,1,"project_cooling_to_suNg " __FILE__,"not implemented with quaternions");
#else
  suNg Ug[3];
  suNg tmp[2];
  int k,l;
  int j,i,ncool;
  double c[NG];
  complex f[2];
  double norm;

  Ug[0]=*g_in;
  
      
  for (j=0; j<NG; j++){					
    c[j]=0.0;

    for (i=0; i<NG; i++){					
      _complex_0(f[1]);
      
      for (k=0; k<j; k++){
	_complex_0(f[0]);
	
	for (l=0; l<NG; l++){
	  _complex_mul_star_assign(f[0],(Ug[0]).c[l*NG+j], (Ug[1]).c[l*NG+k]); 
	}
	_complex_mulcr_assign(f[1],c[k],(Ug[1]).c[i*NG+k],f[0]);

      }
      
      _complex_sub(Ug[1].c[i*NG+j],Ug[0].c[i*NG+j],f[1]);
      _complex_mul_star_assign_re(c[j],Ug[1].c[i*NG+j],Ug[1].c[i*NG+j]);

    }

    c[j]= 1.0/c[j];
  }

  for(i=0;i<NG;i++)
    {
      norm=0.0;
      for(j=0;j<NG;j++){
	_complex_mul_star_assign_re(norm,Ug[1].c[i+NG*j],Ug[1].c[i+NG*j]);
      }
      
      for(j=0;j<NG;j++){
	_complex_mulr(Ug[1].c[i+NG*j],1.0/sqrt(norm),Ug[1].c[i+NG*j]);
	
      }
    }

  
  
  _suNg_dagger(Ug[2],*g_in);
 
  for (ncool=0; ncool<cooling; ncool++) 
    {
      
      _suNg_times_suNg(Ug[0],Ug[2],Ug[1]);
      
      for (i=0; i<NG; i++) {
	for (j=i+1; j<NG; j++) {
	  
	  _complex_add_star(f[0],Ug[0].c[i+NG*i],Ug[0].c[j+NG*j]);
	  _complex_sub_star(f[1],Ug[0].c[j+NG*i],Ug[0].c[i+NG*j]);
	  
	  norm = 1.0/sqrt( _complex_prod_re(f[0],f[0]) + _complex_prod_re(f[1],f[1]) );
	  
	  _complex_mulr(f[0],norm,f[0]);
	  _complex_mulr(f[1],norm,f[1]);
	  
	  _suNg_unit(tmp[0]);
	  
	  _complex_star(tmp[0].c[i+NG*i],f[0]);
	  _complex_star(tmp[0].c[i+NG*j],f[1]);
	  tmp[0].c[j+NG*j]=f[0];
	  _complex_minus(tmp[0].c[j+NG*i],f[1]);
	  
	  
	  
	  _suNg_times_suNg(tmp[1],Ug[1],tmp[0]);
	  Ug[1]=tmp[1];
	  
	  _suNg_times_suNg(tmp[1],Ug[0],tmp[0]);
	  Ug[0]=tmp[1];
	  
	}
      }
    }
  
  *g_out = Ug[1]; 
#endif
}
