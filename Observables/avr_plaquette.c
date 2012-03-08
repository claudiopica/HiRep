/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *   
* All rights reserved.                                                      * 
\***************************************************************************/

/*******************************************************************************
*
* File plaquette.c
*
* Routines for the average plaquette
*
*******************************************************************************/

#include "global.h"
#include "suN.h"
#include "communications.h"
#include "logger.h"

double plaq(int ix,int mu,int nu)
{
  int iy,iz;
  double p;
  suNg *v1,*v2,*v3,*v4,w1,w2,w3;

  iy=iup(ix,mu);
  iz=iup(ix,nu);

  v1=pu_gauge(ix,mu);
  v2=pu_gauge(iy,nu);
  v3=pu_gauge(iz,mu);
  v4=pu_gauge(ix,nu);

  _suNg_times_suNg(w1,(*v1),(*v2));
  _suNg_times_suNg(w2,(*v4),(*v3));
  _suNg_times_suNg_dagger(w3,w1,w2);      

  _suNg_trace_re(p,w3);
  return p;
}


void cplaq(complex *ret,int ix,int mu,int nu)
{
  int iy,iz;;
  suNg *v1,*v2,*v3,*v4,w1,w2,w3;

  iy=iup(ix,mu);
  iz=iup(ix,nu);

  v1=pu_gauge(ix,mu);
  v2=pu_gauge(iy,nu);
  v3=pu_gauge(iz,mu);
  v4=pu_gauge(ix,nu);

  _suNg_times_suNg(w1,(*v1),(*v2));
  _suNg_times_suNg(w2,(*v4),(*v3));
  _suNg_times_suNg_dagger(w3,w1,w2);      
      
   _suNg_trace_re(ret->re,w3);
   _suNg_trace_im(ret->im,w3);
}


double avr_plaquette()
{
  _DECLARE_INT_ITERATOR(ix);
  double pa=0.;

  _PIECE_FOR(&glattice,ix) {
    _SITE_FOR(&glattice,ix) {
      pa+=plaq(ix,1,0);
      pa+=plaq(ix,2,0);
      pa+=plaq(ix,2,1);
      pa+=plaq(ix,3,0);
      pa+=plaq(ix,3,1);
      pa+=plaq(ix,3,2);
    }
    if(_PIECE_INDEX(ix)==0) {
      /* wait for gauge field to be transfered */
      complete_gf_sendrecv(u_gauge);
    }
  }

  global_sum(&pa, 1);

  return pa/(double)(6*GLB_T*GLB_X*GLB_Y*GLB_Z*NG);

}

void full_plaquette()
{
  _DECLARE_INT_ITERATOR(ix);
  int k;
  complex pa[6];

  for(k=0;k<6;k++)
    pa[k].re=pa[k].im=0.;

  _PIECE_FOR(&glattice,ix) {
    _SITE_FOR(&glattice,ix) {
      complex tmp;
		  cplaq(&tmp,ix,1,0); _complex_add_assign(pa[0],tmp);
		  cplaq(&tmp,ix,2,0); _complex_add_assign(pa[1],tmp);
		  cplaq(&tmp,ix,2,1); _complex_add_assign(pa[2],tmp);
		  cplaq(&tmp,ix,3,0); _complex_add_assign(pa[3],tmp);
		  cplaq(&tmp,ix,3,1); _complex_add_assign(pa[4],tmp);
		  cplaq(&tmp,ix,3,2); _complex_add_assign(pa[5],tmp);
    }
    if(_PIECE_INDEX(ix)==0) {
      /* wait for gauge field to be transfered */
      complete_gf_sendrecv(u_gauge);
    }
  }

  global_sum((double*)pa,12);
  for(k=0;k<6;k++) {
    pa[k].re /= GLB_T*GLB_X*GLB_Y*GLB_Z*NG;
    pa[k].im /= GLB_T*GLB_X*GLB_Y*GLB_Z*NG;
  }

  lprintf("PLAQ",0,"Plaq( %d , %d) = ( %.10f , %.10f )\n",1,0,pa[0].re,pa[0].im);
  lprintf("PLAQ",0,"Plaq( %d , %d) = ( %.10f , %.10f )\n",2,0,pa[1].re,pa[1].im);
  lprintf("PLAQ",0,"Plaq( %d , %d) = ( %.10f , %.10f )\n",2,1,pa[2].re,pa[2].im);
  lprintf("PLAQ",0,"Plaq( %d , %d) = ( %.10f , %.10f )\n",3,0,pa[3].re,pa[3].im);
  lprintf("PLAQ",0,"Plaq( %d , %d) = ( %.10f , %.10f )\n",3,1,pa[4].re,pa[4].im);
  lprintf("PLAQ",0,"Plaq( %d , %d) = ( %.10f , %.10f )\n",3,2,pa[5].re,pa[5].im);

}

double local_plaq(int ix)
{
  double pa;

  pa=plaq(ix,1,0);
  pa+=plaq(ix,2,0);
  pa+=plaq(ix,2,1);
  pa+=plaq(ix,3,0);
  pa+=plaq(ix,3,1);
  pa+=plaq(ix,3,2);

  return pa;

}
