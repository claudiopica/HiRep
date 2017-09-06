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
#include "utils.h"
#include "update.h"
#include "memory.h"
#include "random.h"
#include "dirac.h"
#include "representation.h"
#include "linear_algebra.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "logger.h"
#include "communications.h"


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

#ifdef PLAQ_WEIGHTS
  if(plaq_weight==NULL) return p;
  return plaq_weight[ix*16+mu*4+nu]*p;
#else
  return p;
#endif
}

void cplaq(complex *ret,int ix,int mu,int nu)
{
  int iy,iz;
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
#ifdef GAUGE_SON
  ret->im=0;
#else
  _suNg_trace_im(ret->im,w3);
#endif

#ifdef PLAQ_WEIGHTS
  if(plaq_weight!=NULL) {
    ret->re *= plaq_weight[ix*16+mu*4+nu];
    ret->im *= plaq_weight[ix*16+mu*4+nu];
  }
#endif

}


double avr_plaquette()
{
  double pa=0.;

  _PIECE_FOR(&glattice,ixp) {
    if(ixp==glattice.inner_master_pieces) {
      _OMP_PRAGMA( master )
      /* wait for gauge field to be transfered */
      complete_gf_sendrecv(u_gauge);
      _OMP_PRAGMA( barrier )
    }
    _SITE_FOR_SUM(&glattice,ixp,ix,pa) {
      pa+=plaq(ix,1,0);
      pa+=plaq(ix,2,0);
      pa+=plaq(ix,2,1);
      pa+=plaq(ix,3,0);
      pa+=plaq(ix,3,1);
      pa+=plaq(ix,3,2);
    }
  }

  global_sum(&pa, 1);

#ifdef BC_T_OPEN
	pa /= 6.0*NG*GLB_VOLUME*(GLB_T-1)/GLB_T;
#else
	pa /= 6.0*NG*GLB_VOLUME;
#endif

  return pa;

}

void full_plaquette()
{
	complex pa[6];
	double r0re = 0, r0im = 0;
	double r1re = 0, r1im = 0;
	double r2re = 0, r2im = 0;
	double r3re = 0, r3im = 0;
	double r4re = 0, r4im = 0;
	double r5re = 0, r5im = 0;

	_PIECE_FOR(&glattice,ixp)
	{
		if(ixp == glattice.inner_master_pieces)
		{
			_OMP_PRAGMA( master )
			/* wait for gauge field to be transfered */
			complete_gf_sendrecv(u_gauge);
			_OMP_PRAGMA( barrier )
		}

		_SITE_FOR_SUM(&glattice,ixp,ix,r0re,r0im,r1re,r1im,r2re,r2im,r3re,r3im,r4re,r4im,r5re,r5im)
		{
			complex tmp;

			cplaq(&tmp,ix,1,0);
			r0re += tmp.re;
			r0im += tmp.im;

			cplaq(&tmp,ix,2,0);
			r1re += tmp.re;
			r1im += tmp.im;

			cplaq(&tmp,ix,2,1);
			r2re += tmp.re;
			r2im += tmp.im;

			cplaq(&tmp,ix,3,0);
			r3re += tmp.re;
			r3im += tmp.im;

			cplaq(&tmp,ix,3,1);
			r4re += tmp.re;
			r4im += tmp.im;

			cplaq(&tmp,ix,3,2);
			r5re += tmp.re;
			r5im += tmp.im;

		  /*if(twbc_plaq[ix*16+2*4+1]==-1 &&
		    twbc_plaq[ix*16+3*4+1]==-1 &&
		    twbc_plaq[ix*16+3*4+2]==-1) {
		  cplaq(&tmp,ix,1,0);
		  lprintf("LOCPL",0,"Plaq( %d , %d , %d ) = ( %f , %f )\n",t,1,0,tmp.re,tmp.im);
		  cplaq(&tmp,ix,2,0);
		  lprintf("LOCPL",0,"Plaq( %d , %d , %d ) = ( %f , %f )\n",t,2,0,tmp.re,tmp.im);
		  cplaq(&tmp,ix,2,1);
		  lprintf("LOCPL",0,"Plaq( %d , %d , %d ) = ( %f , %f )\n",t,2,1,tmp.re,tmp.im);
		  cplaq(&tmp,ix,3,0);
		  lprintf("LOCPL",0,"Plaq( %d , %d , %d ) = ( %f , %f )\n",t,3,0,tmp.re,tmp.im);
		  cplaq(&tmp,ix,3,1);
		  lprintf("LOCPL",0,"Plaq( %d , %d , %d ) = ( %f , %f )\n",t,3,1,tmp.re,tmp.im);
		  cplaq(&tmp,ix,3,2);
		  lprintf("LOCPL",0,"Plaq( %d , %d , %d ) = ( %f , %f )\n",t,3,2,tmp.re,tmp.im);
		  t++;
		    } */
		}
	}

	pa[0].re=r0re; pa[0].im=r0im;
	pa[1].re=r1re; pa[1].im=r1im;
	pa[2].re=r2re; pa[2].im=r2im;
	pa[3].re=r3re; pa[3].im=r3im;
	pa[4].re=r4re; pa[4].im=r4im;
	pa[5].re=r5re; pa[5].im=r5im;

	global_sum((double*)pa,12);
	for(int k = 0; k < 6; k++)
	{
#ifdef BC_T_OPEN
		pa[k].re /= NG*GLB_VOLUME*(GLB_T-1)/GLB_T;
		pa[k].im /= NG*GLB_VOLUME*(GLB_T-1)/GLB_T;
#else
		pa[k].re /= NG*GLB_VOLUME;
		pa[k].im /= NG*GLB_VOLUME;
#endif
	}

	lprintf("PLAQ",0,"Plaq(%d,%d) = ( %f , %f )\n",1,0,pa[0].re,pa[0].im);
	lprintf("PLAQ",0,"Plaq(%d,%d) = ( %f , %f )\n",2,0,pa[1].re,pa[1].im);
	lprintf("PLAQ",0,"Plaq(%d,%d) = ( %f , %f )\n",2,1,pa[2].re,pa[2].im);
	lprintf("PLAQ",0,"Plaq(%d,%d) = ( %f , %f )\n",3,0,pa[3].re,pa[3].im);
	lprintf("PLAQ",0,"Plaq(%d,%d) = ( %f , %f )\n",3,1,pa[4].re,pa[4].im);
	lprintf("PLAQ",0,"Plaq(%d,%d) = ( %f , %f )\n",3,2,pa[5].re,pa[5].im);
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

void full_momenta(suNg_av_field *momenta){
  scalar_field *la=alloc_sfield(1, &glattice); 
  _MASTER_FOR(&glattice,i) {
    double a=0., tmp;
    /* Momenta */
    for (int j=0;j<4;++j) {
      suNg_algebra_vector *cmom=momenta->ptr+coord_to_index(i,j);
      _algebra_vector_sqnorm_g(tmp,*cmom); 
      a+=tmp; /* this must be positive */
    }
    a*=0.5*_FUND_NORM2;
    *_FIELD_AT(la,i)=a;
  }
  double mom=0;
  _MASTER_FOR_SUM(la->type,i,mom) {
    mom += *_FIELD_AT(la,i);
  }  
  lprintf("MOMENTA",0,"%1.8g\n",mom);
  free_sfield(la);
}
