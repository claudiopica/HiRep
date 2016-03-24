/***************************************************************************\
* Copyright (c) 2013, Vincent Drach                                         *
* All rights reserved.                                                      *
\***************************************************************************/


/*******************************************************************************
 *
 * File avr_rect.c
 *
 * Routines for the average 1x2 Wilson Loops
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

double rect_1x2(int ix,int mu,int nu)
{
				int ix_pmu,ix_pmu_pnu,ix_pnu,ix_pnu_pnu;
				double p;
				suNg *v1,*v2,*v3,*v4,*v5,*v6,w1,w2,w3,w4;

				ix_pmu=iup(ix,mu);
				ix_pmu_pnu=iup(ix_pmu,nu);
				ix_pnu=iup(ix,nu);
				ix_pnu_pnu=iup(ix_pnu,nu);

				v1=pu_gauge(ix,mu);
				v2=pu_gauge(ix_pmu,nu);
				v3=pu_gauge(ix_pmu_pnu,nu);
				v4=pu_gauge(ix_pnu_pnu,mu);
				v5=pu_gauge(ix_pnu,nu);
				v6=pu_gauge(ix,nu);

				_suNg_times_suNg(w1,(*v1),(*v2));
				_suNg_times_suNg(w2,w1,(*v3));/*v1 v2 v3*/
				_suNg_times_suNg(w3,(*v6 ),(*v5));
				_suNg_times_suNg(w4,w3,(*v4));/* v6 v5 v4 */ 
				_suNg_times_suNg_dagger(w3,w2,w4);/* v1 v2 v3 v4d v5d v6d */  

				_suNg_trace_re(p,w3);

#ifdef PLAQ_WEIGHTS
				error("Plaq weights in rect Not implemented whatever it is ... ");  
#else
				return p;
#endif
}
void crect_1x2(complex *ret,int ix,int mu,int nu)
{				
				int ix_pmu,ix_pmu_pnu,ix_pnu,ix_pnu_pnu;
				double p;
				suNg *v1,*v2,*v3,*v4,*v5,*v6,w1,w2,w3,w4;

				ix_pmu=iup(ix,mu);
				ix_pmu_pnu=iup(ix_pmu,nu);
				ix_pnu=iup(ix,nu);
				ix_pnu_pnu=iup(ix_pnu,nu);

				v1=pu_gauge(ix,mu);
				v2=pu_gauge(ix_pmu,nu);
				v3=pu_gauge(ix_pmu_pnu,nu);
				v4=pu_gauge(ix_pnu_pnu,mu);
				v5=pu_gauge(ix_pnu,nu);
				v6=pu_gauge(ix,nu);

				_suNg_times_suNg(w1,(*v1),(*v2));
				_suNg_times_suNg(w2,w1,(*v3));/*v1 v2 v3*/
				_suNg_times_suNg(w3,(*v6 ),(*v5));
				_suNg_times_suNg(w4,w3,(*v4));/* v6 v5 v4 */ 

				_suNg_times_suNg_dagger(w3,w2,w4);/* v1 v2 v3 v4d v5d v6d */  

				_suNg_trace_re(ret->re,w3);
#ifdef GAUGE_SON
				ret->im=0;
#else
				_suNg_trace_im(ret->im,w3);
#endif

#ifdef PLAQ_WEIGHTS
				error("NOT TESTED");
#endif

}


double avr_rect_1x2()
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
      pa+=rect_1x2(ix,1,0);
      pa+=rect_1x2(ix,2,0);
      pa+=rect_1x2(ix,2,1);
      pa+=rect_1x2(ix,3,0);
      pa+=rect_1x2(ix,3,1);
      pa+=rect_1x2(ix,3,2);
      pa+=rect_1x2(ix,0,1);
      pa+=rect_1x2(ix,0,2);
      pa+=rect_1x2(ix,1,2);
      pa+=rect_1x2(ix,0,3);
      pa+=rect_1x2(ix,1,3);
      pa+=rect_1x2(ix,2,3);
    
    }
  }

  global_sum(&pa, 1);

  return pa/(12.*NG)/GLB_VOLUME;

}

void full_rect_1x2()
{
  complex pa[12];
  for(int k=0;k<12;k++)
    {
      pa[k].re=0.;
      pa[k].im=0.;
    }
  
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
	  crect_1x2(&tmp,ix,1,0); _complex_add_assign(pa[0],tmp);
	  crect_1x2(&tmp,ix,2,0); _complex_add_assign(pa[1],tmp);
	  crect_1x2(&tmp,ix,2,1); _complex_add_assign(pa[2],tmp);
	  crect_1x2(&tmp,ix,3,0); _complex_add_assign(pa[3],tmp);
	  crect_1x2(&tmp,ix,3,1); _complex_add_assign(pa[4],tmp);
	  crect_1x2(&tmp,ix,3,2); _complex_add_assign(pa[5],tmp);
	  crect_1x2(&tmp,ix,0,1); _complex_add_assign(pa[6],tmp);
	  crect_1x2(&tmp,ix,0,2); _complex_add_assign(pa[7],tmp);
	  crect_1x2(&tmp,ix,1,2); _complex_add_assign(pa[8],tmp);
	  crect_1x2(&tmp,ix,0,3); _complex_add_assign(pa[9],tmp);
	  crect_1x2(&tmp,ix,1,3); _complex_add_assign(pa[10],tmp);
	  crect_1x2(&tmp,ix,2,3); _complex_add_assign(pa[11],tmp);

		
	}
    }

	

  global_sum((double*)pa,12);
  for(int k = 0; k < 12; k++)
    {
      pa[k].re /= GLB_VOLUME*NG;
      pa[k].im /= GLB_VOLUME*NG;
    }
  
  lprintf("PLAQ",0,"Rect 1x2 ( %d , %d) = ( %f , %f )\n",1,0,pa[0].re,pa[0].im);
  lprintf("PLAQ",0,"Rect 1x2 ( %d , %d) = ( %f , %f )\n",2,0,pa[1].re,pa[1].im);
  lprintf("PLAQ",0,"Rect 1x2( %d , %d) = ( %f , %f )\n",2,1,pa[2].re,pa[2].im);
  lprintf("PLAQ",0,"Rect 1x2( %d , %d) = ( %f , %f )\n",3,0,pa[3].re,pa[3].im);
  lprintf("PLAQ",0,"Rect 1x2( %d , %d) = ( %f , %f )\n",3,1,pa[4].re,pa[4].im);
  lprintf("PLAQ",0,"Rect 1x2( %d , %d) = ( %f , %f )\n",3,2,pa[5].re,pa[5].im);
  
  lprintf("PLAQ",0,"Rect 1x2 ( %d , %d) = ( %f , %f )\n",0,1,pa[6].re,pa[6].im);
  lprintf("PLAQ",0,"Rect 1x2 ( %d , %d) = ( %f , %f )\n",0,2,pa[7].re,pa[7].im);
  lprintf("PLAQ",0,"Rect 1x2( %d , %d) = ( %f , %f )\n",1,2,pa[8].re,pa[8].im);
  lprintf("PLAQ",0,"Rect 1x2( %d , %d) = ( %f , %f )\n",0,3,pa[9].re,pa[9].im);
  lprintf("PLAQ",0,"Rect 1x2( %d , %d) = ( %f , %f )\n",1,3,pa[10].re,pa[10].im);
  lprintf("PLAQ",0,"Rect 1x2( %d , %d) = ( %f , %f )\n",2,3,pa[11].re,pa[11].im);
}


double local_rect_1x2(int ix)
{
  double pa;
  
  pa=rect_1x2(ix,1,0);
  pa+=rect_1x2(ix,2,0);
  pa+=rect_1x2(ix,2,1);
  pa+=rect_1x2(ix,3,0);
  pa+=rect_1x2(ix,3,1);
  pa+=rect_1x2(ix,3,2);
  pa+=rect_1x2(ix,0,1);
  pa+=rect_1x2(ix,0,2);
  pa+=rect_1x2(ix,1,2);
  pa+=rect_1x2(ix,0,3);
  pa+=rect_1x2(ix,1,3);
  pa+=rect_1x2(ix,2,3);
				
  
  return pa;

}

