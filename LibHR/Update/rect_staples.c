/****************************************************************************
* Copyright (c) 2013, Vincent Drach                                         *   
* All rights reserved.                                                      * 
\***************************************************************************/

/*******************************************************************************
* File rect_staples.c
*
* Computation of the "generalized staple" needed to compute forces with a
* rectangular term in the gauge action. 
* ____
*|    |                          
*|    |        ________      ________
*|    |       |        |    |        |       x
*|    |  +    |    ____|  + |____    |   +   |    |    +   x    ____  +   ____x    
*x            x                  x           |    |        |        |    |        |
*                                            |    |        |________|    |________|
*                                            |____|                      
*  Aup            Bup          Cup             Adn             Bdn           Cdn    
* based on staples.c from Claudio Pica
*******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "suN.h"
#include "global.h"
#include "utils.h"
#include "communications.h"
#include "update.h"

static suNg *u1up,*u2up,*u3up,*u4up,*u5up;
static suNg *u1dn,*u2dn,*u3dn,*u4dn,*u5dn;

static suNg staple, tr1, tr2,tr3;

static int ixpmu,ixpnu,ixpnupnu,ixpnupnupmu,ixpmupnu,ixpmupmu;
static int ixmnu,ixmnumnu, ixmnumnupmu,ixpmumnu,ixpmupmumnu,ixpmumnupmu;
static int ixmmu,ixmmupnu,ixmmumnu;

static void add_to_v(suNg *v){
  _suNg_add_assign(*v,staple);
}

  

void rect_staples_1x2(int ix,int mu,suNg *v)
{
	int i,nu;
	
	_suNg_zero(*v);
	
	ixpmu=iup(ix,mu);
	ixmmu=idn(ix,mu);
	ixpmupmu=iup(ixpmu,mu);
	
	for (i=1;i<4;i++)
	{
		nu=(mu+i)&0x3;
		// printf("i=%3d  mu = %3d  nu %3d\n",i,mu,nu);
		ixpnu=iup(ix,nu);
		ixpnupnu=iup(ixpnu,nu);
		ixpmupnu=iup(ixpmu,nu);
		
		ixmnu=idn(ix,nu);
		ixmnumnu=idn(ixmnu,nu);
		ixmnumnupmu=iup(ixmnumnu,mu);
		ixpmumnu=idn(ixpmu,nu);
		ixpmupmumnu=idn(ixpmupmu,nu);
		ixpmumnupmu=iup(ixpmumnu,mu);
		
		ixmmupnu=iup(ixmmu,nu);
		ixmmumnu=idn(ixmmu,nu);

#ifdef PLAQ_WEIGHTS
		error(1,1,"rect_staples_1x2","Plaq weights not implemented for 1x2 loops");
#endif
		
		// type Aup
		u1up=pu_gauge(ix,nu);
		u2up=pu_gauge(ixpnu,nu);
		u3up=pu_gauge(ixpnupnu,mu);
		u4up=pu_gauge(ixpmupnu,nu);
		u5up=pu_gauge(ixpmu,nu);
		
		// links product
		_suNg_times_suNg(tr1,*u1up,*u2up);
		_suNg_times_suNg(tr2,tr1,*u3up);
		_suNg_times_suNg(tr1,*u5up,*u4up);
		_suNg_dagger(tr3,tr1);
		_suNg_times_suNg(staple,tr2,tr3);
		
		add_to_v(v);
		
		// type Adn
		u1dn=pu_gauge(ixmnu,nu);
		u2dn=pu_gauge(ixmnumnu,nu);
		u3dn=pu_gauge(ixmnumnu,mu);
		u4dn=pu_gauge(ixmnumnupmu,nu);
		u5dn=pu_gauge(ixpmumnu,nu);
		
		// links product
		_suNg_times_suNg(tr1,*u3dn,*u4dn);
		_suNg_times_suNg(tr2,tr1,*u5dn);
		_suNg_times_suNg(tr1,*u2dn,*u1dn);
		_suNg_dagger(tr3,tr1);
		_suNg_times_suNg(staple,tr3,tr2);
		
		add_to_v(v);
		
		// Bup
		u1up=pu_gauge(ix,nu);
		u2up=pu_gauge(ixpnu,mu);
		u3up=pu_gauge(ixpmupnu,mu);
		u4up=pu_gauge(ixpmupmu,nu);
		u5up=pu_gauge(ixpmu,mu);
		
		// links product
		_suNg_times_suNg(tr1,*u1up,*u2up);
		_suNg_times_suNg(tr2,tr1,*u3up);
		_suNg_times_suNg(tr1,*u5up,*u4up);
		_suNg_dagger(tr3,tr1);
		_suNg_times_suNg(staple,tr2,tr3);
		
		add_to_v(v);
		
		
		//  Bdn
		u1dn=pu_gauge(ixmnu,nu);
		u2dn=pu_gauge(ixmnu,mu);
		u3dn=pu_gauge(ixpmumnu,mu);
		u4dn=pu_gauge(ixpmumnupmu,nu);
		u5dn=pu_gauge(ixpmu,mu);
		
		// links product
		
		_suNg_dagger(tr1,*u1dn);
		_suNg_times_suNg(tr2,tr1,*u2dn);
		_suNg_times_suNg(tr1,tr2,*u3dn);
		_suNg_times_suNg(tr2,tr1,*u4dn);
		_suNg_dagger(tr3,*u5dn);
		_suNg_times_suNg(staple,tr2,tr3);
		
		add_to_v(v);
		
		// Cup
		u1up=pu_gauge(ixmmu,mu);
		u2up=pu_gauge(ixmmu,nu);
		u3up=pu_gauge(ixmmupnu,mu);
		u4up=pu_gauge(ixpnu,mu);
		u5up=pu_gauge(ixpmu,nu);
		
		// links product
		_suNg_dagger(tr1,*u1up);
		_suNg_times_suNg(tr2,tr1,*u2up);
		_suNg_times_suNg(tr1,tr2,*u3up);
		_suNg_times_suNg(tr2,tr1,*u4up);
		_suNg_dagger(tr3,*u5up);
		_suNg_times_suNg(staple,tr2,tr3);
		
		
		add_to_v(v);
		
		// type Cdn
		u1dn=pu_gauge(ixmmu,mu);
		u2dn=pu_gauge(ixmmumnu,nu);
		u3dn=pu_gauge(ixmmumnu,mu);
		u4dn=pu_gauge(ixmnu,mu);
		u5dn=pu_gauge(ixpmumnu,nu);
		
		// links product
		_suNg_times_suNg(tr1,*u3dn,*u4dn);
		_suNg_times_suNg(tr2,tr1,*u5dn);
		_suNg_times_suNg(tr1,*u2dn,*u1dn);
		_suNg_dagger(tr3,tr1);
		_suNg_times_suNg(staple,tr3,tr2);
		
		
		add_to_v(v);
		
	}
	
	
}

#include "observables.h"
#include "logger.h"
void test_rect_staples_1x2()
{
  int mu;
  suNg s, res;
  double pa=0.0;
  double ps=0.0, pl=0.;
  double tr;
  // _DECLARE_INT_ITERATOR(ix);
	
  pa = avr_rect_1x2();
  _MASTER_FOR(&glattice,ix){
    for (mu=0; mu<4; ++mu){
      rect_staples_1x2(ix, mu, &res);
      //													 _suNg_dagger(s, res);
      s=res;
      _suNg_times_suNg_dagger(res, *(pu_gauge(ix, mu)), s);
      _suNg_trace_re(tr,res);
      ps += tr;
    }
    pl+=local_rect_1x2(ix);
  }
  global_sum(&ps,1);
  global_sum(&pl,1);

  ps /= 3*24*NG*GLB_VOLUME;
  pl /= 2*6*NG*GLB_VOLUME;
  
  lprintf("TESTING",0,"rect_1x2 staple test: ");
  if (fabs(pa-ps)<1.e-6)
    lprintf("TESTING",0,"PASSED.");
  else
    lprintf("TESTING",0,"FAILED.");
  lprintf("TESTING",0," [diff1 = %1.8e][diff2 = %1.8e][diff3 = %1.8e]\n", pa-ps, pl-ps, pa-pl);
  lprintf("TESTING",0," [avr_rect_1x2 pa = %1.8e][rect_staples ps =  %1.8e][local rect pl = %1.8e]\n", pa, ps, pl);
}

