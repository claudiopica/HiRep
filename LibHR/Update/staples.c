/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *   
* All rights reserved.                                                      * 
\***************************************************************************/

/*******************************************************************************
*
* File staples.c
*
* Computation of the "staples"
*
*******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "suN.h"
#include "global.h"
#include "communications.h"
#include "update.h"

void staples(int ix,int mu,suNg *v)
{
  suNg staple, tr1, tr2;
  
  int ixpmu=iup(ix,mu);

  _suNg_zero(*v);
  
  for (int i=1;i<4;i++) {
    int nu=(mu+i)&0x3;
    int ixpnu=iup(ix,nu);
    int ixmnu=idn(ix,nu);
    int ixpmumnu=idn(ixpmu,nu);

    //Up Staple
    _suNg_times_suNg(tr2,*pu_gauge(ix,nu),*pu_gauge(ixpnu,mu));
    _suNg_dagger(tr1,*pu_gauge(ixpmu,nu));
    _suNg_times_suNg(staple,tr2,tr1);
#ifdef PLAQ_WEIGHTS
    if(plaq_weight!=NULL) {
      _suNg_mul(staple,plaq_weight[ix*16+nu*4+mu],staple);
    }
#endif
    _suNg_add_assign(*v,staple);
    
    //Down Staple
    _suNg_times_suNg(tr2,*pu_gauge(ixmnu,mu),*pu_gauge(ixpmumnu,nu));
    _suNg_dagger(tr1,*pu_gauge(ixmnu,nu));
    _suNg_times_suNg(staple,tr1,tr2);
#ifdef PLAQ_WEIGHTS
    if(plaq_weight!=NULL) {
      _suNg_mul(staple,plaq_weight[ixmnu*16+mu*4+nu],staple);
    }
#endif
    _suNg_add_assign(*v,staple);
  }
}

#if 0
#include "observables.h"
#include "logger.h"
void test_staples()
{
  double pa=0.0;
  double ps=0.0, pl=0.;
  _DECLARE_INT_ITERATOR(ix);

  pa = avr_plaquette();
  _MASTER_FOR_SUM(&glattice,ix,ps,pl){
    suNg s, res;
    for (int mu=0; mu<4; ++mu){
      double tr;
      staples(ix, mu, &res);
      /* _suNg_dagger(s, res); */
      s=res;
      _suNg_times_suNg_dagger(res, *(pu_gauge(ix, mu)), s);
      _suNg_trace_re(tr,res);
      ps += tr;
    }
    pl+=local_plaq(ix);
  }
  global_sum(&ps,1);
  global_sum(&pl,1);

  ps /= 24*NG*GLB_VOLUME;
  pl /= 6*NG*GLB_VOLUME;
	lprintf("TESTING",50,"Staple test: ");
  if (fabs(pa-ps)<1.e-6)
		lprintf("TESTING",50,"PASSED.");
  else
		lprintf("TESTING",50,"FAILED.");
  lprintf("TESTING",50," [diff1 = %1.8e][diff2 = %1.8e][diff3 = %1.8e]\n", pa-ps, pl-ps, pa-pl);
}
#endif