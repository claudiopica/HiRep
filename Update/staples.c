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

static suNg *u1up,*u2up,*u3up;
static suNg *u1dn,*u2dn,*u3dn;

static suNg staple, tr1, tr2;

static int ixpmu,ixpnu,ixmnu,ixpmumnu;
#ifdef PLAQ_WEIGHTS
static int IX,MU,NU;
#endif

static void add_to_v(suNg *v){
  _suNg_add_assign(*v,staple);
}


static void up_staple(void)
{
  _suNg_times_suNg(tr2,*u1up,*u2up);
  _suNg_dagger(tr1,*u3up);
  _suNg_times_suNg(staple,tr2,tr1);
#ifdef PLAQ_WEIGHTS
  if(plaq_weight!=NULL) {
    _suNg_mul(staple,plaq_weight[IX*16+NU*4+MU],staple);
  }
#endif
}


static void dn_staple(void)
{
  _suNg_times_suNg(tr2,*u2dn,*u3dn);
  _suNg_dagger(tr1,*u1dn);
  _suNg_times_suNg(staple,tr1,tr2);
#ifdef PLAQ_WEIGHTS
  if(plaq_weight!=NULL) {
    _suNg_mul(staple,plaq_weight[ixmnu*16+MU*4+NU],staple);
  }
#endif
}
   

void staples(int ix,int mu,suNg *v)
{
   int i,nu;

   _suNg_zero(*v);

   ixpmu=iup(ix,mu);
   

   for (i=1;i<4;i++)
     {
       nu=(mu+i)&0x3;
       ixpnu=iup(ix,nu);
       ixmnu=idn(ix,nu);
       ixpmumnu=idn(ixpmu,nu);
#ifdef PLAQ_WEIGHTS
       IX=ix;MU=mu;NU=nu;
#endif
       
       u1up=pu_gauge(ix,nu);
       u2up=pu_gauge(ixpnu,mu);
       u3up=pu_gauge(ixpmu,nu);   
       
       u1dn=pu_gauge(ixmnu,nu);
       u2dn=pu_gauge(ixmnu,mu);
       u3dn=pu_gauge(ixpmumnu,nu);   
       
       up_staple();
       add_to_v(v);
       dn_staple();
       add_to_v(v);
     }
}

#include "observables.h"
#include "logger.h"
void test_staples()
{
  int mu;
  suNg s, res;
  double pa=0.0;
  double ps=0.0, pl=0.;
  double tr;
  int fl;
  _DECLARE_INT_ITERATOR(ix);



  fl = 0;
  pa = avr_plaquette();
  _MASTER_FOR(&glattice,ix){
    for (mu=0; mu<4; ++mu){
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

  ps /= 4.0*(double)(6*GLB_T*GLB_X*GLB_Y*GLB_Z*NG);
  pl/=(double)(6*GLB_T*GLB_X*GLB_Y*GLB_Z*NG);
	lprintf("TESTING",50,"Staple test: ");
  if (fabs(pa-ps)<1.e-6)
		lprintf("TESTING",50,"PASSED.");
  else
		lprintf("TESTING",50,"FAILED.");
  lprintf("TESTING",50," [diff1 = %1.8e][diff2 = %1.8e][diff3 = %1.8e]\n", pa-ps, pl-ps, pa-pl);
}
