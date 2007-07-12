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


static suNg *u1up,*u2up,*u3up;
static suNg *u1dn,*u2dn,*u3dn;


/*
static suNg_vector psi1,psi2,chi1,chi2;

static void add_to_v(suNg *v)
{
   (*v).c1_1.re+=psi1.c1.re;
   (*v).c1_1.im+=psi1.c1.im;
   (*v).c2_1.re+=psi1.c2.re;
   (*v).c2_1.im+=psi1.c2.im;
   (*v).c3_1.re+=psi1.c3.re;
   (*v).c3_1.im+=psi1.c3.im;

   (*v).c1_2.re+=psi2.c1.re;
   (*v).c1_2.im+=psi2.c1.im;
   (*v).c2_2.re+=psi2.c2.re;
   (*v).c2_2.im+=psi2.c2.im;
   (*v).c3_2.re+=psi2.c3.re;
   (*v).c3_2.im+=psi2.c3.im;

   (*v).c1_3.re+=
      (psi1.c2.re*psi2.c3.re-psi1.c2.im*psi2.c3.im)-
      (psi2.c2.re*psi1.c3.re-psi2.c2.im*psi1.c3.im);
   (*v).c1_3.im+=
      (psi2.c2.re*psi1.c3.im+psi2.c2.im*psi1.c3.re)-
      (psi1.c2.re*psi2.c3.im+psi1.c2.im*psi2.c3.re);
   (*v).c2_3.re+=
      (psi1.c3.re*psi2.c1.re-psi1.c3.im*psi2.c1.im)-
      (psi2.c3.re*psi1.c1.re-psi2.c3.im*psi1.c1.im);
   (*v).c2_3.im+=
      (psi2.c3.re*psi1.c1.im+psi2.c3.im*psi1.c1.re)-
      (psi1.c3.re*psi2.c1.im+psi1.c3.im*psi2.c1.re);
   (*v).c3_3.re+=
      (psi1.c1.re*psi2.c2.re-psi1.c1.im*psi2.c2.im)-
      (psi2.c1.re*psi1.c2.re-psi2.c1.im*psi1.c2.im);
   (*v).c3_3.im+=
      (psi2.c1.re*psi1.c2.im+psi2.c1.im*psi1.c2.re)-
      (psi1.c1.re*psi2.c2.im+psi1.c1.im*psi2.c2.re);
}


static void up_staple(void)
{
   chi1.c1.re=
      (*u2up).c1_1.re*(*u3up).c1_1.re+(*u2up).c1_1.im*(*u3up).c1_1.im+
      (*u2up).c1_2.re*(*u3up).c1_2.re+(*u2up).c1_2.im*(*u3up).c1_2.im+ 
      (*u2up).c1_3.re*(*u3up).c1_3.re+(*u2up).c1_3.im*(*u3up).c1_3.im;
   chi1.c1.im=
      (*u2up).c1_1.im*(*u3up).c1_1.re-(*u2up).c1_1.re*(*u3up).c1_1.im+
      (*u2up).c1_2.im*(*u3up).c1_2.re-(*u2up).c1_2.re*(*u3up).c1_2.im+ 
      (*u2up).c1_3.im*(*u3up).c1_3.re-(*u2up).c1_3.re*(*u3up).c1_3.im;
   chi1.c2.re=
      (*u2up).c2_1.re*(*u3up).c1_1.re+(*u2up).c2_1.im*(*u3up).c1_1.im+ 
      (*u2up).c2_2.re*(*u3up).c1_2.re+(*u2up).c2_2.im*(*u3up).c1_2.im+
      (*u2up).c2_3.re*(*u3up).c1_3.re+(*u2up).c2_3.im*(*u3up).c1_3.im;
   chi1.c2.im=
      (*u2up).c2_1.im*(*u3up).c1_1.re-(*u2up).c2_1.re*(*u3up).c1_1.im+
      (*u2up).c2_2.im*(*u3up).c1_2.re-(*u2up).c2_2.re*(*u3up).c1_2.im+
      (*u2up).c2_3.im*(*u3up).c1_3.re-(*u2up).c2_3.re*(*u3up).c1_3.im;
   chi1.c3.re=
      (*u2up).c3_1.re*(*u3up).c1_1.re+(*u2up).c3_1.im*(*u3up).c1_1.im+
      (*u2up).c3_2.re*(*u3up).c1_2.re+(*u2up).c3_2.im*(*u3up).c1_2.im+ 
      (*u2up).c3_3.re*(*u3up).c1_3.re+(*u2up).c3_3.im*(*u3up).c1_3.im;
   chi1.c3.im=
      (*u2up).c3_1.im*(*u3up).c1_1.re-(*u2up).c3_1.re*(*u3up).c1_1.im+ 
      (*u2up).c3_2.im*(*u3up).c1_2.re-(*u2up).c3_2.re*(*u3up).c1_2.im+
      (*u2up).c3_3.im*(*u3up).c1_3.re-(*u2up).c3_3.re*(*u3up).c1_3.im;

   chi2.c1.re=
      (*u2up).c1_1.re*(*u3up).c2_1.re+(*u2up).c1_1.im*(*u3up).c2_1.im+
      (*u2up).c1_2.re*(*u3up).c2_2.re+(*u2up).c1_2.im*(*u3up).c2_2.im+ 
      (*u2up).c1_3.re*(*u3up).c2_3.re+(*u2up).c1_3.im*(*u3up).c2_3.im;
   chi2.c1.im=
      (*u2up).c1_1.im*(*u3up).c2_1.re-(*u2up).c1_1.re*(*u3up).c2_1.im+
      (*u2up).c1_2.im*(*u3up).c2_2.re-(*u2up).c1_2.re*(*u3up).c2_2.im+ 
      (*u2up).c1_3.im*(*u3up).c2_3.re-(*u2up).c1_3.re*(*u3up).c2_3.im;
   chi2.c2.re=
      (*u2up).c2_1.re*(*u3up).c2_1.re+(*u2up).c2_1.im*(*u3up).c2_1.im+ 
      (*u2up).c2_2.re*(*u3up).c2_2.re+(*u2up).c2_2.im*(*u3up).c2_2.im+
      (*u2up).c2_3.re*(*u3up).c2_3.re+(*u2up).c2_3.im*(*u3up).c2_3.im;
   chi2.c2.im=
      (*u2up).c2_1.im*(*u3up).c2_1.re-(*u2up).c2_1.re*(*u3up).c2_1.im+
      (*u2up).c2_2.im*(*u3up).c2_2.re-(*u2up).c2_2.re*(*u3up).c2_2.im+
      (*u2up).c2_3.im*(*u3up).c2_3.re-(*u2up).c2_3.re*(*u3up).c2_3.im;
   chi2.c3.re=
      (*u2up).c3_1.re*(*u3up).c2_1.re+(*u2up).c3_1.im*(*u3up).c2_1.im+
      (*u2up).c3_2.re*(*u3up).c2_2.re+(*u2up).c3_2.im*(*u3up).c2_2.im+ 
      (*u2up).c3_3.re*(*u3up).c2_3.re+(*u2up).c3_3.im*(*u3up).c2_3.im;
   chi2.c3.im=
      (*u2up).c3_1.im*(*u3up).c2_1.re-(*u2up).c3_1.re*(*u3up).c2_1.im+ 
      (*u2up).c3_2.im*(*u3up).c2_2.re-(*u2up).c3_2.re*(*u3up).c2_2.im+
      (*u2up).c3_3.im*(*u3up).c2_3.re-(*u2up).c3_3.re*(*u3up).c2_3.im;

   _suNg_multiply(psi1,(*u1up),chi1);
   _suNg_multiply(psi2,(*u1up),chi2);   
}


static void dn_staple(void)
{
   chi1.c1.re=
      (*u2dn).c1_1.re*(*u3dn).c1_1.re-(*u2dn).c1_1.im*(*u3dn).c1_1.im+
      (*u2dn).c1_2.re*(*u3dn).c2_1.re-(*u2dn).c1_2.im*(*u3dn).c2_1.im+ 
      (*u2dn).c1_3.re*(*u3dn).c3_1.re-(*u2dn).c1_3.im*(*u3dn).c3_1.im;
   chi1.c1.im=
      (*u2dn).c1_1.im*(*u3dn).c1_1.re+(*u2dn).c1_1.re*(*u3dn).c1_1.im+
      (*u2dn).c1_2.im*(*u3dn).c2_1.re+(*u2dn).c1_2.re*(*u3dn).c2_1.im+ 
      (*u2dn).c1_3.im*(*u3dn).c3_1.re+(*u2dn).c1_3.re*(*u3dn).c3_1.im;
   chi1.c2.re=
      (*u2dn).c2_1.re*(*u3dn).c1_1.re-(*u2dn).c2_1.im*(*u3dn).c1_1.im+ 
      (*u2dn).c2_2.re*(*u3dn).c2_1.re-(*u2dn).c2_2.im*(*u3dn).c2_1.im+
      (*u2dn).c2_3.re*(*u3dn).c3_1.re-(*u2dn).c2_3.im*(*u3dn).c3_1.im;
   chi1.c2.im=
      (*u2dn).c2_1.im*(*u3dn).c1_1.re+(*u2dn).c2_1.re*(*u3dn).c1_1.im+
      (*u2dn).c2_2.im*(*u3dn).c2_1.re+(*u2dn).c2_2.re*(*u3dn).c2_1.im+
      (*u2dn).c2_3.im*(*u3dn).c3_1.re+(*u2dn).c2_3.re*(*u3dn).c3_1.im;
   chi1.c3.re=
      (*u2dn).c3_1.re*(*u3dn).c1_1.re-(*u2dn).c3_1.im*(*u3dn).c1_1.im+
      (*u2dn).c3_2.re*(*u3dn).c2_1.re-(*u2dn).c3_2.im*(*u3dn).c2_1.im+ 
      (*u2dn).c3_3.re*(*u3dn).c3_1.re-(*u2dn).c3_3.im*(*u3dn).c3_1.im;
   chi1.c3.im=
      (*u2dn).c3_1.im*(*u3dn).c1_1.re+(*u2dn).c3_1.re*(*u3dn).c1_1.im+ 
      (*u2dn).c3_2.im*(*u3dn).c2_1.re+(*u2dn).c3_2.re*(*u3dn).c2_1.im+
      (*u2dn).c3_3.im*(*u3dn).c3_1.re+(*u2dn).c3_3.re*(*u3dn).c3_1.im;

   chi2.c1.re=
      (*u2dn).c1_1.re*(*u3dn).c1_2.re-(*u2dn).c1_1.im*(*u3dn).c1_2.im+
      (*u2dn).c1_2.re*(*u3dn).c2_2.re-(*u2dn).c1_2.im*(*u3dn).c2_2.im+ 
      (*u2dn).c1_3.re*(*u3dn).c3_2.re-(*u2dn).c1_3.im*(*u3dn).c3_2.im;
   chi2.c1.im=
      (*u2dn).c1_1.im*(*u3dn).c1_2.re+(*u2dn).c1_1.re*(*u3dn).c1_2.im+
      (*u2dn).c1_2.im*(*u3dn).c2_2.re+(*u2dn).c1_2.re*(*u3dn).c2_2.im+ 
      (*u2dn).c1_3.im*(*u3dn).c3_2.re+(*u2dn).c1_3.re*(*u3dn).c3_2.im;
   chi2.c2.re=
      (*u2dn).c2_1.re*(*u3dn).c1_2.re-(*u2dn).c2_1.im*(*u3dn).c1_2.im+ 
      (*u2dn).c2_2.re*(*u3dn).c2_2.re-(*u2dn).c2_2.im*(*u3dn).c2_2.im+
      (*u2dn).c2_3.re*(*u3dn).c3_2.re-(*u2dn).c2_3.im*(*u3dn).c3_2.im;
   chi2.c2.im=
      (*u2dn).c2_1.im*(*u3dn).c1_2.re+(*u2dn).c2_1.re*(*u3dn).c1_2.im+
      (*u2dn).c2_2.im*(*u3dn).c2_2.re+(*u2dn).c2_2.re*(*u3dn).c2_2.im+
      (*u2dn).c2_3.im*(*u3dn).c3_2.re+(*u2dn).c2_3.re*(*u3dn).c3_2.im;
   chi2.c3.re=
      (*u2dn).c3_1.re*(*u3dn).c1_2.re-(*u2dn).c3_1.im*(*u3dn).c1_2.im+
      (*u2dn).c3_2.re*(*u3dn).c2_2.re-(*u2dn).c3_2.im*(*u3dn).c2_2.im+ 
      (*u2dn).c3_3.re*(*u3dn).c3_2.re-(*u2dn).c3_3.im*(*u3dn).c3_2.im;
   chi2.c3.im=
      (*u2dn).c3_1.im*(*u3dn).c1_2.re+(*u2dn).c3_1.re*(*u3dn).c1_2.im+ 
      (*u2dn).c3_2.im*(*u3dn).c2_2.re+(*u2dn).c3_2.re*(*u3dn).c2_2.im+
      (*u2dn).c3_3.im*(*u3dn).c3_2.re+(*u2dn).c3_3.re*(*u3dn).c3_2.im;

   _suNg_inverse_multiply(psi1,(*u1dn),chi1);
   _suNg_inverse_multiply(psi2,(*u1dn),chi2);      
}
*/


static suNg staple, tr1, tr2;

static void add_to_v(suNg *v){
  _suNg_add_assign(*v,staple);
}

static void up_staple(void)
{
  _suNg_times_suNg(tr2,*u1up,*u2up);
  _suNg_dagger(tr1,*u3up);
  _suNg_times_suNg(staple,tr2,tr1);
}


static void dn_staple(void)
{
  _suNg_times_suNg(tr2,*u2dn,*u3dn);
  _suNg_dagger(tr1,*u1dn);
  _suNg_times_suNg(staple,tr1,tr2);
}
   

void staples(int ix,int mu,suNg *v)
{
   int i,nu,ixpmu,ixpnu,ixmnu,ixpmumnu;

   _suNg_zero(*v);

   ixpmu=iup[ix][mu];

   for (i=1;i<4;i++)
   {
      nu=(mu+i)&0x3;
      ixpnu=iup[ix][nu];
      ixmnu=idn[ix][nu];
      ixpmumnu=idn[ixpmu][nu];
      
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
  int ix,mu;
  suNg s, res;
  double pa=0.0;
  double ps=0.0, pl=0.;
  int fl;

  fl = 0;
  pa = avr_plaquette();
  for (ix=0;ix<VOLUME;++ix){
    for (mu=0; mu<4; ++mu){
      staples(ix, mu, &res);
      /* _suNg_dagger(s, res); */
      s=res;
      _suNg_times_suNg_dagger(res, *(pu_gauge(ix, mu)), s);
      ps += _suNg_trace_re(res);
    }
    pl+=local_plaq(ix);
  }

  ps /= 4.0*(double)(6*VOLUME*NG);
  pl/=(double)(6*VOLUME*NG);
	lprintf("TESTING",50,"Staple test: ");
  if (fabs(pa-ps)<1.e-6)
		lprintf("TESTING",50,"PASSED.");
  else
		lprintf("TESTING",50,"FAILED.");
	lprintf("TESTING",50," [diff1 = %1.8e][diff2 = %1.8e]\n", pa-ps, pl-ps);
}
