/*******************************************************************************
*
* File plaquette.c
*
* Routines for the average plaquette
*
*******************************************************************************/

#include "global.h"
#include "suN.h"

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
   /*   _suNg_dagger(w3,w2); */
   _suNg_times_suNg_dagger(w3,w1,w2);      
      
   _suNg_trace_re(p,w3);
	 return p;
}

double avr_plaquette()
{
  int ix, mu, nu;
  double pa=0.0;
  
	for (ix=0;ix<VOLUME;ix++)
		for (mu=1;mu<4;mu++)
			for (nu=0;nu<mu;nu++)
				pa+=(double)(plaq(ix,mu,nu));

  return pa/(double)(6*VOLUME*NG);

}

double local_plaq(int ix)
{
  int mu, nu;
  double pa=0.0;
  
  for (mu=1;mu<4;++mu)
     for (nu=0;nu<mu;++nu)
			 pa+=(double)(plaq(ix,mu,nu));
  
  return pa;

}
