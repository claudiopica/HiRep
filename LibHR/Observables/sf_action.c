/*! \file
 * \brief Routines for the SF gauge observable
 *
 * Routines for the Schrodinger Functional gauge observable.
 *
 */

/*! \defgroup sfgauge SF Gauge Observable
 * \ingroup obs
 */
#ifdef GAUGE_SUN


#include "global.h"
#include "suN.h"
#include "error.h"
#include "communications.h"
#include "geometry.h"
#include "observables.h"
#include <stdio.h>

/*! \ingroup sfgauge
 * Bottom E_8 component
 *
 * \f$ E_k^8({\bf x})= {\frac{1}{a^2}}\mathrm{Re}\, \mathrm{tr}\, \left\{ i \lambda_8 U(x,k) U(x+a\hat k,0) U(x+a\hat 0,k)^{-1} U(x,0)^{-1} \right\}_{x_0=0}\f$
 *
 * The matrix \f$\lambda_8\f$ depends on the size of the gauge group, for SU(3) it is given by \f$\lambda_8={\rm diag}(1,-1/2,-1/2)\f$
 *
 */

static double E_8(int ix,int k)
{
   double p;
   suNg *v1,*v2,*v3,*v4,w1,w2,w3,w4,ilambda8;
   int iy, iz;

#if NG==4
   _suNg_zero(ilambda8);
   ((ilambda8).c[0]).re = 0.0;
   ((ilambda8).c[0]).im = 0.5;
   ((ilambda8).c[5]).re = 0.0;
   ((ilambda8).c[5]).im = 0.5;
   ((ilambda8).c[10]).re = 0.0;
   ((ilambda8).c[10]).im = -0.5;
   ((ilambda8).c[15]).re = 0.0;
   ((ilambda8).c[15]).im = -0.5;
#elif NG==3
   _suNg_zero(ilambda8);
   ((ilambda8).c[0]).re = 0.0;
   ((ilambda8).c[0]).im = 1.0;
   ((ilambda8).c[4]).re = 0.0;
   ((ilambda8).c[4]).im = -0.5;
   ((ilambda8).c[8]).re = 0.0;
   ((ilambda8).c[8]).im = -0.5;
#elif NG==2
   _suNg_zero(ilambda8);
#ifdef WITH_QUATERNIONS
  (ilambda8).c[3]=1.0;
#else
   ((ilambda8).c[0]).re = 0.0;
   ((ilambda8).c[0]).im = 1.0;
   ((ilambda8).c[3]).re = 0.0;
   ((ilambda8).c[3]).im = -1.0;
#endif
#else
#error  "No explicit generator for observable for given NG"
#endif 
   
   iy=iup(ix,k);
   iz=iup(ix,0);
   
   v1=pu_gauge(ix,k);
   v2=pu_gauge(iy,0);
   v3=pu_gauge(iz,k);
   v4=pu_gauge(ix,0);

   _suNg_times_suNg(w1,(*v1),(*v2));
   _suNg_times_suNg(w2,(*v4),(*v3));
   _suNg_times_suNg_dagger(w3,w1,w2);      
   _suNg_times_suNg(w4,ilambda8,w3);

   _suNg_trace_re(p,w4);

   return p;
}

/*! \ingroup sfgauge
 * Top E_8 component
 *
 * \f$ (E_k^8)'({\bf x})= -{\frac{1}{a^2}}\mathrm{Re}\, \mathrm{tr}\, \left\{ i \lambda_8 U(x+a\hat 0,k)^{-1} U(x,0)^{-1} U(x,k) U(x+a\hat k,0) \right\}_{x_0=T-2a}\f$
 *
 * The matrix \f$\lambda_8\f$ depends on the size of the gauge group, for SU(3) it is given by \f$\lambda_8={\rm diag}(1,-1/2,-1/2)\f$

 */

static double E_8_top(int ix,int k)
{
  double p;
  suNg *v1,*v2,*v3,*v4,w1,w2,w3,w4,ilambda8;
  int iy, iz;
#if NG==4
  _suNg_zero(ilambda8);
  ((ilambda8).c[0]).re = 0.0;
  ((ilambda8).c[0]).im = 0.5;
  ((ilambda8).c[5]).re = 0.0;
  ((ilambda8).c[5]).im = 0.5;
  ((ilambda8).c[10]).re = 0.0;
  ((ilambda8).c[10]).im = -0.5;
  ((ilambda8).c[15]).re = 0.0;
  ((ilambda8).c[15]).im = -0.5;
#elif NG==3
  _suNg_zero(ilambda8);
  ((ilambda8).c[0]).re = 0.0;
  ((ilambda8).c[0]).im = 1.0;
  ((ilambda8).c[4]).re = 0.0;
  ((ilambda8).c[4]).im = -0.5;
  ((ilambda8).c[8]).re = 0.0;
  ((ilambda8).c[8]).im = -0.5;
#elif NG==2
  _suNg_zero(ilambda8);
#ifdef WITH_QUATERNIONS
  (ilambda8).c[3]=1.0;
#else
  _suNg_zero(ilambda8);
  ((ilambda8).c[0]).re = 0.0;
  ((ilambda8).c[0]).im = 1.0;
  ((ilambda8).c[3]).re = 0.0;
  ((ilambda8).c[3]).im = -1.0;
#endif
#else
#error "No explicit generator for observable for given NG"
#endif
  iy=iup(ix,k);
  iz=iup(ix,0);
  
  v1=pu_gauge(ix,k);
  v2=pu_gauge(iy,0);
  v3=pu_gauge(iz,k);
  v4=pu_gauge(ix,0);
  
  _suNg_times_suNg(w1,(*v1),(*v2));
  _suNg_times_suNg(w2,ilambda8,(*v3));
  _suNg_times_suNg_dagger(w3,w2,w1);      
  _suNg_times_suNg(w4,(*v4),w3);

   _suNg_trace_re(p,w4);

   return -p;
}

/*! \ingroup sfgauge
 *  SF Gauge Observable
 *
 * \f$\frac{\partial S_g}{\partial\eta}= - {2 \over g_0^2 \, L} a^3 \sum_{\bf x}\left\{E_k^8({\bf x}) - (E_k^8)'({\bf x}) \right\}\f$
 *
 * where
 *
 * \f$ E_k^8({\bf x})= {\frac{1}{a^2}}\mathrm{Re}\, \mathrm{tr}\, \left\{ i \lambda_8 U(x,k) U(x+a\hat k,0) U(x+a\hat 0,k)^{-1} U(x,0)^{-1} \right\}_{x_0=0}\f$
 *
 * \f$ (E_k^8)'({\bf x})= -{\frac{1}{a^2}}\mathrm{Re}\, \mathrm{tr}\, \left\{ i \lambda_8 U(x+a\hat 0,k)^{-1} U(x,0)^{-1} U(x,k) U(x+a\hat k,0) \right\}_{x_0=T-2a}\f$
 *
 * The matrix \f$\lambda_8\f$ depends on the size of the gauge group, for SU(3) it is given by \f$\lambda_8={\rm diag}(1,-1/2,-1/2)\f$

 */

double SF_action(double beta)
{
  double pa=0.;
  int ix, iy, iz,index;
  complete_gf_sendrecv(u_gauge);
  if(COORD[0] == 0) {
    for (ix=0;ix<X;++ix) for (iy=0;iy<Y;++iy) for (iz=0;iz<Z;++iz){
      index=ipt(1,ix,iy,iz);
      pa+=(double)(E_8(index,1));
      pa+=(double)(E_8(index,2));
      pa+=(double)(E_8(index,3));
    }
  }
  if(COORD[0] == NP_T-1) {
    for (ix=0;ix<X;++ix) for (iy=0;iy<Y;++iy) for (iz=0;iz<Z;++iz){
      index=ipt(T-2,ix,iy,iz);
      pa+=(double)(E_8_top(index,1));
      pa+=(double)(E_8_top(index,2));
      pa+=(double)(E_8_top(index,3));
    }
  }
  global_sum(&pa, 1);
  return pa*(double)(beta/(NG*GLB_X));
}

#endif
