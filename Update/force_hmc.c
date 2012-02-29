/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *   
* All rights reserved.                                                      * 
\***************************************************************************/

#include "global.h"
#include "update.h"
#include "suN.h"
#include "utils.h"
#include "dirac.h"
#include "inverters.h"
#include "rational_functions.h"
#include "representation.h"
#include "logger.h"
#include "linear_algebra.h"
#include "memory.h"
#include "communications.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

/* declared in update_rhmc.c */
extern rhmc_par _update_par;

#define _print_avect(a) printf("(%3.5e,%3.5e,%3.5e,%3.5e,%3.5e,%3.5e,%3.5e,%3.5e)\n",(a).c1,(a).c2,(a).c3,(a).c4,(a).c5,(a).c6,(a).c7,(a).c8)

#define _print_mat(a) printf("(%3.5f,%3.5f,%3.5f)\n(%3.5f,%3.5f,%3.5f)\n(%3.5f,%3.5f,%3.5f)\n",(a).c1_1.re,(a).c1_2.re,(a).c1_3.re,(a).c2_1.re,(a).c2_2.re,(a).c2_3.re,(a).c3_1.re,(a).c3_2.re,(a).c3_3.re);printf("(%3.5f,%3.5f,%3.5f)\n(%3.5f,%3.5f,%3.5f)\n(%3.5f,%3.5f,%3.5f)\n",(a).c1_1.im,(a).c1_2.im,(a).c1_3.im,(a).c2_1.im,(a).c2_2.im,(a).c2_3.im,(a).c3_1.im,(a).c3_2.im,(a).c3_3.im)

/* we need to compute  Tr  U(x,mu) g_5*(1-g_mu) chi2 # chi1^+
 * where # indicates the tensor product and Tr is the trace on Lorentz space.
 * the strategy is the following:
 * given the form of g_5(1-g_mu) one can compute only the first two lorentz
 * components of the spinor; so we first apply g_5(1-g_mu) to chi2 to find the first
 * two components; then we multiply these two vectors by U(x,mu) and
 * store the result in p.c[0], p.c[1]; when computing the trace we can factorize p.c[0] and p.c[1]
 * as they both multiply two components of chi1^+; we store these factors in p.c[2] and p.c[3].
 * the tensor product is performed by the macro 
 * _suNf_FMAT(u,p): u = p.c[0] # p.c[2]^+ + p.c[1] # p.c[3]^+
 */

/* these macros use the variables ptmp, p */


#define _F_DIR0(u,chi1,chi2)				      \
  _vector_add_f(ptmp,(chi2)->c[0],(chi2)->c[2]);		      \
  _suNf_multiply(p.c[0],*(pu_gauge_f(x,0)),ptmp);		      \
  _vector_add_f(ptmp,(chi2)->c[1],(chi2)->c[3]);		      \
  _suNf_multiply(p.c[1],*(pu_gauge_f(x,0)),ptmp);		      \
  _vector_sub_f(p.c[2],(chi1)->c[0],(chi1)->c[2]);	      \
  _vector_sub_f(p.c[3],(chi1)->c[1],(chi1)->c[3]);	      \
  _suNf_FMAT((u),p)

#define _F_DIR1(u,chi1,chi2)				      \
  _vector_i_add_f(ptmp,(chi2)->c[0],(chi2)->c[3]);		      \
  _suNf_multiply(p.c[0],*(pu_gauge_f(x,1)),ptmp);		      \
  _vector_i_add_f(ptmp,(chi2)->c[1],(chi2)->c[2]);		      \
  _suNf_multiply(p.c[1],*(pu_gauge_f(x,1)),ptmp);		      \
  _vector_i_sub_f(p.c[2],(chi1)->c[0],(chi1)->c[3]);	      \
  _vector_i_sub_f(p.c[3],(chi1)->c[1],(chi1)->c[2]);	      \
  _suNf_FMAT((u),p)

#define _F_DIR2(u,chi1,chi2)				      \
  _vector_add_f(ptmp,(chi2)->c[0],(chi2)->c[3]);		      \
  _suNf_multiply(p.c[0],*(pu_gauge_f(x,2)),ptmp);		      \
  _vector_sub_f(ptmp,(chi2)->c[1],(chi2)->c[2]);		      \
  _suNf_multiply(p.c[1],*(pu_gauge_f(x,2)),ptmp);		      \
  _vector_sub_f(p.c[2],(chi1)->c[0],(chi1)->c[3]);	      \
  _vector_add_f(p.c[3],(chi1)->c[1],(chi1)->c[2]);	      \
  _suNf_FMAT((u),p)

#define _F_DIR3(u,chi1,chi2)				      \
  _vector_i_add_f(ptmp,(chi2)->c[0],(chi2)->c[2]);		      \
  _suNf_multiply(p.c[0],*(pu_gauge_f(x,3)),ptmp);		      \
  _vector_i_sub_f(ptmp,(chi2)->c[1],(chi2)->c[3]);		      \
  _suNf_multiply(p.c[1],*(pu_gauge_f(x,3)),ptmp);		      \
  _vector_i_sub_f(p.c[2],(chi1)->c[0],(chi1)->c[2]);	      \
  _vector_i_add_f(p.c[3],(chi1)->c[1],(chi1)->c[3]);	      \
  _suNf_FMAT((u),p)



/*
typedef struct {
  double mass;
  double has_a;
  double has_b;
  float mass_flt;
  float has_a_flt;
  float has_b_flt;
} dirac_par;
*/

void D(spinor_field *out, spinor_field *in){
#ifdef UPDATE_EO
    Dphi_eopre(_update_par.mass, out, in);
#else
    Dphi(_update_par.mass, out, in);
#endif
}


void D_flt(spinor_field_flt *out, spinor_field_flt *in){
#ifdef UPDATE_EO
    Dphi_eopre_flt((float)(_update_par.mass), out, in);
#else
    Dphi_flt((float)(_update_par.mass), out, in);
#endif
}


void Force_hmc_f(double dt, suNg_av_field *force, void *vpar){
  _DECLARE_INT_ITERATOR(x);
  int mu,  k;
  static suNg_algebra_vector f;
  static suNf_vector ptmp;
  static suNf_spinor p;
  static suNf s1;
  static MINRES_par inv_par;
  spinor_field *X, *Y, *eta;
  spinor_field_flt *eta_flt;
#ifdef UPDATE_EO
  spinor_field Xe, Xo, Ye, Yo;
#endif
  double forcestat[2]; /* used for computation of avr and max force */
  double nsq;

  force_hmc_par *par = (force_hmc_par*)vpar;
  spinor_field *pf = par->pf;
  
  /* check input types */
#ifndef CHECK_SPINOR_MATCHING
   _TWO_SPINORS_MATCHING(u_gauge,force);
#endif

  /* allocate spinors */
#ifndef UPDATE_EO
  X = alloc_spinor_field_f(3,&glattice);
  Y = X+1;
  eta = Y+1;
  eta_flt = alloc_spinor_field_f_flt(2,&glattice);
#else
  X = alloc_spinor_field_f(2,&glattice);
  Y = X+1;
  eta = alloc_spinor_field_f(1,&glat_even);
  eta_flt = alloc_spinor_field_f_flt(2,&glat_even);
#endif

  inv_par.max_iter=0;

  assign_ud2u_f();
  
  for (k=0; k<_update_par.n_pf; ++k) {

#ifndef UPDATE_EO
    /* X = H^{-1} pf[k] = D^{-1} g5 pf[k] */
    g5QMR_fltacc_par mpar;
    mpar.err2 = _update_par.force_prec;
    mpar.max_iter = 0;
    mpar.err2_flt = _update_par.force_prec_flt;
    mpar.max_iter_flt = 0;
    spinor_field_zero_f(X);
    spinor_field_g5_assign_f(&pf[k]);
    g5QMR_fltacc(&mpar, &D, &D_flt, &pf[k], X);
    spinor_field_g5_assign_f(&pf[k]);

    /*
    spinor_field_zero_f_flt(&eta_flt[0]);
    assign_sd2s(&eta_flt[1],&pf[k]);
    inv_par.err2 = _update_par.force_prec_flt;
    HBiCGstab_flt(&inv_par, &H_flt, &eta_flt[1], &eta_flt[0]);
    assign_s2sd(X,&eta_flt[0]);
    inv_par.err2= _update_par.force_prec;
    HBiCGstab(&inv_par, &H, &pf[k], X);
    */
        
    /* Y = H^{-1} ( a g5 pf[k] + b X ) = D^{-1} ( a pf[k] + b g5 X ) */
    if(par->hasenbusch != 2) {
      spinor_field_g5_f(eta,X);
    } else {
      spinor_field_g5_f(eta,X);
      spinor_field_mul_f(eta,par->bY,eta);
      spinor_field_mul_add_assign_f(eta,par->aY,&pf[k]);
    }
    spinor_field_zero_f(Y);
    g5QMR_fltacc(&mpar, &D, &D_flt, eta, Y);

    /*
    if(par->hasenbusch != 2) {
      spinor_field_copy_f(eta,X);
    } else {
      spinor_field_g5_f(eta,&pf[k]);
      spinor_field_mul_f(eta,par->aY,eta);
      spinor_field_mul_add_assign_f(eta,par->bY,X);
    }
    spinor_field_zero_f_flt(&eta_flt[0]);
    assign_sd2s(&eta_flt[1],eta);
    inv_par.err2 = _update_par.force_prec_flt;
    HBiCGstab_flt(&inv_par, &H_flt, &eta_flt[1], &eta_flt[0]);
    assign_s2sd(Y,&eta_flt[0]);
    inv_par.err2= _update_par.force_prec;
    HBiCGstab(&inv_par, &H, eta, Y);
    */
    
#else
    /* X_e = H^{-1} pf[k] */
    /* X_o = D_{oe} X_e = D_{oe} H^{-1} pf[k] */
    Xe=*X; Xe.type=&glat_even;
    Xo=*X; Xo.type=&glat_odd;

    /*
    spinor_field_zero_f_flt(&eta_flt[0]);
    assign_sd2s(&eta_flt[1],&pf[k]);
    inv_par.err2 = _update_par.force_prec_flt;
    HBiCGstab_flt(&inv_par, &H_flt, &eta_flt[1], &eta_flt[0]);
    assign_s2sd(&Xe,&eta_flt[0]);
    inv_par.err2= _update_par.force_prec;
    HBiCGstab(&inv_par, &H, &pf[k], &Xe);
    Dphi_(&Xo,&Xe);
    */

    g5QMR_fltacc_par mpar;
    mpar.err2 = _update_par.force_prec;
    mpar.max_iter = 0;
    mpar.err2_flt = _update_par.force_prec_flt;
    mpar.max_iter_flt = 0;
    spinor_field_zero_f(&Xe);
    /* H^{-1} pf = D^{-1} g5 pf */
    spinor_field_g5_assign_f(&pf[k]);
    g5QMR_fltacc(&mpar, &D, &D_flt, &pf[k], &Xe);
    spinor_field_g5_assign_f(&pf[k]);
    Dphi_(&Xo,&Xe);
    
    /* Y_e = H^{-1} ( a g5 pf[k] + b X_e ) */
    /* Y_o = D_oe H^{-1} ( a g5 pf[k] + b X_e ) */
    Ye=*Y; Ye.type=&glat_even;
    Yo=*Y; Yo.type=&glat_odd;

    if(par->hasenbusch != 2) {
      spinor_field_copy_f(eta,&Xe);
    } else {
      spinor_field_g5_f(eta,&pf[k]);
      spinor_field_mul_f(eta,par->aY,eta);
      spinor_field_mul_add_assign_f(eta,par->bY,&Xe);
    }
    
    /*
    spinor_field_zero_f_flt(&eta_flt[0]);
    assign_sd2s(&eta_flt[1],eta);
    inv_par.err2 = _update_par.force_prec_flt;
    HBiCGstab_flt(&inv_par, &H_flt, &eta_flt[1], &eta_flt[0]);
    assign_s2sd(&Ye,&eta_flt[0]);
    inv_par.err2= _update_par.force_prec;
    HBiCGstab(&inv_par, &H, eta, &Ye);
    Dphi_(&Yo,&Ye);
    */

    spinor_field_zero_f(&Ye);
    spinor_field_g5_assign_f(eta);
    g5QMR_fltacc(&mpar, &D, &D_flt, eta, &Ye);
    spinor_field_g5_assign_f(eta);
    Dphi_(&Yo,&Ye);
    
#endif


    lprintf("FORCE",50,"|X| = %1.8e |Y| = %1.8e\n",
	    sqrt(spinor_field_sqnorm_f(X)),
	    sqrt(spinor_field_sqnorm_f(Y))
	    );
    
    /* reset force stat counters */
    start_sf_sendrecv(X);
    start_sf_sendrecv(Y);
    
    forcestat[1]=forcestat[0]=0.;
    
    _PIECE_FOR(&glattice,x) { 
      _SITE_FOR(&glattice,x) {
	
      	for (mu=0; mu<4; ++mu) {
      	  int y;
      	  suNf_spinor *chi1, *chi2;
      	  _suNf_zero(s1);
      	  switch (mu) {
      	  case 0:
      	    y=iup(x,0);
      	    chi1=_FIELD_AT(X,x);
      	    chi2=_FIELD_AT(Y,y);
      	    _F_DIR0(s1,chi1,chi2);
      	    chi1=_FIELD_AT(Y,x);
      	    chi2=_FIELD_AT(X,y);
      	    _F_DIR0(s1,chi1,chi2);
      	    break;
      	  case 1:
      	    y=iup(x,1);
      	    chi1=_FIELD_AT(X,x);
      	    chi2=_FIELD_AT(Y,y);
      	    _F_DIR1(s1,chi1,chi2);
      	    chi1=_FIELD_AT(Y,x);
      	    chi2=_FIELD_AT(X,y);
      	    _F_DIR1(s1,chi1,chi2);
      	    break;
      	  case 2:
      	    y=iup(x,2);
      	    chi1=_FIELD_AT(X,x);
      	    chi2=_FIELD_AT(Y,y);
      	    _F_DIR2(s1,chi1,chi2);
      	    chi1=_FIELD_AT(Y,x);
      	    chi2=_FIELD_AT(X,y);
      	    _F_DIR2(s1,chi1,chi2);
      	    break;
      	  default: /* DIR 3 */
      	    y=iup(x,3);
      	    chi1=_FIELD_AT(X,x);
      	    chi2=_FIELD_AT(Y,y);
      	    _F_DIR3(s1,chi1,chi2);
      	    chi1=_FIELD_AT(Y,x);
      	    chi2=_FIELD_AT(X,y);
      	    _F_DIR3(s1,chi1,chi2);
      	  }
	  
      	  _algebra_project(f,s1);


#ifdef UPDATE_EO
          if(par->hasenbusch != 1) {
      	    _algebra_vector_mul_add_assign_g(*_4FIELD_AT(force,x,mu),-dt*(_REPR_NORM2/_FUND_NORM2),f);
      	  } else {
      	    _algebra_vector_mul_add_assign_g(*_4FIELD_AT(force,x,mu),-par->bD*dt*(_REPR_NORM2/_FUND_NORM2),f);
          }
#else
          if(par->hasenbusch != 1) {
        	  _algebra_vector_mul_add_assign_g(*_4FIELD_AT(force,x,mu),dt*(_REPR_NORM2/_FUND_NORM2),f);	
      	  } else {
      	    _algebra_vector_mul_add_assign_g(*_4FIELD_AT(force,x,mu),par->bD*dt*(_REPR_NORM2/_FUND_NORM2),f);
          }
#endif
	  
    	    _algebra_vector_sqnorm_g(nsq,f);
    	    forcestat[0]+=sqrt(nsq);
    	    for(y=0;y<NG*NG-1;++y){
    	      if(forcestat[1]<fabs(*(((double*)&f)+y))) forcestat[1]=fabs(*(((double*)&f)+y));
    	    }
      	}
      }
      if(_PIECE_INDEX(x)==0) {
      	complete_sf_sendrecv(X);
       	complete_sf_sendrecv(Y);
      }
    }

    global_sum(forcestat,2);
    forcestat[0]*=dt*(_REPR_NORM2/_FUND_NORM2)/((double)(4*GLB_T*GLB_X*GLB_Y*GLB_Z));
    forcestat[1]*=dt*(_REPR_NORM2/_FUND_NORM2);
    lprintf("FORCE",50," avr |force| = %1.8e maxforce = %1.8e \n",forcestat[0],forcestat[1]);
    
#if defined(BASIC_SF) || defined(ROTATED_SF)
    SF_force_bcs(force);
#endif /* BASIC_SF || ROTATED_SF*/
    
  }
  
  free_spinor_field_flt(eta_flt);
  free_spinor_field(X);
#ifdef UPDATE_EO
  free_spinor_field(eta);
#endif

}

#undef _F_DIR0
#undef _F_DIR1
#undef _F_DIR2
#undef _F_DIR3
