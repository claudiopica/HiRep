/***************************************************************************\
* Copyright (c) 2008, Claudio Pica, Vincent Drach and Ari Hietanen          *   
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

#define P99_PROTECT(...) __VA_ARGS__ 

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

#ifdef BC_T_THETA
#define _T_theta_mulc(r) _vector_mulc_f(ptmp,eitheta[0],(r)); (r)=ptmp
#else
#define _T_theta_mulc(r)
#endif
#ifdef BC_X_THETA
#define _X_theta_mulc(r) _vector_mulc_f(ptmp,eitheta[1],(r)); (r)=ptmp
#else
#define _X_theta_mulc(r)
#endif
#ifdef BC_Y_THETA
#define _Y_theta_mulc(r) _vector_mulc_f(ptmp,eitheta[2],(r)); (r)=ptmp
#else
#define _Y_theta_mulc(r)
#endif
#ifdef BC_Z_THETA
#define _Z_theta_mulc(r) _vector_mulc_f(ptmp,eitheta[3],(r)); (r)=ptmp
#else
#define _Z_theta_mulc(r)
#endif

#define _F_DIR0(u,chi1,chi2)				      \
  _vector_add_f(ptmp,(chi2)->c[0],(chi2)->c[2]);		      \
  _suNf_multiply(p.c[0],*(pu_gauge_f(x,0)),ptmp);		      \
  _T_theta_mulc(p.c[0]);                                      \
  _vector_add_f(ptmp,(chi2)->c[1],(chi2)->c[3]);		      \
  _suNf_multiply(p.c[1],*(pu_gauge_f(x,0)),ptmp);		      \
  _T_theta_mulc(p.c[1]);                                      \
  _vector_sub_f(p.c[2],(chi1)->c[0],(chi1)->c[2]);	      \
  _vector_sub_f(p.c[3],(chi1)->c[1],(chi1)->c[3]);	      \
  _suNf_FMAT((u),p)

#define _F_DIR1(u,chi1,chi2)				      \
  _vector_i_add_f(ptmp,(chi2)->c[0],(chi2)->c[3]);		      \
  _suNf_multiply(p.c[0],*(pu_gauge_f(x,1)),ptmp);		      \
  _X_theta_mulc(p.c[0]);                                      \
  _vector_i_add_f(ptmp,(chi2)->c[1],(chi2)->c[2]);		      \
  _suNf_multiply(p.c[1],*(pu_gauge_f(x,1)),ptmp);		      \
  _X_theta_mulc(p.c[1]);                                      \
  _vector_i_sub_f(p.c[2],(chi1)->c[0],(chi1)->c[3]);	      \
  _vector_i_sub_f(p.c[3],(chi1)->c[1],(chi1)->c[2]);	      \
  _suNf_FMAT((u),p)

#define _F_DIR2(u,chi1,chi2)				      \
  _vector_add_f(ptmp,(chi2)->c[0],(chi2)->c[3]);		      \
  _suNf_multiply(p.c[0],*(pu_gauge_f(x,2)),ptmp);		      \
  _Y_theta_mulc(p.c[0]);                                      \
  _vector_sub_f(ptmp,(chi2)->c[1],(chi2)->c[2]);		      \
  _suNf_multiply(p.c[1],*(pu_gauge_f(x,2)),ptmp);		      \
  _Y_theta_mulc(p.c[1]);                                      \
  _vector_sub_f(p.c[2],(chi1)->c[0],(chi1)->c[3]);	      \
  _vector_add_f(p.c[3],(chi1)->c[1],(chi1)->c[2]);	      \
  _suNf_FMAT((u),p)

#define _F_DIR3(u,chi1,chi2)				      \
  _vector_i_add_f(ptmp,(chi2)->c[0],(chi2)->c[2]);		      \
  _suNf_multiply(p.c[0],*(pu_gauge_f(x,3)),ptmp);		      \
  _Z_theta_mulc(p.c[0]);                                      \
  _vector_i_sub_f(ptmp,(chi2)->c[1],(chi2)->c[3]);		      \
  _suNf_multiply(p.c[1],*(pu_gauge_f(x,3)),ptmp);		      \
  _Z_theta_mulc(p.c[1]);                                      \
  _vector_i_sub_f(p.c[2],(chi1)->c[0],(chi1)->c[2]);	      \
  _vector_i_add_f(p.c[3],(chi1)->c[1],(chi1)->c[3]);	      \
  _suNf_FMAT((u),p)


suNg_av_field *force_tmp=NULL;


void force_fermion_core(spinor_field* Xs, spinor_field* Ys, suNg_av_field* force, double dt, double* forcestat, int type){
#ifdef MEASURE_FORCEHMC
  lprintf("FORCE",10,"|Xs| = %1.8e |Ys| = %1.8e\n",
          sqrt(spinor_field_sqnorm_f(Xs)),
          sqrt(spinor_field_sqnorm_f(Ys))
          );
  if (type==1 && force_tmp==NULL){
    force_tmp = alloc_avfield(&glattice);
  }
#endif

    /* reset force stat counters */
    start_sf_sendrecv(Xs);
    start_sf_sendrecv(Ys);

    _PIECE_FOR(&glattice,xp) {
      suNg_algebra_vector f;
      suNf_vector ptmp;
      suNf_spinor p;
      suNf_FMAT s1;

      if (xp==glattice.inner_master_pieces) {
        _OMP_PRAGMA( master )
        {
          complete_sf_sendrecv(Xs);
          complete_sf_sendrecv(Ys);
        }
        _OMP_PRAGMA( barrier )
      }

#ifdef MEASURE_FORCEHMC
      //      _SITE_FOR_SUM(&glattice,xp,x,forcestat[0],forcestat[1]) {
      _SITE_FOR_SUM(&glattice,xp,x,forcestat[0],forcestat[1]) {
#else
      _SITE_FOR(&glattice,xp,x) {
#endif
      	for (int mu=0; mu<4; ++mu) {
      	  int y;
      	  suNf_spinor *chi1, *chi2;
      	  _suNf_FMAT_zero(s1);
      	  switch (mu) {
            case 0:
              y=iup(x,0);
              chi1=_FIELD_AT(Xs,x);
              chi2=_FIELD_AT(Ys,y);
              _F_DIR0(s1,chi1,chi2);
              chi1=_FIELD_AT(Ys,x);
              chi2=_FIELD_AT(Xs,y);
              _F_DIR0(s1,chi1,chi2);
              break;
            case 1:
              y=iup(x,1);
              chi1=_FIELD_AT(Xs,x);
              chi2=_FIELD_AT(Ys,y);
              _F_DIR1(s1,chi1,chi2);
              chi1=_FIELD_AT(Ys,x);
              chi2=_FIELD_AT(Xs,y);
              _F_DIR1(s1,chi1,chi2);
              break;
            case 2:
              y=iup(x,2);
              chi1=_FIELD_AT(Xs,x);
              chi2=_FIELD_AT(Ys,y);
              _F_DIR2(s1,chi1,chi2);
              chi1=_FIELD_AT(Ys,x);
              chi2=_FIELD_AT(Xs,y);
              _F_DIR2(s1,chi1,chi2);
              break;
            default: /* DIR 3 */
              y=iup(x,3);
              chi1=_FIELD_AT(Xs,x);
              chi2=_FIELD_AT(Ys,y);
              _F_DIR3(s1,chi1,chi2);
              chi1=_FIELD_AT(Ys,x);
              chi2=_FIELD_AT(Xs,y);
              _F_DIR3(s1,chi1,chi2);
      	  }
          _algebra_project(f,s1);
          _algebra_vector_mul_add_assign_g(*_4FIELD_AT(force,x,mu),dt*(_REPR_NORM2/_FUND_NORM2),f);

#ifdef MEASURE_FORCEHMC
          if (type==1){
            *_4FIELD_AT(force_tmp,x,mu)=f;
          }
          else{
            if (type==2){
              _algebra_vector_sub_assign_g(f,*_4FIELD_AT(force_tmp,x,mu));
            }
            double nsq;
            _algebra_vector_sqnorm_g(nsq,f);
            nsq=sqrt(nsq);
            forcestat[0]+=nsq;
            if (nsq>forcestat[1]) forcestat[1]=nsq;
          }
#endif
      	} //directions for
      } //SITE_FOR
      } //PIECE FOR
      

#if defined(BASIC_SF) || defined(ROTATED_SF)
    SF_force_bcs(force);
#endif /* BASIC_SF || ROTATED_SF*/

}

#undef _F_DIR0
#undef _F_DIR1
#undef _F_DIR2
#undef _F_DIR3
