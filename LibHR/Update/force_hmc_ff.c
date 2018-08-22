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



//Temporary storage fields from force_hmc.c
extern spinor_field *Xs, *Ys, *eta;
extern spinor_field_flt *eta_flt;

void init_force_hmc();
void free_force_hmc();

void spinor_sigma_pi_rho_div_assign(spinor_field *out,scalar_field *sigma,scalar_field *pi,double rho, spinor_field *in);
void spinor_sigma_pi_dagger_rho_div_assign(spinor_field *out,scalar_field *sigma,scalar_field *pi,double rho, spinor_field *in);
void force_hmc_auxfields_fermion(double dt, void *vpar, scalar_field *sigma_mom, scalar_field *pi_mom,spinor_field *Xs, spinor_field *Ys, int hasenbusch );

void force_hmc_ff(double dt, void *vpar){

#ifdef UPDATE_EO
  spinor_field Xe, Xo, Ye, Yo;
#endif

  init_force_hmc();
  force_hmc_par *par = (force_hmc_par*)vpar;
  suNg_av_field *force = *par->momenta;
  spinor_field *pf = par->pf;
  int n_iters = 0;
  (void) n_iters;
  double mass = par->mass;
  set_ff_dirac_mass(par->mass);

  /* check input types */
  _TWO_SPINORS_MATCHING(u_gauge,force);

  for (int k=0; k<par->n_pf; ++k) {

#ifndef UPDATE_EO
	  // Not implemented
#else

    /*    g5QMR_fltacc_par mpar;
    mpar.err2 = par->inv_err2;
    mpar.max_iter = 0;
    mpar.err2_flt = par->inv_err2_flt;
    mpar.max_iter_flt = 0;*/
    double tmp;
    mshift_par mpar;
    mpar.err2 = par->inv_err2;
    mpar.max_iter=0;
    mpar.n = 1;
    mpar.shift = &tmp; // mpar.n = 1
    mpar.shift[0] = 0;

    double rho = 4.+mass;
    spinor_field *tmp_spinor_field = alloc_spinor_field_f(1,&glattice);


    /* Y_e = (D^ D)^{-1} pf[k]  */
    Ye=*Ys; Ye.type=&glat_even;
    Yo=*Ys; Yo.type=&glat_odd; Yo.ptr+=glat_odd.master_shift;

    /* X_e = (D^)^(-1) pf[k] = g5 D Y_e */
    Xe=*Xs; Xe.type=&glat_even;
    Xo=*Xs; Xo.type=&glat_odd; Xo.ptr+=glat_odd.master_shift;

    set_ff_dirac_mass(par->mass);

    if(par->hasenbusch == 2) {
      /*  Y_e = (D^ D)^{-1} (D+b) pf[k] */
      set_ff_dirac_shift(par->b);
      Dff_dagger(&Xe,&pf[k]);
      set_ff_dirac_shift(0.);
    } else {
      spinor_field_copy_f(&Xe,&pf[k]);
    }

    lprintf("FORCE",0," id %d mass %e shift %e \n",
	   par->id, par->mass, par->b );


    spinor_field_zero_f(&Ye);
    mre_guess( &par->mpar, 0, &Ye, &Dff_sq, &Xe);
    n_iters+=cg_mshift( &mpar, &Dff_sq, &Xe, &Ye );
    mre_store( &par->mpar, 0, &Ye);
    
    Dff(&Xe,&Ye);

    if(par->hasenbusch == 2) {
      /*  X_e = D (D^ D)^{-1} (D+b) pf[k] - pf[k] */
      spinor_field_mul_add_assign_f(&Xe,-1,&pf[k]);
    }
    
    spinor_field_g5_assign_f(&Xe);

    /* Y_o = A_o^{-1} D_oe (D^ D)^{-1} pf[k]  */
    Dphi_(&Yo,&Ye);
    spinor_sigma_pi_rho_div_assign(&Yo,ff_sigma,ff_pi,rho, &Yo);

    /* X_o = (A_o)^^{-1} D_{oe} X_e = (A_o^)^{-1} D_{oe} g5 (D^)^{-1} pf[k] */
    Dphi_(&Xo,&Xe);
    spinor_sigma_pi_dagger_rho_div_assign(&Xo,ff_sigma,ff_pi,rho, &Xo);
    
    //Add the force of the fermion fields on the auxiliary fields,
    //from the derivative of the Dirac operator.
    force_hmc_auxfields_fermion( dt, vpar, ff_sigma_mom, ff_pi_mom, Xs, Ys, par->hasenbusch);
    
    free_spinor_field_f(tmp_spinor_field);
        
#endif

    
#ifdef MEASURE_FORCEHMC
    lprintf("FORCE",50,"|Xs| = %1.8e |Ys| = %1.8e\n",
	    sqrt(spinor_field_sqnorm_f(Xs)),
	    sqrt(spinor_field_sqnorm_f(Ys))
	    );
    double forcestat[2]={0.,0.}; /* used for computation of avr and max force */
#endif

    if( gauge_field_active ){
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

      	  _algebra_project_FMAT(f,s1);


#ifdef UPDATE_EO
          if(par->hasenbusch != 2) {
      	    _algebra_vector_mul_add_assign_g(*_4FIELD_AT(force,x,mu),-dt*(_REPR_NORM2/_FUND_NORM2),f);
      	  } else {
      	    _algebra_vector_mul_add_assign_g(*_4FIELD_AT(force,x,mu),-dt*(_REPR_NORM2/_FUND_NORM2),f);
          }
#else
          if(par->hasenbusch != 2) {
            _algebra_vector_mul_add_assign_g(*_4FIELD_AT(force,x,mu),dt*(_REPR_NORM2/_FUND_NORM2),f);
      	  } else {
      	    _algebra_vector_mul_add_assign_g(*_4FIELD_AT(force,x,mu),dt*(_REPR_NORM2/_FUND_NORM2),f);
          }
#endif

#ifdef MEASURE_FORCEHMC
          double nsq;
	  _algebra_vector_sqnorm_g(nsq,f);
	  forcestat[0]+=sqrt(nsq);
	  for(y=0;y<NG*NG-1;++y){
	    if(forcestat[1]<fabs(*(((double*)&f)+y))) forcestat[1]=fabs(*(((double*)&f)+y));
	  }
#endif
      	} //directions for
      } //SITE_FOR
    } //PIECE FOR
    } //if

#ifdef MEASURE_FORCEHMC
    global_sum(forcestat,1);
    forcestat[0]*=dt*(_REPR_NORM2/_FUND_NORM2)/((double)(4*GLB_T*GLB_X*GLB_Y*GLB_Z));
    global_max(forcestat+1,1);
    forcestat[1]*=dt*(_REPR_NORM2/_FUND_NORM2);
    force_ave[par->id+1]+=forcestat[0];
    force_max[par->id+1]+=forcestat[1];
    n_inv_iter[par->id]+=n_iters;
    lprintf("FORCE_HMC",20,"avr dt |force| = %1.8e dt maxforce = %1.8e, dt = %1.8e \n",forcestat[0],forcestat[1],dt);
#endif

  }

}

#undef _F_DIR0
#undef _F_DIR1
#undef _F_DIR2
#undef _F_DIR3

