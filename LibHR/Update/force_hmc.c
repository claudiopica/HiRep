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


static double static_mass=0.;
static double static_mu=0.;



static spinor_field *Xs=NULL, *Ys=NULL, *eta=NULL;

//static spinor_field_flt *eta_flt=NULL;

#ifdef UPDATE_EO
static spinor_field *xi=NULL;
#endif 


void free_force_hmc();

void init_force_hmc() {
  static int init=0;
  if (!init){
#ifndef UPDATE_EO
    Xs = alloc_spinor_field_f(4,&glattice);
    Ys = Xs+1;
    eta = Ys+1;
#else
    Xs = alloc_spinor_field_f(2,&glattice);
    Ys = Xs+1;
    eta = alloc_spinor_field_f(1,&glat_even);
    xi = alloc_spinor_field_f(1,&glat_odd);
#endif
    init=1;
    atexit(free_force_hmc);
  }
}

void free_force_hmc() {
  free_spinor_field_f(Xs);
#ifdef UPDATE_EO
  free_spinor_field_f(eta);
  free_spinor_field_f(xi);
#endif
}


void force_hmc(double dt, suNg_av_field *force, void *vpar){
#ifdef UPDATE_EO
  spinor_field Xe, Xo, Ye, Yo;
#endif
#ifdef MEASURE_FORCEHMC
  double forcestat[2]={0.,0.}; /* used for computation of avr and max force */
#else
  double* forcestat=NULL;
#endif

  force_hmc_par *par = (force_hmc_par*)vpar;
  spinor_field *pf = par->pf;
  int n_iters = 0;
  (void) n_iters;
  init_force_hmc();
  set_dirac_mass(par->mass);
  set_twisted_mass(par->mu);


  /* check input types */
  _TWO_SPINORS_MATCHING(u_gauge,force);

  for (int k=0; k<par->n_pf; ++k) {

#ifndef UPDATE_EO
    double tmp;
    mshift_par mpar_ms;
    mpar_ms.err2 = par->inv_err2;
    mpar_ms.max_iter=0;
    mpar_ms.n = 1;
    mpar_ms.shift = &tmp; // mpar.n = 1
    mpar_ms.shift[0] = 0;


    if ( par->mu == 0 || par->hasenbusch != 0){
      /* X = H^{-1} pf[k] = D^{-1} g5 pf[k] */    
      spinor_field_zero_f(Xs);
      spinor_field_g5_assign_f(&pf[k]);
      n_iters += g5QMR_mshift(&mpar_ms, &D, &pf[k], Xs);
      spinor_field_g5_assign_f(&pf[k]);
      if(par->hasenbusch == 0) {
        /* Y  D^{-1} (  g5 X ) */
        spinor_field_g5_f(eta,Xs);
      } 
      else if (par->hasenbusch == 1) {
        /* Y = H^{-1} ( g5 pf[k] + b X ) = D^{-1} ( pf[k] + b g5 X ) */
        spinor_field_g5_f(eta,Xs);
        spinor_field_mul_f(eta,par->b,eta);
        spinor_field_add_assign_f(eta,&pf[k]);
      }
      else if (par->hasenbusch == 2) {//Twisted mass
        /* Y= -i D^{-1}(pf[k] + imu g5 X)*/
        double mu1 = par->mu;
        double mu2 = par->mu+par->b;
        double muS = mu2*mu2-mu1*mu1;
        lprintf("FORCE_HMC",50,"par->mu: %g, par->b: %g, muS: %g\n",static_mu,par->b,muS);
        spinor_field_g5_f(eta,Xs);
        spinor_field_mul_f(Xs,muS,Xs);
      }
      spinor_field_zero_f(Ys);
      n_iters += g5QMR_mshift(&mpar_ms, &D, eta, Ys);
    }
    else{
      n_iters += cg_mshift(&mpar_ms,QpQm_tm_alt,&pf[k],Xs);
      Qtm_p_alt(Ys,Xs);
    }
    
    if (par->hasenbusch != 1){
      force_fermion_core(Xs,Ys,force,dt,forcestat,0);
    }
    else{
      force_fermion_core(Xs,Ys,force,dt*par->b,forcestat,0);
    }

#else

    double tmp;
    mshift_par mpar;
    mpar.err2 = par->inv_err2;
    mpar.max_iter=0;
    mpar.n = 1;
    mpar.shift = &tmp; // mpar.n = 1
    mpar.shift[0] = 0;

    /* xi termporary on odd and eta temporary on even */
    if (par->mu==0){ 
      /* X_e = H^{-1} pf[k] */
      /* X_o = D_{oe} X_e = D_{oe} H^{-1} pf[k] */
      Xe=*Xs; Xe.type=&glat_even;
      Xo=*Xs; Xo.type=&glat_odd; Xo.ptr+=glat_odd.master_shift;
      spinor_field_zero_f(&Xe);
      /* H^{-1} pf = D^{-1} g5 pf */
      spinor_field_g5_assign_f(&pf[k]);
      
      mre_guess(&par->mpar,0, &Xe, &D, &pf[k]);
      n_iters+=g5QMR_mshift(&mpar, &D, &pf[k], &Xe);
      mre_store(&par->mpar,0, &Xe);
      spinor_field_g5_assign_f(&pf[k]);
      Dphi_(&Xo,&Xe);
      
      /* Y_e = H^{-1} ( g5 pf[k] + b X_e ) */
      /* Y_o = D_oe H^{-1} ( g5 pf[k] + b X_e ) */
      Ye=*Ys; Ye.type=&glat_even;
      Yo=*Ys; Yo.type=&glat_odd; Yo.ptr+=glat_odd.master_shift;      
      if(par->hasenbusch != 1) {
        spinor_field_copy_f(eta,&Xe);
      } else {
        spinor_field_g5_f(eta,&pf[k]);
        spinor_field_mul_add_assign_f(eta,par->b,&Xe);
      }
      
      spinor_field_zero_f(&Ye);
      spinor_field_g5_assign_f(eta);
      mre_guess(&par->mpar,1, &Ye, &D, eta);
      n_iters+=g5QMR_mshift(&mpar, &D, eta, &Ye);
      mre_store(&par->mpar,1, &Ye);
      spinor_field_g5_assign_f(eta);
      Dphi_(&Yo,&Ye);

      if (par->hasenbusch == 2){
        double muS = par->b*par->b;
        spinor_field_mul_f(Xs,muS,Xs);
      }
    }
    else{
      Ye=*Ys; Ye.type=&glat_even;
      Yo=*Ys; Yo.type=&glat_odd; Yo.ptr+=glat_odd.master_shift;
      spinor_field_zero_f(&Ye);
      /* Ye = 1/(QpQm+mu^2) \phi */
      mre_guess(&par->mpar,0, &Ye, QpQm_tm_alt, &pf[k]);
      n_iters += 2*cg_mshift(&mpar,QpQm_tm_alt, &pf[k],&Ye);
       mre_store(&par->mpar,0, &Ye);
      /* Yo = D_oe Y_e */
      Dphi_(&Yo,&Ye);

      Xe=*Xs; Xe.type=&glat_even;
      Xo=*Xs; Xo.type=&glat_odd; Xo.ptr+=glat_odd.master_shift;
      
      /* Xe = Qtm_m Ye */
      Qtm_m_alt(&Xe,&Ye);
      /* X_o = D_{oe} X_e = D_{oe} H^{-1} pf[k] */
      Dphi_(&Xo,&Xe);

      if (par->hasenbusch == 2){
        double mu1 = static_mu;
        double mu2 = static_mu+par->b;
        double muS = mu2*mu2-mu1*mu1;
        spinor_field_mul_f(Xs,muS,Xs);
      }

    }

    if (par->hasenbusch != 1){
      force_fermion_core(Xs,Ys,force,-dt,forcestat,0);
    }
    else{
      force_fermion_core(Xs,Ys,force,-dt*par->b,forcestat,0);
    }
#endif
  }
#ifdef MEASURE_FORCEHMC
  global_sum(forcestat,1);
  lprintf("FORCE_HMC",20,"HB %d: dt=%g, b=%g\n",par->hasenbusch,dt,par->b);
  if (par->hasenbusch !=1){
    forcestat[0]*=fabs(dt)*(_REPR_NORM2/_FUND_NORM2)/((double)(4*GLB_T*GLB_X*GLB_Y*GLB_Z));
    global_max(forcestat+1,1);
    forcestat[1]*=fabs(dt)*(_REPR_NORM2/_FUND_NORM2);
  }
  else{
    forcestat[0]*=fabs(dt*par->b)*(_REPR_NORM2/_FUND_NORM2)/((double)(4*GLB_T*GLB_X*GLB_Y*GLB_Z));
    global_max(forcestat+1,1);
    forcestat[1]*=fabs(dt*par->b)*(_REPR_NORM2/_FUND_NORM2);
  }
  force_ave[par->id]+=forcestat[0];
  force_max[par->id]+=forcestat[1];
  lprintf("FORCE_HMC",20,"avr dt |force| = %1.8e dt maxforce = %1.8e, dt = %1.8e \n",forcestat[0],forcestat[1],dt);
  n_inv_iter[par->id-1]+=n_iters;
#endif
}

