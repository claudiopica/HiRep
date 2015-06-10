/***************************************************************************\
 * Copyright (c) 2008, Vincent Drach and Ari Hietanen                        *   
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



static void Mee_inv(spinor_field *out, double mass, double mu, spinor_field *in){
  /* (Mee^+)^-1 = (4+m+imu g5)^-(1)
     = (4+m-imu g5)/((4+m)^2+mu^2) */
  lprintf("FORCE_HMC_TM",50,"mass=%g, mu=%g\n",mass,mu);

  double norm = (4+mass)*(4+mass)+mu*mu;
  double rho = (4+mass)/norm;
  complex imu;
  imu.re=0;
  imu.im=-mu/norm;
  spinor_field_mul_f(out,rho,in);
  spinor_field_g5_mulc_add_assign_f(out,imu,in);
  //  spinor_field_mul_f(out,(4+mass),out);
}

static spinor_field *Xs=NULL, *Ys=NULL,*xi=NULL,*eta=NULL;

static void free_force_hmc_tm();

static void init_force_hmc_tm() {
  static int init=0;
  if (!init){
    Xs = alloc_spinor_field_f(2,&glattice);
    Ys = Xs+1;
    eta = alloc_spinor_field_f(1,&glat_even);
    xi = alloc_spinor_field_f(1,&glat_odd);
    atexit(&free_force_hmc_tm);
    init = 1;
  }
}

static void free_force_hmc_tm() {
  free_spinor_field_f(Xs);
  free_spinor_field_f(eta);
  free_spinor_field_f(xi);
}


void force_hmc_tm(double dt, suNg_av_field *force, void *vpar){
#ifndef UPDATE_EO
  error(1,1,"FORCE_HMC_TM","Use only with even odd preconditioned case\n");
#endif
  double forcestat[2]={0.,0.}; /* used for computation of avr and max force */
#ifdef MEASURE_FORCEHMC
  //  double forcestat[2]={0.,0.}; /* used for computation of avr and max force */
#endif

  init_force_hmc_tm();
  double tmp =0;
  spinor_field Xe, Xo, Ye, Yo;

  force_hmc_par *par = (force_hmc_par*)vpar;
  spinor_field *pf = par->pf;
  int n_iters = 0;
  (void) n_iters;
  set_dirac_mass(par->mass);
  set_twisted_mass(par->mu);

  mshift_par mpar;
  mpar.err2 = par->inv_err2;
  mpar.max_iter=0;
  mpar.n = 1;
  mpar.shift = &tmp; // mpar.n = 1
  mpar.shift[0] = 0;
  /* check input types */
  _TWO_SPINORS_MATCHING(u_gauge,force);

  for (int k=0; k<par->n_pf; ++k) {
    if (par->hasenbusch==0){// 1/( Q^2+mu^2)
      /* Ye = (\hat{Q}_+ \hat{Q}_-)^(-1)\phi */
      Ye=*Ys; Ye.type=&glat_even;
      Yo=*Ys; Yo.type=&glat_odd; Yo.ptr+=glat_odd.master_shift;
      spinor_field_zero_f(&Ye);
      mre_guess(&par->mpar,0, &Ye, &QpQm_tm, &pf[k]);
      n_iters += 2*cg_mshift(&mpar,QpQm_tm,&pf[k],&Ye);
      mre_store(&par->mpar,0, &Ye);
      /* Xe = (\hat{Q}+)^-1\phi = \hat{Q}_- * Ye */
      Xe=*Xs; Xe.type=&glat_even;
      Xo=*Xs; Xo.type=&glat_odd; Xo.ptr+=glat_odd.master_shift;
      Qtm_m(&Xe,&Ye);
      
      /* Yo = (M_ee^-)^-1 * M_eo Ye */
      Dphi_(xi,&Ye);
      Mee_inv(&Yo,par->mass,-par->mu,xi);
      
      /* Xo = (M_ee^+)^-1 M_eo Xe */
      Dphi_(xi,&Xe);
      Mee_inv(&Xo,par->mass,par->mu,xi);
      force_fermion_core(Xs,Ys,force,-dt,forcestat,0);
    }
    else{//   (Q^2+mu^2)/(Q^2)
      Ye=*Ys; Ye.type=&glat_even;
      Yo=*Ys; Yo.type=&glat_odd; Yo.ptr+=glat_odd.master_shift;
      Xe=*Xs; Xe.type=&glat_even;
      Xo=*Xs; Xo.type=&glat_odd; Xo.ptr+=glat_odd.master_shift;
      //First contribution to force
      //Xe = (Q+Q-)^{-1} W+pf[k]
      set_twisted_mass(par->mu+par->b);
      Qtm_p(&Ye,&pf[k]);
      set_twisted_mass(par->mu);
      spinor_field_zero_f(&Xe);
      mre_guess(&par->mpar,0, &Xe, &QpQm_tm, &pf[k]);
      n_iters += 2*cg_mshift(&mpar,QpQm_tm,&Ye,&Xe);
      mre_store(&par->mpar,0, &Xe);
      // Ye = Q_- Xe
      Qtm_m(&Ye,&Xe);
      /* Yo = (M_ee^+)^-1 * M_eo Ye */
      Dphi_(xi,&Ye);
      Mee_inv(&Yo,par->mass,par->mu,xi);
      /* Xo = (M_ee^-)^-1 M_eo Xe */
      Dphi_(xi,&Xe);
      Mee_inv(&Xo,par->mass,-par->mu,xi);
      force_fermion_core(Xs,Ys,force,-dt,forcestat,1);
      //Second contribution to force
      // Ye = pf[k];
      spinor_field_copy_f(&Ye,&pf[k]);
      /* Yo = (M_ee^+)^-1 * M_eo pf[k] */
      Dphi_(xi,&pf[k]);
      Mee_inv(&Yo,par->mass,par->mu+par->b,xi);
      /* Xo = (M_ee^-)^-1 M_eo Xe */
      Dphi_(xi,&Xe);
      Mee_inv(&Xo,par->mass,-par->mu-par->b,xi);
      force_fermion_core(Xs,Ys,force,dt,forcestat,2);
    }
  }
#ifdef MEASURE_FORCEHMC
    global_sum(forcestat,1);
    forcestat[0]*=fabs(dt)*((double) _REPR_NORM2/_FUND_NORM2)/((double)(4*GLB_T*GLB_X*GLB_Y*GLB_Z));
    global_max(forcestat+1,1);
    forcestat[1]*=fabs(dt)*((double) _REPR_NORM2/_FUND_NORM2);
    force_ave[par->id]+=fabs(forcestat[0]);
    force_max[par->id]+=forcestat[1];
    lprintf("FORCE_HMC",20,"avr dt |force| = %1.8e dt maxforce = %1.8e, dt = %1.8e \n",forcestat[0],forcestat[1],dt);
    n_inv_iter[par->id-1]+=n_iters;
#endif
}

