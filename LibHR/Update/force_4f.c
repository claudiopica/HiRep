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

/*  Update the momentum of the auxiliary fields with 
 *  the auxfield part of the action.
 */

void force_hmc_auxfields(double dt, void *vpar){
  force_auxfield_par *par = (force_auxfield_par*)vpar;
  double g = par->gamma*par->gamma*4.0;
  _MASTER_FOR(ff_sigma->type,i) {
    double ts = *_FIELD_AT(ff_sigma,i);
    double dm=0; 
    
    dm = 2.*ts/g;
    *_FIELD_AT(ff_sigma_mom,i)-=dt*dm;
   // if(i==500) printf("force %g %f \n",dm,g);
  }
  _MASTER_FOR(ff_pi->type,i) {
    double ts = *_FIELD_AT(ff_pi,i);
    double dm=0;
    dm = 2.*ts/g;
    *_FIELD_AT(ff_pi_mom,i)-=dt*dm;
  }

  
}



/* From force_hmc_4f.c: */
/* Y_e = (D^ D)^{-1} pf[k] */
/* Y_o = A_o^{-1} D_oe (D^ D)^{-1} pf[k]  */
/* X_e = g5 (D^)^(-1) pf[k]  */
/* X_o = (A)_o^{-1} g5 (D^)^(-1) pf[k] */

/* dA/dsigma = 1 and dA/dpi = i g5 
 * Thus, the force for sigma is 2*ReTr(Y_e X_e^)
 * and for pi 2*ReTr(i g5 Y_e X_e^)
 */
void force_hmc_auxfields_fermion(double dt, void *vpar, scalar_field *sigma_mom, scalar_field *pi_mom,spinor_field *Xs, spinor_field *Ys, int hasenbusch ){
  force_hmc_par *par = (force_hmc_par*)vpar;
#ifdef UPDATE_EO
  spinor_field Xe, Xo, Yo;
  Yo=*Ys; Yo.type=&glat_odd; Yo.ptr+=glat_odd.master_shift;
  Xe=*Xs; Xe.type=&glat_even;
  Xo=*Xs; Xo.type=&glat_odd; Xo.ptr+=glat_odd.master_shift;
  _MASTER_FOR(Xe.type,i) {
     suNf_spinor x,y;
     y = *_FIELD_AT(Ys,i);
     x = *_FIELD_AT(Xs,i);
     double prod=0;
     _spinor_g5_assign_f(y);
     _spinor_prod_re_f(prod,x,y);
     *_FIELD_AT(sigma_mom,i) += 2.*dt*prod;

     _spinor_g5_assign_f(y);
     _spinor_prod_im_f(prod,x,y);
     *_FIELD_AT(pi_mom,i) -= 2.*dt*prod;
  }
   _MASTER_FOR(Xo.type,i) {
     suNf_spinor x,y;
     y = *_FIELD_AT(Ys,i);
     x = *_FIELD_AT(Xs,i);
     double prod=0;
     _spinor_g5_assign_f(y);
     _spinor_prod_re_f(prod,x,y);
     *_FIELD_AT(sigma_mom,i) += 2.*dt*prod;
     //if(i==500) printf("fforce %g\n",-2.*prod);

     _spinor_g5_assign_f(y);
     _spinor_prod_im_f(prod,x,y);
     *_FIELD_AT(pi_mom,i) -= 2.*dt*prod;
  }
#else
  // Not implemented
#endif

#ifdef UPDATE_EO
 if(hasenbusch!=2) {
   double mass = par->mass;
  /* If EO preconditioning is used, the odd diagonal part of the 
   * fermion determinant is not included in the pseudo-fermion action.
   * Det(A_o) = exp( -Trlog(A_o^ A_o) ) */
   _MASTER_FOR(&glat_odd,i) {
     double ts = *_FIELD_AT(ff_sigma,i);
     double tp = *_FIELD_AT(ff_pi,i);
     double rho = 4. + mass + ts;
     double Ao2 = 1./(rho*rho+tp*tp);

     int Nd=4;
     int Ng = NG;
#ifdef REPR_ADJOINT
     Ng=Ng*Ng-1;
#endif

     //The action is a=-Nd*Ng*log(rho*rho + tp*tp), and the derivatives
     // da/dsigma = -2*Nd*Ng*rho/(rho*rho + tp*tp) and
     // da/dpi = -2*Nd*Ng*pi/(rho*rho + tp*tp),
     double dms = -2.*Nd*Ng*rho*Ao2;
     double dmp = -2.*Nd*Ng*tp*Ao2;

     *_FIELD_AT(ff_sigma_mom,i)-=dt*dms;
     *_FIELD_AT(ff_pi_mom,i)-=dt*dmp;
     //if(i==4000) printf("force %g mass %f\n",dt*dms,mass);
   }
 }
#endif

}



//Update the auxiliary fields
void update_auxfields(double dt, void *vpar){
  _MASTER_FOR(ff_sigma->type,i) {
    double tm = *_FIELD_AT(ff_sigma_mom,i);
    *_FIELD_AT(ff_sigma,i)+=dt*tm;
  }
  _MASTER_FOR(ff_pi->type,i) {
    double tm = *_FIELD_AT(ff_pi_mom,i);
    *_FIELD_AT(ff_pi,i)+=dt*tm;
  }
}



