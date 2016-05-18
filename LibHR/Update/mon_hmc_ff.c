/***************************************************************************\
 * Copyright (c) 2015, Martin Hansen, Claudio Pica, Jarno Rantaharju      *
 * All rights reserved.                                                   *
 \***************************************************************************/

//Monomial defining a fermion field that interacts with the gauge, if it exists,
//and trough a four fermions interaction. The monomial four_fermion needs to be
//included to define the auxiliary fields.

#include "global.h"
#include "update.h"
#include "logger.h"
#include "memory.h"
#include "dirac.h"
#include "linear_algebra.h"
#include "inverters.h"
#include <stdlib.h>
#include <math.h>

static spinor_field *tmp_pf = NULL;

void hmc_ff_init_traj(const struct _monomial *m) {
   static int init = 1;
   
   /* allocate temporary memory */
   if (init) {
#ifdef UPDATE_EO
      tmp_pf = alloc_spinor_field_f(1, &glat_even); /* even lattice for preconditioned dynamics */
#else
      tmp_pf = alloc_spinor_field_f(1, &glattice); /* global lattice */
#endif
      
      init = 0;
   }
}

void hmc_ff_gaussian_pf(const struct _monomial *m) {
   mon_hmc_par *par = (mon_hmc_par*)(m->data.par);
   gaussian_spinor_field(par->pf);
}

void hmc_ff_correct_pf(const struct _monomial *m) {
   mon_hmc_par *par = (mon_hmc_par*)(m->data.par);
   
   /* compute H2^{1/2}*pf = H*pf */
   spinor_field_g5_f(tmp_pf, par->pf);
   set_ff_dirac_mass(par->mass);
   Dff_dagger(par->pf, tmp_pf);
}

void hmc_ff_correct_la_pf(const struct _monomial *m) {
   mon_hmc_par *par = (mon_hmc_par*)(m->data.par);
   double shift;

   mshift_par cgpar;
   shift = 0.;
   cgpar.n=1;
   cgpar.shift = &shift;
   cgpar.max_iter=0;
   cgpar.err2 = m->data.MT_prec;
  
   set_ff_dirac_mass(par->mass);
   spinor_field_zero_f(tmp_pf);
   cg_mshift( &cgpar, &Dff_sq, par->pf, tmp_pf );
   Dff(par->pf,tmp_pf);
   spinor_field_g5_assign_f(par->pf);
}

const spinor_field* hmc_ff_pseudofermion(const struct _monomial *m) {
   mon_hmc_par *par = (mon_hmc_par*)(m->data.par);
   return par->pf;
}

void hmc_ff_add_local_action(const struct _monomial *m, scalar_field *loc_action) {
   mon_hmc_par *par = (mon_hmc_par*)(m->data.par);
   /* pseudo fermion action = phi^2 */
   pf_local_action(loc_action, par->pf);

#ifdef UPDATE_EO
  /* If EO preconditioning is used, there the odd diagonal part of the 
   * fermion determinant is not included in the pseudo-fermion action.
   * Det(A_o) = -exp( Trlog(A_o^ A_o) ) */
   _MASTER_FOR(&glat_odd,i) {
     double a=0.;
     double ts = *_FIELD_AT(ff_sigma,i);
     double tp = *_FIELD_AT(ff_pi,i);
     double rho = 4. + par->mass + ts;

     int Nd=4;  //Number of dirac d.o.f.

     //Trace over site -> N_d*N_c. N_f == 2.
     a=-Nd*NF*log(rho*rho + tp*tp);
     *_FIELD_AT(loc_action,i)+=a;
     //if(i==4000) printf(" TrlogA, odd: %g   %g  mass %g \n", a, *_FIELD_AT(loc_action,i), par->mass); 
   }
#endif
}

void hmc_ff_free(struct _monomial *m) {
  mon_hmc_par *par = (mon_hmc_par*)m->data.par;
  if (par->pf!=NULL) {  free_spinor_field_f(par->pf); }
  free(par);
  free(m);
  /* It does NOT deallocate temporary spinor as it is shared by all mon */
}




struct _monomial* hmc_ff_create(const monomial_data *data) {
  monomial *m = malloc(sizeof(*m));
  mon_hmc_par *par = (mon_hmc_par*)data->par;
  
  // Copy data structure
  m->data = *data;
  
  // Allocate memory for spinor field
#ifdef UPDATE_EO
  par->pf=alloc_spinor_field_f(1, &glat_even);
#else
  par->pf=alloc_spinor_field_f(1, &glattice);
#endif
  
  // Setup force parameters
  par->fpar.id=data->id;
  par->fpar.n_pf = 1;
  par->fpar.pf = par->pf;
  par->fpar.inv_err2 = data->force_prec;
  par->fpar.inv_err2_flt = 1e-6;
  par->fpar.mass = par->mass;
  par->fpar.b = 0;
  par->fpar.hasenbusch = 0;
  par->fpar.mu = 0;
  
  // Setup chronological inverter
  mre_init(&(par->fpar.mpar), par->mre_past, data->force_prec);
  
  // Setup pointers to update functions
  m->free = &hmc_ff_free;
  
  m->force_f = &force_hmc_ff;
  m->force_par = &par->fpar;

  m->pseudofermion = &hmc_ff_pseudofermion;
  m->init_traj = &hmc_ff_init_traj;
  m->gaussian_pf = &hmc_ff_gaussian_pf;
  m->correct_pf = &hmc_ff_correct_pf;
  m->correct_la_pf = &hmc_ff_correct_la_pf;
  m->add_local_action = &hmc_ff_add_local_action;
 
  return m;
}
