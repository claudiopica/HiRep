/***************************************************************************\
 * Copyright (c) 2015, Martin Hansen, Claudio Pica                        *
 * All rights reserved.                                                   *
 \***************************************************************************/

#include "global.h"
#include "update.h"
#include "logger.h"
#include "memory.h"
#include "dirac.h"
#include "linear_algebra.h"
#include "inverters.h"
#include <stdlib.h>

static spinor_field *tmp_pf = NULL;

void hasen_tm_alt_init_traj(const struct _monomial *m) {
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

void hasen_tm_alt_gaussian_pf(const struct _monomial *m) {
   mon_hasenbusch_tm_par *par = (mon_hasenbusch_tm_par*)(m->data.par);
   gaussian_spinor_field(par->pf);
}

void hasen_tm_alt_correct_pf(const struct _monomial *m) {
   mon_hasenbusch_tm_par *par = (mon_hasenbusch_tm_par*)(m->data.par);
   double shift;
   mshift_par mpar;

   mpar.err2 = m->data.MT_prec;
   mpar.max_iter = 0;
   mpar.n = 1;
   mpar.shift = &shift;
   mpar.shift[0] = 0;

   /* S = | (g5 D + i mu ) D^{-1} g5 psi |^2 */
   /* ( g5 D + i mu ) D^{-1} g5 psi = A */
   /* psi = g5 D ( g5 D + i mu )^{-1}  A */
   set_dirac_mass(par->mass);
   set_twisted_mass(par->mu + par->dmu);
   spinor_field_zero_f(tmp_pf);
   tm_invert_alt(tmp_pf, par->pf, &mpar);
   set_twisted_mass(par->mu);
   set_twisted_mass(par->mu);
   Qtm_p_alt(par->pf, tmp_pf);
}

void hasen_tm_alt_correct_la_pf(const struct _monomial *m) {
   mon_hasenbusch_tm_par *par = (mon_hasenbusch_tm_par*)(m->data.par);
   double shift;
   mshift_par mpar;

   mpar.err2 = m->data.MT_prec;
   mpar.max_iter = 0;
   mpar.n = 1;
   mpar.shift = &shift;
   mpar.shift[0] = 0;

   /* S = | (g5 D + i mu ) D^{-1} g5 psi |^2 */ 
   spinor_field_copy_f(tmp_pf, par->pf);
   spinor_field_zero_f(par->pf);
   set_dirac_mass(par->mass);
   set_twisted_mass(par->mu);
   tm_invert_alt(par->pf, tmp_pf, &mpar);
   set_twisted_mass(par->mu + par->dmu);
   Qtm_p_alt(tmp_pf, par->pf);
   spinor_field_copy_f(par->pf, tmp_pf);
}

const spinor_field* hasen_tm_alt_pseudofermion(const struct _monomial *m) {
   mon_hasenbusch_tm_par *par = (mon_hasenbusch_tm_par*)(m->data.par);
   return par->pf;
}

void hasen_tm_alt_add_local_action(const struct _monomial *m, scalar_field *loc_action) {
   mon_hasenbusch_tm_par *par = (mon_hasenbusch_tm_par*)(m->data.par);
   /* pseudo fermion action = phi^2 */
   pf_local_action(loc_action, par->pf);
}

void hasen_tm_alt_free(struct _monomial *m) {
  mon_hasenbusch_tm_par *par = (mon_hasenbusch_tm_par*)m->data.par;
  if (par->pf!=NULL) {  free_spinor_field_f(par->pf); }
  free(par);
  free(m);
  /* It does NOT deallocate temporary spinor as it is shared by all mon */
}

struct _monomial* hasen_tm_alt_create(const monomial_data *data) {
  monomial *m = malloc(sizeof(*m));
  mon_hasenbusch_tm_par *par = (mon_hasenbusch_tm_par*)data->par;
  
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
  par->fpar.mu = par->mu;
  par->fpar.b = par->dmu;
  par->fpar.hasenbusch = 2;
  
  // Setup chronological inverter
  mre_init(&(par->fpar.mpar), par->mre_past, data->force_prec);
  
  // Setup pointers to update functions
  m->free = &hasen_tm_alt_free;
  
  m->force_f = &force_hmc;
  m->force_par = &par->fpar;

  m->pseudofermion = &hasen_tm_alt_pseudofermion;
  m->init_traj = &hasen_tm_alt_init_traj;
  m->gaussian_pf = &hasen_tm_alt_gaussian_pf;
  m->correct_pf = &hasen_tm_alt_correct_pf;
  m->correct_la_pf = &hasen_tm_alt_correct_la_pf;
  m->add_local_action = &hasen_tm_alt_add_local_action;
  
  return m;
}
