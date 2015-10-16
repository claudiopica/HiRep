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
#include "rational_functions.h"
#include <stdlib.h>

static spinor_field *tmp_pf = NULL;
static rational_app r_APP = {0};  /* used for computing HB and MT for RHMC monomials */

static int gcd(int a, int b) {
   while (b!=0){ int t=b; b=a%t; a=t; }
   return a;
}

static void reduce_fraction(int *a, int *b){
   int f=gcd(abs(*a),abs(*b));
   if (*b!=0 && f!=1){ *a/=f; *b/=f; }
}

void rhmc_init_traj(const struct _monomial *m) {
   mon_rhmc_par *par = (mon_rhmc_par*)(m->data.par);
   double minev, maxev; /* min and max eigenvalue of H^2 */
   static int init = 1;

   /* allocate temporary memory */
   if (init) {
#ifdef UPDATE_EO
      tmp_pf = alloc_spinor_field_f(1, &glat_even); /* even lattice for preconditioned dynamics */
#else
      tmp_pf = alloc_spinor_field_f(1, &glattice); /* global lattice */
#endif

      /* Pre-allocate rational approximation */
      r_APP.order = 16;
      r_app_alloc(&r_APP);

      init = 0;
   }

   /* find min/max and best rational approximation */
   set_dirac_mass(par->mass);
   find_spec_H2(&maxev, &minev); /* find spectral interval of H^2 */
   r_app_set(&(par->ratio), minev, maxev);
}

void rhmc_gaussian_pf(const struct _monomial *m) {
   mon_rhmc_par *par = (mon_rhmc_par*)(m->data.par);
   gaussian_spinor_field(par->pf);
}

void rhmc_correct_pf(const struct _monomial *m) {
   mon_rhmc_par *par = (mon_rhmc_par*)(m->data.par);

   /* r_APP = x^{-n/(2*d)} */
   r_APP.rel_error = m->data.MT_prec;
   /* use n=-n and d=2*d respect to the r_app used for the MD */
   r_APP.n = -par->ratio.n;
   r_APP.d = 2*par->ratio.d;
   reduce_fraction(&r_APP.n, &r_APP.d);
   r_app_set(&r_APP, par->ratio.min, par->ratio.max);

   set_dirac_mass(par->mass);
   rational_func(&r_APP, &H2, par->pf, par->pf);
}

void rhmc_correct_la_pf(const struct _monomial *m) {
   mon_rhmc_par *par = (mon_rhmc_par*)(m->data.par);

   /* r_APP = x^{n/(2*d)} */
   r_APP.rel_error = m->data.MT_prec;
   /* use n=n and d=2*d respect to the r_app used for the MD */
   r_APP.n = par->ratio.n;
   r_APP.d = 2*par->ratio.d;
   reduce_fraction(&r_APP.n, &r_APP.d);
   r_app_set(&r_APP, par->ratio.min, par->ratio.max);

   set_dirac_mass(par->mass);
   rational_func(&r_APP, &H2, par->pf, par->pf);
}

const spinor_field* rhmc_pseudofermion(const struct _monomial *m) {
   mon_rhmc_par *par = (mon_rhmc_par*)(m->data.par);
   return par->pf;
}

void rhmc_add_local_action(const struct _monomial *m, scalar_field *loc_action) {
   mon_rhmc_par *par = (mon_rhmc_par*)(m->data.par);
   /* pseudo fermion action = phi^2 */
   pf_local_action(loc_action, par->pf);
}


void rhmc_free(struct _monomial *m) {
  mon_rhmc_par *par = (mon_rhmc_par*)m->data.par;
  if (par->pf!=NULL) {  free_spinor_field_f(par->pf); }
  r_app_free(&par->ratio);
  free(par);
  free(m);
  /* It does NOT deallocate temporary spinor as it is shared by all mon */
}


struct _monomial* rhmc_create(const monomial_data *data) {
  monomial *m = malloc(sizeof(*m));
  mon_rhmc_par *par = (mon_rhmc_par*)data->par;

  // Copy data structure
  m->data = *data;

  // Allocate memory for mon_rhmc_par
#ifdef UPDATE_EO
  par->pf=alloc_spinor_field_f(1, &glat_even);
#else
  par->pf=alloc_spinor_field_f(1, &glattice);
#endif

  // Setup force parameters
  //  par->fpar.id=data->id;
  par->fpar.n_pf = 1;
  par->fpar.pf = par->pf;
  par->fpar.inv_err2 = data->force_prec;
  par->fpar.mass = par->mass;
  par->ratio.rel_error = data->MD_prec;
  par->fpar.ratio = &par->ratio;
  par->ratio.order = 16;
  r_app_alloc(&par->ratio);

  // Setup pointers to update functions
  m->free = &rhmc_free;

  m->force_f = &force_rhmc;
  m->force_par = &par->fpar;

  m->pseudofermion = &rhmc_pseudofermion;
  m->init_traj = &rhmc_init_traj;
  m->gaussian_pf = &rhmc_gaussian_pf;
  m->correct_pf = &rhmc_correct_pf;
  m->correct_la_pf = &rhmc_correct_la_pf;
  m->add_local_action = &rhmc_add_local_action;

  return m;
}
