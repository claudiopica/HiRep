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
#include "observables.h"
#include <stdlib.h>

void pg_init_traj(const struct _monomial *m) {
  /* empty */
}

void pg_gaussian_pf(const struct _monomial *m) {
  /* empty */
}

void pg_correct_pf(const struct _monomial *m) {
  /* empty */
}

void pg_correct_la_pf(const struct _monomial *m) {
  /* empty */
}

const spinor_field* pg_pseudofermion(const struct _monomial *m) {
  return NULL;
}

void pg_add_local_action(const struct _monomial *m, scalar_field *loc_action) {
  mon_pg_par *par = (mon_pg_par*)(m->data.par);
  
  /* Gauge action */
  _MASTER_FOR(&glattice,i) {
    *_FIELD_AT(loc_action,i) += -(par->beta/((double)NG))*local_plaq(i);
  }
}

void pg_free(struct _monomial *m) {
  mon_pg_par *par = (mon_pg_par*)m->data.par;
  free(par);
  free(m);
}

struct _monomial* pg_create(const monomial_data *data) {
  monomial *m = malloc(sizeof(*m));
  mon_pg_par *par = (mon_pg_par*)(data->par);
  
  // Copy data structure
  m->data = *data;
  
  // Allocate memory for spinor field
  /* empty */
  
  // Setup force parameters
  /* empty */
  
  // Setup pointers to update functions
  m->free = &pg_free;
  
  m->force_f = &force0;
  m->force_par = &par->beta;

  m->pseudofermion = &pg_pseudofermion;
  m->init_traj = &pg_init_traj;
  m->gaussian_pf = &pg_gaussian_pf;
  m->correct_pf = &pg_correct_pf;
  m->correct_la_pf = &pg_correct_la_pf;
  m->add_local_action = &pg_add_local_action;
  
  return m;
}
