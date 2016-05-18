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
#include <string.h>

void write_ff_scalar_fields(char filename[]);
void read_ff_scalar_fields(char filename[]);



void ff_init_traj(const struct _monomial *m) {
  //Generate momentum field. Clearly different structure from gauge, change?
  gaussian_scalar_field( ff_sigma_mom );
  gaussian_scalar_field( ff_pi_mom );
}

void ff_gaussian_pf(const struct _monomial *m) {
  /* empty */
}

void ff_correct_pf(const struct _monomial *m) {
  /* empty */
}

void ff_correct_la_pf(const struct _monomial *m) {
  /* empty */
}

const spinor_field* ff_pseudofermion(const struct _monomial *m) {
  return NULL;
}


void ff_add_local_action(const struct _monomial *m, scalar_field *loc_action) {
  mon_ff_par *par = (mon_ff_par*)(m->data.par);
  /* Gauge action */
  _MASTER_FOR(&glattice,i) {

    double a=0.;
    double ts = *_FIELD_AT(ff_sigma,i);
    double tp = *_FIELD_AT(ff_pi,i);
    
    double g = par->gamma*par->gamma*4.0;
    a  = tp*tp / g;
    a += ts*ts / g;

    ts = *_FIELD_AT(ff_sigma_mom,i);
    tp = *_FIELD_AT(ff_pi_mom,i);
    a += 0.5*(tp*tp + ts*ts);
    *_FIELD_AT(loc_action,i)+=a;
    //if(i==500) printf("action %g %f \n",a,g);
  }
}

void ff_free(struct _monomial *m) {
  mon_ff_par *par = (mon_ff_par*)m->data.par;

  free(ff_sigma);
  free(ff_pi);
  free(ff_sigma_mom);
  free(ff_pi_mom);
  free(par);
  free(m);
}

struct _monomial* ff_create(const monomial_data *data) {
  monomial *m = malloc(sizeof(*m));
  mon_ff_par *par = (mon_ff_par*)(data->par);

  four_fermion_active = 1;
  
  // Copy data structure
  m->data = *data;
  
  // Allocate auxiliary field
  ff_sigma = alloc_sfield(1,&glattice);
  ff_pi = alloc_sfield(1,&glattice);
  ff_sigma_mom = alloc_sfield(1,&glattice);
  ff_pi_mom = alloc_sfield(1,&glattice);

  //Read or set auxiliary field configuration or set cold value
  if(strcmp(par->start_config, "cold") == 0){
    set_scalar_field( ff_sigma , par->start_value );
    set_scalar_field( ff_pi , 0 );
  }

  // Allocate memory for spinor field
  /* empty */
  
  // Setup force parameters
  (par->fpar).gamma = par->gamma;

  // Setup pointers to update functions
  m->free = &ff_free;
  
  m->force_f = &force_hmc_auxfields;
  m->force_par = &par->fpar;

  m->pseudofermion = &ff_pseudofermion;
  m->init_traj = &ff_init_traj;
  m->gaussian_pf = &ff_gaussian_pf;
  m->correct_pf = &ff_correct_pf;
  m->correct_la_pf = &ff_correct_la_pf;
  m->add_local_action = &ff_add_local_action;
 
  return m;
}

