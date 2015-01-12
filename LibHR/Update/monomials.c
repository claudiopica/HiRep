/***************************************************************************\
 * Copyright (c) 2014, Martin Hansen, Claudio Pica                        *
 * All rights reserved.                                                      *
 \***************************************************************************/

#include "global.h"
#include "update.h"
#include "logger.h"
#include "memory.h"
#include "linear_algebra.h"
#include "rational_functions.h"
#include <stdlib.h>

typedef struct _mon_list {
  monomial m;
  struct _mon_list *next;
} mon_list;

static int nmon=0;
static mon_list *action=NULL;

void free_mon(mon_list *mon) {
  if (mon!=NULL) {
    if (mon->m.par!=NULL) {
      switch (mon->m.type) {
        case HMC:
          free_spinor_field_f(((mon_hmc_par*)mon->m.par)->pf);
          break;
        case RHMC:
          free_spinor_field_f(((mon_rhmc_par*)mon->m.par)->pf);
          break;
        case Hasenbusch:
          free_spinor_field_f(((mon_hasenbusch_par*)mon->m.par)->pf);
          break;
        default:
          lprintf("MONOMIAL",0,"WARNING: unknown type!\n");
          break;
      }
      free(mon->m.par);
    }
    free(mon);
  }
}

void free_mon_list(mon_list *action) {
  while (action!=NULL) {
    mon_list *next=action->next;
    free_mon(action);
    action=next;
  }
}

static spinor_field* alloc_1_sf(){
#ifdef UPDATE_EO
	return alloc_spinor_field_f(1, &glat_even);
#else
	return alloc_spinor_field_f(1, &glattice);
#endif
}

mon_list *alloc_mon(mon_type type) {
  mon_list *new_mon=malloc(sizeof(mon_list));
  switch (type) {
    case PureGauge:
      new_mon->m.par = malloc(sizeof(mon_pg_par));
      break;
    case HMC:
      new_mon->m.par = malloc(sizeof(mon_hmc_par));
      ((mon_hmc_par*)new_mon->m.par)->pf=alloc_1_sf();
      break;
    case RHMC:
      new_mon->m.par = malloc(sizeof(mon_rhmc_par));
      ((mon_rhmc_par*)new_mon->m.par)->pf=alloc_1_sf();
      break;
    case Hasenbusch:
      new_mon->m.par = malloc(sizeof(mon_hasenbusch_par));
      ((mon_hasenbusch_par*)new_mon->m.par)->pf=alloc_1_sf();
      break;
    default:
      lprintf("MONOMIAL",0,"WARNING: unknown type!\n");
      break;
  }
  
  return new_mon;
}

void mon_copy(monomial *new, monomial *old) {
  new->id=old->id;
  new->type=old->type;
  new->MT_prec=old->MT_prec;
  new->MD_prec=old->MD_prec;
  new->force_prec=old->force_prec;
  switch (new->type) {
    case PureGauge:
	  {
		  mon_pg_par *par_old = (mon_pg_par*)old->par;
		  mon_pg_par *par_new = (mon_pg_par*)new->par;
		  par_new->beta = par_old->beta;
	  }
      break;
    case HMC:
	  {
		  mon_hmc_par *par_old = (mon_hmc_par*)old->par;
		  mon_hmc_par *par_new = (mon_hmc_par*)new->par;
		  par_new->mass = par_old->mass;
		  par_new->fpar = par_old->fpar;
		  par_new->fpar.pf = par_new->pf;
		  spinor_field_copy_f(par_new->pf, par_old->pf);
	  }
      break;
    case RHMC:
	  {
		  mon_rhmc_par *par_old = (mon_rhmc_par*)old->par;
		  mon_rhmc_par *par_new = (mon_rhmc_par*)new->par;
		  par_new->mass = par_old->mass;
		  par_new->ratio = par_old->ratio;
		  r_app_alloc(&par_new->ratio);
		  par_new->fpar = par_old->fpar;
		  par_new->fpar.pf = par_new->pf;
		  par_new->fpar.ratio = &par_new->ratio;
		  spinor_field_copy_f(par_new->pf, par_old->pf);
	  }
      break;
    case Hasenbusch:
	  {
		  mon_hasenbusch_par *par_old = (mon_hasenbusch_par*)old->par;
		  mon_hasenbusch_par *par_new = (mon_hasenbusch_par*)new->par;
		  par_new->mass = par_old->mass;
		  par_new->dm = par_old->dm;
		  par_new->fpar = par_old->fpar;
  		  par_new->fpar.pf = par_new->pf;
		  spinor_field_copy_f(par_new->pf, par_old->pf);
	  }
      break;
    default:
      lprintf("MONOMIAL",0,"WARNING: unknown type!\n");
      break;
  }
}

static void free_mem() {
  free_mon(action);
}

const monomial *add_mon(monomial *mon) {
  if(nmon == 0) { atexit(free_mem);}
  mon_list *new_mon = alloc_mon(mon->type);
  mon_copy(&(new_mon->m), mon);
  new_mon->next = action;
  action = new_mon;
  nmon++;
  return &(new_mon->m);
}


int num_mon() {
  return nmon;
}

const monomial *mon_n(int i) {
  mon_list *curr = action;
#ifndef NDEBUG
  error((i<0 || i>=nmon),1,"mon_n","Wrong monomial number\n");
#endif
  while (i>0) { curr=curr->next; i--; }
  return &(curr->m);
}






