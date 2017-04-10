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
  monomial *m;
  struct _mon_list *next;
} mon_list;

static int nmon=0;
static mon_list *action=NULL;

static void free_mon(mon_list *mon) {
  if (mon!=NULL) {
    mon->m->free(mon->m);
    free(mon);
  }
}

static void free_mon_list(mon_list *action) {
  while (action!=NULL) {
    mon_list *next=action->next;
    free_mon(action);
    action=next;
  }
}


static void free_mem() {
  free_mon_list(action);
}

const monomial *add_mon(monomial_data *mon_dat) {
  if(nmon == 0) { atexit(free_mem);} /* is it needed */
  mon_list *new_mon = malloc(sizeof(*new_mon));
  new_mon->next = action;
  action = new_mon;
  nmon++;
  
  /* create new monomial */
  switch (mon_dat->type) {
    case PureGauge:
      new_mon->m = pg_create(mon_dat);
      break;
    case LuscherWeisz:
      new_mon->m = lw_create(mon_dat);
      break;
    case HMC:
      new_mon->m = hmc_create(mon_dat);
      break;
    case RHMC:
      new_mon->m = rhmc_create(mon_dat);
      break;
    case TM:
      new_mon->m = tm_create(mon_dat);
      break;
    case TM_alt:
      new_mon->m = tm_alt_create(mon_dat);
      break;
    case Hasenbusch:
      new_mon->m = hasen_create(mon_dat);
      break;
    case Hasenbusch_tm:
      new_mon->m = hasen_tm_create(mon_dat);
      break;
    case Hasenbusch_tm_alt:
      new_mon->m = hasen_tm_alt_create(mon_dat);
      break;
    case FourFermion:
      new_mon->m = ff_create(mon_dat);
      break;
    case HMC_ff:
      new_mon->m = hmc_ff_create(mon_dat);
      break;
    case Hasenbusch_ff:
      new_mon->m = hasen_ff_create(mon_dat);
      break;
    case Scalar:
      new_mon->m = scalar_create(mon_dat);
      break;
    default:
      lprintf("MONOMIAL",0,"WARNING: unknown type!\n");
      break;
  }

  return new_mon->m;
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
  return curr->m;
}






