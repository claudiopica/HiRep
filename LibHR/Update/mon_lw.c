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

void lw_gaussian_pf(const struct _monomial *m)
{
	/* empty */
}

void lw_correct_pf(const struct _monomial *m)
{
	/* empty */
}

void lw_correct_la_pf(const struct _monomial *m)
{
	/* empty */
}

const spinor_field* lw_pseudofermion(const struct _monomial *m)
{
	return NULL;
}

void lw_add_local_action(const struct _monomial *m, scalar_field *loc_action)
{
	mon_lw_par *par = (mon_lw_par*)(m->data.par);
	lw_local_action(loc_action, par->beta, par->c0, par->c1);
}

void lw_free(struct _monomial *m)
{
	mon_lw_par *par = (mon_lw_par*)m->data.par;
	free(par);
	free(m);
}

struct _monomial* lw_create(const monomial_data *data)
{
	monomial *m = malloc(sizeof(*m));
	mon_lw_par *par = (mon_lw_par*)(data->par);

	// Copy data structure
	m->data = *data;
	par->c1 = (1.-par->c0)/8.;
	gauge_field_active = 1;

	// Setup force parameters
	par->force_par.beta = par->beta;
	par->force_par.c0 = par->c0;
	par->force_par.c1 = par->c1;
	par->force_par.momenta = &suN_momenta;
	par->field_par.field = &u_gauge;
	par->field_par.momenta = &suN_momenta;

	// Setup pointers to update functions
	m->free = &lw_free;
	m->update_force = &lw_force;
	m->force_par = &par->force_par;
	m->update_field = &update_gauge_field;
	m->field_par = &par->field_par;

	m->pseudofermion = &lw_pseudofermion;
	m->gaussian_pf = &lw_gaussian_pf;
	m->correct_pf = &lw_correct_pf;
	m->correct_la_pf = &lw_correct_la_pf;
	m->add_local_action = &lw_add_local_action;

	return m;
}
