/***************************************************************************\
 * Copyright (c) 2017                                                     *
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

void scalar_blank(const struct _monomial *m)
{
	/* empty */
}

const spinor_field* scalar_pseudofermion(const struct _monomial *m)
{
	return NULL;
}

void scalar_add_local_action(const struct _monomial *m, scalar_field *loc_action)
{
	mon_scalar_par *par = (mon_scalar_par*)m->data.par;
	double Msq = par->mass;
	Msq = Msq*Msq + 8.0;
	double lambda = par->lambda;

	_MASTER_FOR(&glattice,ix)
	{
		suNg_vector *Sx, *Sup, UtimesSup;
		double SUSup, Ssq;
		suNg *U;

		Sx = _FIELD_AT(u_scalar,ix);
		_vector_prod_re_g(Ssq, *Sx, *Sx);
		*_FIELD_AT(loc_action,ix) += Msq*Ssq + lambda*Ssq*Ssq;

		for(int mu = 0; mu < 4; mu++)
		{
			Sup = _FIELD_AT(u_scalar, iup(ix,mu));
			U = _4FIELD_AT(u_gauge, ix, mu);
			_suNg_multiply(UtimesSup, *U, *Sup);
			_vector_prod_re_g(SUSup, *Sx, UtimesSup);
			*_FIELD_AT(loc_action,ix) -= 2.0*SUSup;
		}
	}
}

void scalar_free(struct _monomial *m)
{
	mon_scalar_par *par = (mon_scalar_par*)m->data.par;
	free(par);
	free(m);
}

struct _monomial* scalar_create(const monomial_data *data)
{
	monomial *m = malloc(sizeof(*m));
	mon_scalar_par *par = (mon_scalar_par*)data->par;

	// Allocate global field
	if(u_scalar == NULL)
	{
		u_scalar = alloc_scalar_field(&glattice);
	}

	// Copy data structure
	m->data = *data;

	//Setup force parameters
	par->force_par.mass = par->mass;
	par->force_par.lambda = par->lambda;
	par->force_par.momenta = &scalar_momenta;
	par->force_par.g_momenta = &suN_momenta;
	par->field_par.field = &u_scalar;
	par->field_par.momenta = &scalar_momenta;

	// Setup pointers to update functions
	m->free = &scalar_free;
	m->update_force = &force_scalar;
	m->force_par = &par->force_par;
	m->update_field = &update_scalar_field;
	m->field_par = &par->field_par;	

	m->pseudofermion = &scalar_pseudofermion;
	m->gaussian_pf = &scalar_blank;
	m->correct_pf = &scalar_blank;
	m->correct_la_pf = &scalar_blank;
	m->add_local_action = &scalar_add_local_action;

	return m;
}

