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
static int mon_init = 1;

void hasen_gaussian_pf(const struct _monomial *m)
{
	mon_hasenbusch_par *par = (mon_hasenbusch_par*)(m->data.par);
	gaussian_spinor_field(par->pf);
}

void hasen_correct_pf(const struct _monomial *m)
{
	mon_hasenbusch_par *par = (mon_hasenbusch_par*)(m->data.par);
	double shift;

	mshift_par mpar;
	mpar.err2 = m->data.MT_prec;
	mpar.max_iter = 0;
	mpar.n = 1;
	mpar.shift = &shift;
	mpar.shift[0] = 0;

	/* compute H(m)D^{-1}(m+dm)*pf */
	set_dirac_mass(par->mass + par->dm);
	spinor_field_zero_f(tmp_pf); /* mshift inverter uses this as initial guess for 1 shift */
	g5QMR_mshift(&mpar, &D, par->pf, tmp_pf);
	set_dirac_mass(par->mass);
	H(par->pf, tmp_pf);
}

void hasen_correct_la_pf(const struct _monomial *m)
{
	mon_hasenbusch_par *par = (mon_hasenbusch_par*)(m->data.par);
	double shift;

	mshift_par mpar;
	mpar.err2 = m->data.MT_prec;
	mpar.max_iter = 0;
	mpar.n = 1;
	mpar.shift = &shift;
	mpar.shift[0] = 0;

	/* compute D(m+dm)D^{-1}(m)*g5*pf */
	spinor_field_g5_assign_f(par->pf);
	set_dirac_mass(par->mass);
	spinor_field_zero_f(tmp_pf);
	g5QMR_mshift(&mpar, &D, par->pf, tmp_pf);
	set_dirac_mass(par->mass + par->dm);
	D(par->pf, tmp_pf);
}

const spinor_field* hasen_pseudofermion(const struct _monomial *m)
{
	mon_hasenbusch_par *par = (mon_hasenbusch_par*)(m->data.par);
	return par->pf;
}

void hasen_add_local_action(const struct _monomial *m, scalar_field *loc_action)
{
	mon_hasenbusch_par *par = (mon_hasenbusch_par*)(m->data.par);
	pf_local_action(loc_action, par->pf);
}

void hasen_free(struct _monomial *m)
{
	mon_hasenbusch_par *par = (mon_hasenbusch_par*)m->data.par;

	if(par->pf != NULL)
	{
		free_spinor_field_f(par->pf);
	}

	free(par);
	free(m);
}

struct _monomial* hasen_create(const monomial_data *data)
{
	monomial *m = malloc(sizeof(*m));
	mon_hasenbusch_par *par = (mon_hasenbusch_par*)data->par;
  
	// Copy data structure
	m->data = *data;
  
	// Allocate memory for spinor field
	if(mon_init)
	{
		tmp_pf = alloc_spinor_field_f(1, &glat_default);
		mon_init = 0;
	}
	par->pf = alloc_spinor_field_f(1, &glat_default);

	// Setup force parameters
	par->fpar.id = data->id;
	par->fpar.n_pf = 1;
	par->fpar.pf = par->pf;
	par->fpar.inv_err2 = data->force_prec;
	par->fpar.inv_err2_flt = 1e-6;
	par->fpar.mass = par->mass;
#ifdef UPDATE_EO
	par->fpar.b = (4+par->mass+par->dm)*(4+par->mass+par->dm)-(4+par->mass)*(4+par->mass);
#else
	par->fpar.b = par->dm;
#endif
	par->fpar.hasenbusch = 1;
	par->fpar.mu = 0;
	par->fpar.logdet = 0;
	par->fpar.momenta = &suN_momenta;

	// Setup chronological inverter
	mre_init(&(par->fpar.mpar), par->mre_past, data->force_prec);

	// Setup pointers to update functions
	m->free = &hasen_free;
	m->update_force = &force_hmc;
	m->force_par = &par->fpar;
	m->update_field = 0;
	m->field_par = 0;

	m->pseudofermion = &hasen_pseudofermion;
	m->gaussian_pf = &hasen_gaussian_pf;
	m->correct_pf = &hasen_correct_pf;
	m->correct_la_pf = &hasen_correct_la_pf;
	m->add_local_action = &hasen_add_local_action;

	return m;
}
