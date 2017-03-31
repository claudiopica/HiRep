/***************************************************************************\
 * Copyright (c) 2016, Jarno Rantaharju                                     *
 * All rights reserved.                                                     *
 \***************************************************************************/

//Monomial defining a fermion field that interacts with the gauge, if it exists,
//and through a four fermions interaction. The monomial four_fermion needs to be
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
static int mon_init = 1;

void hasen_ff_gaussian_pf(const struct _monomial *m)
{
	mon_hasenbusch_par *par = (mon_hasenbusch_par*)(m->data.par);
	gaussian_spinor_field(par->pf);
}

//Note: the g5 is here just to change the
//force calculation to resemple the standard one.
//It has no effect on the action
//
/* S = | (a D +b) D^{-1} g5 psi |^2 */
/* (a D +b) D^{-1} g5 psi = A */
/* psi = g5 D (a D + b)^{-1} A */
void hasen_ff_correct_pf(const struct _monomial *m)
{
	mon_hasenbusch_par *par = (mon_hasenbusch_par*)(m->data.par);
   double shift;

	mshift_par cgpar;
	shift = 0.;
	cgpar.n=1;
	cgpar.shift = &shift;
	cgpar.max_iter=0;
	cgpar.err2 = m->data.MT_prec;

	set_ff_dirac_mass(par->mass);
	set_ff_dirac_shift(par->dm);
	spinor_field_zero_f(tmp_pf);
	cg_mshift( &cgpar, &Dff_sq, par->pf, tmp_pf );
	Dff(par->pf, tmp_pf);
	spinor_field_g5_f(tmp_pf,par->pf);

	set_ff_dirac_shift(0.);
	spinor_field_g5_assign_f(tmp_pf);
	Dff_dagger(par->pf,tmp_pf);
}

void hasen_ff_correct_la_pf(const struct _monomial *m)
{
	mon_hasenbusch_par *par = (mon_hasenbusch_par*)(m->data.par);
	double shift;

	mshift_par cgpar;
	shift = 0.;
	cgpar.n=1;
	cgpar.shift = &shift;
	cgpar.max_iter=0;
	cgpar.err2 = m->data.MT_prec;

	spinor_field_zero_f(tmp_pf);
	set_ff_dirac_mass(par->mass);
	cg_mshift( &cgpar, &Dff_sq, par->pf, tmp_pf );
	Dff(par->pf,tmp_pf);
	spinor_field_g5_f(tmp_pf,par->pf);

	/* S = | (D+b) D^{-1} g5 psi |^2 */
	set_ff_dirac_shift(par->dm);
	spinor_field_g5_assign_f(tmp_pf);
	Dff_dagger(par->pf,tmp_pf);

	set_ff_dirac_shift(0.);
}

const spinor_field* hasen_ff_pseudofermion(const struct _monomial *m)
{
	mon_hasenbusch_par *par = (mon_hasenbusch_par*)(m->data.par);
	return par->pf;
}

void hasen_ff_add_local_action(const struct _monomial *m, scalar_field *loc_action)
{
	mon_hasenbusch_par *par = (mon_hasenbusch_par*)(m->data.par);
	pf_local_action(loc_action, par->pf);
}

void hasen_ff_free(struct _monomial *m)
{
	mon_hasenbusch_par *par = (mon_hasenbusch_par*)m->data.par;

	if(par->pf != NULL)
	{
		free_spinor_field_f(par->pf);
	}

	free(par);
	free(m);
}



struct _monomial* hasen_ff_create(const monomial_data *data)
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
	par->fpar.b = par->dm;
	par->fpar.hasenbusch = 2;
	par->fpar.mu = 0;
	par->fpar.momenta = &suN_momenta;

	// Setup chronological inverter
	mre_init(&(par->fpar.mpar), par->mre_past, data->force_prec);

	// Setup pointers to update functions
	m->free = &hasen_ff_free;
	m->update_force = &force_hmc_ff;
	m->force_par = &par->fpar;
	m->update_field = 0;
	m->field_par = 0;

	m->pseudofermion = &hasen_ff_pseudofermion;
	m->gaussian_pf = &hasen_ff_gaussian_pf;
	m->correct_pf = &hasen_ff_correct_pf;
	m->correct_la_pf = &hasen_ff_correct_la_pf;
	m->add_local_action = &hasen_ff_add_local_action;

	return m;
}
