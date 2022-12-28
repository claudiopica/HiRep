/***************************************************************************\
 * Copyright (c) 2015, Martin Hansen, Claudio Pica                        *
 * All rights reserved.                                                   *
 \***************************************************************************/

#include "update.h"
#include "libhr_core.h"
#include "memory.h"
#include "random.h"
#include "Inverters/linear_algebra.h"

static spinor_field *tmp_pf = NULL;
static int mon_init = 1;

static void hasen_tm_gaussian_pf(const struct _monomial *m)
{
	mon_hasenbusch_tm_par *par = (mon_hasenbusch_tm_par*)(m->data.par);
	gaussian_spinor_field(par->pf);
}

static void hasen_tm_correct_pf(const struct _monomial *m)
{
	mon_hasenbusch_tm_par *par = (mon_hasenbusch_tm_par*)(m->data.par);
	double shift;

	mshift_par mpar;
	mpar.err2 = m->data.MT_prec;
	mpar.max_iter = 0;
	mpar.n = 1;
	mpar.shift = &shift;
	mpar.shift[0] = 0;

	/* Note the different convetion between hasenbuch_tm and hasenbusch_tm_alt*/
	/* S = | D^{-1} g5  (g5 D + i mu )  psi |^2 */
	/* D^{-1} g5 ( g5 D + i mu )  psi = A */
	/* psi = ( g5 D + i mu )^{-1}  g5 D  A */
	spinor_field_copy_f(tmp_pf, par->pf);
	set_dirac_mass(par->mass);
	set_twisted_mass(par->mu);
	Qtm_p(par->pf, tmp_pf);

	set_twisted_mass(par->mu + par->dmu);
	spinor_field_zero_f(tmp_pf);
	tm_invert(tmp_pf, par->pf, &mpar);
	spinor_field_copy_f(par->pf, tmp_pf);
}

static void hasen_tm_correct_la_pf(const struct _monomial *m)
{
	mon_hasenbusch_tm_par *par = (mon_hasenbusch_tm_par*)(m->data.par);
	double shift;

	mshift_par mpar;
	mpar.err2 = m->data.MT_prec;
	mpar.max_iter = 0;
	mpar.n = 1;
	mpar.shift = &shift;
	mpar.shift[0] = 0;

	/* S = |  D^{-1} g5 (g5 D + i mu ) psi |^2 */
	set_dirac_mass(par->mass);
	set_twisted_mass(par->mu + par->dmu);
	Qtm_p(tmp_pf, par->pf);
	spinor_field_zero_f(par->pf);
	set_twisted_mass(par->mu);
	tm_invert(par->pf, tmp_pf, &mpar);
}

static const spinor_field* hasen_tm_pseudofermion(const struct _monomial *m)
{
	mon_hasenbusch_tm_par *par = (mon_hasenbusch_tm_par*)(m->data.par);
	return par->pf;
}

static void hasen_tm_add_local_action(const struct _monomial *m, scalar_field *loc_action)
{
	mon_hasenbusch_tm_par *par = (mon_hasenbusch_tm_par*)(m->data.par);
	pf_local_action(loc_action, par->pf);
}

static void hasen_tm_free(struct _monomial *m)
{
	mon_hasenbusch_tm_par *par = (mon_hasenbusch_tm_par*)m->data.par;

	if(par->pf != NULL)
	{
		free_spinor_field_f(par->pf);
	}

	free(par);
	free(m);
}

struct _monomial* hasen_tm_create(const monomial_data *data)
{
	monomial *m = malloc(sizeof(*m));
	mon_hasenbusch_tm_par *par = (mon_hasenbusch_tm_par*)data->par;
  
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
	par->fpar.mu = par->mu;
	par->fpar.b = par->dmu;
	par->fpar.hasenbusch = 3;
	par->fpar.logdet = 0;
	par->fpar.momenta = &suN_momenta;

	// Setup chronological inverter
	mre_init(&(par->fpar.mpar), par->mre_past, data->force_prec);

	// Setup pointers to update functions
	m->free = &hasen_tm_free;
	m->update_force = &force_hmc_tm;
	m->force_par = &par->fpar;
	m->update_field = 0;
	m->field_par = 0;

	m->pseudofermion = &hasen_tm_pseudofermion;
	m->gaussian_pf = &hasen_tm_gaussian_pf;
	m->correct_pf = &hasen_tm_correct_pf;
	m->correct_la_pf = &hasen_tm_correct_la_pf;
	m->add_local_action = &hasen_tm_add_local_action;

	return m;
}
