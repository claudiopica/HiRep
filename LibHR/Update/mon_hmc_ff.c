/***************************************************************************\
 * Copyright (c) 2015, Martin Hansen, Claudio Pica, Jarno Rantaharju      *
 * All rights reserved.                                                   *
 \***************************************************************************/

//Monomial defining a fermion field that interacts with the gauge, if it exists,
//and trough a four fermions interaction. The monomial four_fermion needs to be
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

void hmc_ff_gaussian_pf(const struct _monomial *m)
{
	mon_hmc_par *par = (mon_hmc_par*)(m->data.par);
	gaussian_spinor_field(par->pf);
}

void hmc_ff_correct_pf(const struct _monomial *m)
{
	mon_hmc_par *par = (mon_hmc_par*)(m->data.par);

	/* compute H2^{1/2}*pf = H*pf */
	spinor_field_g5_f(tmp_pf, par->pf);
	set_ff_dirac_mass(par->mass);
	Dff_dagger(par->pf, tmp_pf);
}

void hmc_ff_correct_la_pf(const struct _monomial *m)
{
	mon_hmc_par *par = (mon_hmc_par*)(m->data.par);
	double shift;

	mshift_par cgpar;
	shift = 0.;
	cgpar.n=1;
	cgpar.shift = &shift;
	cgpar.max_iter=0;
	cgpar.err2 = m->data.MT_prec;

	set_ff_dirac_mass(par->mass);
	spinor_field_zero_f(tmp_pf);
	cg_mshift( &cgpar, &Dff_sq, par->pf, tmp_pf );
	Dff(par->pf,tmp_pf);
	spinor_field_g5_assign_f(par->pf);
}

const spinor_field* hmc_ff_pseudofermion(const struct _monomial *m)
{
	mon_hmc_par *par = (mon_hmc_par*)(m->data.par);
	return par->pf;
}

void hmc_ff_add_local_action(const struct _monomial *m, scalar_field *loc_action)
{
	mon_hmc_par *par = (mon_hmc_par*)(m->data.par);
	pf_local_action(loc_action, par->pf);

#ifdef UPDATE_EO
	/* If EO preconditioning is used, there the odd diagonal part of the
	 * fermion determinant is not included in the pseudo-fermion action.
	 * Det(A_o) = -exp( Trlog(A_o^ A_o) ) */
	_MASTER_FOR(&glat_odd,ix)
	{
		double a=0.;
		double ts = *_FIELD_AT(ff_sigma,ix);
		double tp = *_FIELD_AT(ff_pi,ix);
		double rho = 4. + par->mass + ts;
		int Nd=4;  //Number of dirac d.o.f.

		//Trace over site -> N_d*N_c. N_f == 2.
		a=-Nd*NF*log(rho*rho + tp*tp);
		*_FIELD_AT(loc_action,ix)+=a;
	}
#endif
}

void hmc_ff_free(struct _monomial *m)
{
	mon_hmc_par *par = (mon_hmc_par*)m->data.par;

	if(par->pf != NULL)
	{
		free_spinor_field_f(par->pf);
	}

	free(par);
	free(m);
}

struct _monomial* hmc_ff_create(const monomial_data *data)
{
	monomial *m = malloc(sizeof(*m));
	mon_hmc_par *par = (mon_hmc_par*)data->par;
  
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
	par->fpar.b = 0;
	par->fpar.hasenbusch = 0;
	par->fpar.mu = 0;
	par->fpar.momenta = &suN_momenta;

	// Setup chronological inverter
	mre_init(&(par->fpar.mpar), par->mre_past, data->force_prec);

	// Setup pointers to update functions
	m->free = &hmc_ff_free;
	m->update_force = &force_hmc_ff;
	m->force_par = &par->fpar;
	m->update_field = 0;
	m->field_par = 0;

	m->pseudofermion = &hmc_ff_pseudofermion;
	m->gaussian_pf = &hmc_ff_gaussian_pf;
	m->correct_pf = &hmc_ff_correct_pf;
	m->correct_la_pf = &hmc_ff_correct_la_pf;
	m->add_local_action = &hmc_ff_add_local_action;

	return m;
}
