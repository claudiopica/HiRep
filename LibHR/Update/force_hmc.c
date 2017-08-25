/***************************************************************************\
* Copyright (c) 2008, Claudio Pica, Martin Hansen                           *
* All rights reserved.                                                      * 
\***************************************************************************/

#include "global.h"
#include "update.h"
#include "suN.h"
#include "utils.h"
#include "dirac.h"
#include "inverters.h"
#include "rational_functions.h"
#include "representation.h"
#include "logger.h"
#include "linear_algebra.h"
#include "memory.h"
#include "clover_tools.h"
#include "communications.h"

spinor_field *Xs=NULL;
spinor_field *Ys=NULL;
spinor_field *eta=NULL;
#ifdef UPDATE_EO
static spinor_field *xi=NULL;
#endif 

void free_force_hmc()
{
	free_spinor_field_f(Xs);
#ifdef UPDATE_EO
	free_spinor_field_f(eta);
	free_spinor_field_f(xi);
#endif
}

void init_force_hmc()
{
	static int init = 0;
	if(init == 0)
	{
#ifndef UPDATE_EO
		Xs = alloc_spinor_field_f(3, &glattice);
		Ys = Xs+1;
		eta = Ys+1;
#else
		Xs = alloc_spinor_field_f(2, &glattice);
		Ys = Xs+1;
		eta = alloc_spinor_field_f(1, &glat_even);
		xi = alloc_spinor_field_f(1, &glat_odd);
		Xs->type = &glat_even;
		Ys->type = &glat_even;
#endif
		init = 1;
		atexit(free_force_hmc);
	}
}

void force_hmc(double dt, void *vpar)
{
	int n_iters = 0;
	force_hmc_par *par = (force_hmc_par*)vpar;
	suNg_av_field *force = *par->momenta;
	spinor_field *pf = par->pf;

	init_force_hmc();
	set_dirac_mass(par->mass);
	set_twisted_mass(par->mu);
	fermion_force_begin();

	double tmp;
	mshift_par mpar;
	mpar.err2 = par->inv_err2;
	mpar.max_iter = 0;
	mpar.n = 1;
	mpar.shift = &tmp;
	mpar.shift[0] = 0;

#ifndef UPDATE_EO

	if(par->mu == 0 || par->hasenbusch != 0)
	{
		/* X = H^{-1} pf = D^{-1} g5 pf */
		spinor_field_zero_f(Xs);
		spinor_field_g5_assign_f(pf);
		n_iters += g5QMR_mshift(&mpar, &D, pf, Xs);
		spinor_field_g5_assign_f(pf);

		if(par->hasenbusch == 0)
		{
			/* Y  D^{-1} ( g5 X ) */
			spinor_field_g5_f(eta, Xs);
		}
		else if(par->hasenbusch == 1)
		{
			/* Y = H^{-1} ( g5 pf[k] + b X ) = D^{-1}( pf + b g5 X ) */
			spinor_field_g5_f(eta, Xs);
			spinor_field_mul_f(eta, par->b, eta);
			spinor_field_add_assign_f(eta, pf);
		}
		else if(par->hasenbusch == 2)
		{
			/* Y= -i D^{-1} ( pf[k] + imu g5 X )*/
			double mu1 = par->mu;
			double mu2 = par->mu + par->b;
			double muS = mu2*mu2 - mu1*mu1;
			spinor_field_g5_f(eta, Xs);
			spinor_field_mul_f(Xs, muS, Xs);
		}

		spinor_field_zero_f(Ys);
		n_iters += g5QMR_mshift(&mpar, &D, eta, Ys);
	}
	else
	{
		n_iters += cg_mshift(&mpar, QpQm_tm_alt, pf, Xs);
		Qtm_p_alt(Ys,Xs);
	}

#else

	if(par->mu == 0)
	{
		/* X_e = H^{-1} pf */
		/* X_o = D_{oe} X_e = D_{oe} H^{-1} pf */
		spinor_field_g5_assign_f(pf);
		mre_guess(&par->mpar, 0, Xs, &D, pf);
		n_iters += g5QMR_mshift(&mpar, &D, pf, Xs);
		mre_store(&par->mpar, 0, Xs);
		spinor_field_g5_assign_f(pf);

		/* Y_e = H^{-1} ( g5 pf + b X_e ) */
		/* Y_o = D_oe H^{-1} ( g5 pf + b X_e ) */
		if(par->hasenbusch != 1)
		{
			spinor_field_copy_f(eta, Xs);
		}
		else
		{
			spinor_field_g5_f(eta, pf);
			spinor_field_mul_add_assign_f(eta, par->b, Xs);
		}

		spinor_field_g5_assign_f(eta);
		mre_guess(&par->mpar, 1, Ys, &D, eta);
		n_iters += g5QMR_mshift(&mpar, &D, eta, Ys);
		mre_store(&par->mpar, 1, Ys);
		spinor_field_g5_assign_f(eta);

		if(par->hasenbusch == 2)
		{
			double muS = par->b*par->b;
			spinor_field_mul_f(Xs, muS, Xs);
		}
	}
	else
	{
		/* Ye = 1/(QpQm+mu^2) \phi */
		mre_guess(&par->mpar, 0, Ys, QpQm_tm_alt, pf);
		n_iters += 2*cg_mshift(&mpar, QpQm_tm_alt, pf, Ys);
		mre_store(&par->mpar, 0, Ys);
		Qtm_m_alt(Xs, Ys);

		if(par->hasenbusch == 2)
		{
			double mu1 = par->mu;
			double mu2 = par->mu + par->b;
			double muS = mu2*mu2 - mu1*mu1;
			spinor_field_mul_f(Xs, muS, Xs);
		}
	}

#endif

	if(par->hasenbusch != 1)
	{
		force_fermion_core(Xs, Ys, 1, dt, 1.);
	}
	else
	{
		force_fermion_core(Xs, Ys, 1, dt, par->b);
	}

#ifdef WITH_CLOVER_EO

	if(par->logdet)
	{
		force_clover_logdet(par->mass, 2.); // 2 = # of flavors
	}

#endif

	fermion_force_end(dt, force);
}
