/***************************************************************************\
 * Copyright (c) 2008, Vincent Drach and Ari Hietanen                        *
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
#include "communications.h"

static spinor_field *Xs = NULL;
static spinor_field *Ys = NULL;
static spinor_field *xi = NULL;
static spinor_field *eta = NULL;

static void Mee_inv(spinor_field *out, double mass, double mu, spinor_field *in)
{
	/* (Mee^+)^-1 = (4+m+imu g5)^-(1) = (4+m-imu g5)/((4+m)^2+mu^2) */
	lprintf("FORCE_HMC_TM", 50, "mass=%g, mu=%g\n", mass, mu);
	double norm = (4+mass)*(4+mass)+mu*mu;
	double rho = (4+mass)/norm;
	complex imu = {0, -mu/norm};
	spinor_field_mul_f(out, rho, in);
	spinor_field_g5_mulc_add_assign_f(out, imu, in);
}

static void free_force_hmc_tm()
{
	free_spinor_field_f(Xs);
	free_spinor_field_f(eta);
	free_spinor_field_f(xi);
}

static void init_force_hmc_tm()
{
	static int init = 0;
	if(init == 0)
	{
		Xs = alloc_spinor_field_f(2, &glattice);
		Ys = Xs+1;
		Xs->type = &glat_even;
		Ys->type = &glat_even;
		eta = alloc_spinor_field_f(1, &glat_even);
		xi = alloc_spinor_field_f(1, &glat_odd);
		atexit(&free_force_hmc_tm);
		init = 1;
	}
}

void force_hmc_tm(double dt, void *vpar)
{
#ifndef UPDATE_EO
	error(1, 1, "FORCE_HMC_TM", "Use only with even odd preconditioned case\n");
#endif
	
	init_force_hmc_tm();
	fermion_force_begin();
	spinor_field Xo, Yo;
	
	Yo = *Ys;
	Yo.type = &glat_odd;
	Yo.ptr += glat_odd.master_shift;
	
	Xo = *Xs;
	Xo.type = &glat_odd;
	Xo.ptr += glat_odd.master_shift;
	
	int n_iters = 0;
	force_hmc_par *par = (force_hmc_par*)vpar;
	suNg_av_field *force = *par->momenta;
	spinor_field *pf = par->pf;
	set_dirac_mass(par->mass);
	set_twisted_mass(par->mu);

	double tmp;
	mshift_par mpar;
	mpar.err2 = par->inv_err2;
	mpar.max_iter = 0;
	mpar.n = 1;
	mpar.shift = &tmp;
	mpar.shift[0] = 0;
	
	if(par->hasenbusch == 0) // 1/(Q^2+mu^2)
	{
		/* Ye = (\hat{Q}_+ \hat{Q}_-)^(-1)\phi */
		mre_guess(&par->mpar, 0, Ys, &QpQm_tm, pf);
		n_iters += 2*cg_mshift(&mpar, QpQm_tm, pf, Ys);
		mre_store(&par->mpar, 0, Ys);

		/* Xe = (\hat{Q}+)^-1\phi = \hat{Q}_- * Ye */
		Qtm_m(Xs, Ys);

		/* Yo = (M_ee^-)^-1 * M_eo Ye */
		Dphi_(xi, Ys);
		Mee_inv(&Yo, par->mass, -par->mu, xi);

		/* Xo = (M_ee^+)^-1 M_eo Xe */
		Dphi_(xi, Xs);
		Mee_inv(&Xo, par->mass, par->mu, xi);
		force_fermion_core(Xs, Ys, 0, dt, 1.);
	}
	else // (Q^2+mu^2)/(Q^2)
	{
		//First contribution to force
		//Xe = (Q+Q-)^{-1} W + pf
		set_twisted_mass(par->mu+par->b);
		Qtm_p(Ys,pf);
		set_twisted_mass(par->mu);
		spinor_field_zero_f(Xs);
		mre_guess(&par->mpar, 0, Xs, &QpQm_tm, pf);
		n_iters += 2*cg_mshift(&mpar, QpQm_tm, Ys, Xs);
		mre_store(&par->mpar, 0, Xs);

		// Ye = Q_- Xe
		Qtm_m(Ys, Xs);

		/* Yo = (M_ee^+)^-1 * M_eo Ye */
		Dphi_(xi, Ys);
		Mee_inv(&Yo, par->mass, par->mu, xi);

		/* Xo = (M_ee^-)^-1 M_eo Xe */
		Dphi_(xi,Xs);
		Mee_inv(&Xo, par->mass, -par->mu, xi);
		force_fermion_core(Xs, Ys, 0, dt, 1.);

		//Second contribution to force
		// Ye = pf
		spinor_field_copy_f(Ys, pf);

		/* Yo = (M_ee^+)^-1 * M_eo pf */
		Dphi_(xi, pf);
		Mee_inv(&Yo, par->mass, par->mu+par->b, xi);

		/* Xo = (M_ee^-)^-1 M_eo Xe */
		Dphi_(xi, Xs);
		Mee_inv(&Xo, par->mass, -par->mu-par->b, xi);
		force_fermion_core(Xs, Ys, 0, -dt, 1.);
	}

	fermion_force_end(dt, force);
}

