/***************************************************************************\
* Copyright (c) 2017
* All rights reserved.                                                      * 
\***************************************************************************/

#include "global.h"
#include "update.h"
#include "suN.h"
#include "utils.h"
#include "representation.h"
#include "logger.h"
#include "communications.h"
#include <stdio.h>
#include <math.h>

#ifndef WITH_QUATERNIONS

static void force_scalar_s(double dt, void *vpar)
{
	force_scalar_par *par = (force_scalar_par*)vpar;
	suNg_scalar_field *force = *par->momenta;

	_MASTER_FOR(&glattice,ix)
	{
		suNg_vector f;
		_vector_zero_g(f);

		for(int mu = 0; mu < 4; mu++)
		{
			suNg U_xmu = *pu_gauge(ix,mu); // U_mu(x)
			suNg U_down = *pu_gauge(idn(ix,mu),mu); // U_mu(x-mu)
			suNg_vector Splus = *pu_scalar(iup(ix,mu)); // S(x+mu)
			suNg_vector Sminus = *pu_scalar(idn(ix,mu)); // S(x-mu)
			suNg_vector SSum, Down, Up;
			suNg_vector USStar;

			_suNg_inverse_multiply(Down, U_down, Sminus); // U_mu(x-mu)^+S(x-mu)
			_suNg_multiply(Up, U_xmu, Splus); // U_mu(x)S(x+mu)
			_vector_add_g(SSum, Up, Down);
			vector_star(&USStar, &SSum);
			_vector_sub_assign_g(f, USStar);
		}

		double S2;
		double Msq = par->mass;
		double lambda = par->lambda;
		Msq = Msq*Msq + 8.0;

		suNg_vector S = *pu_scalar(ix);
		suNg_vector SStar;
		vector_star(&SStar, &S);
		_vector_prod_re_g(S2, S, S);
		_vector_lc_add_assign_g(f, Msq, SStar, 2.0*lambda*S2, SStar);
		_vector_mul_add_assign_g(*_FIELD_AT(force,ix), -dt, f);
	}
}

static void outer_product(suNg *u, suNg_vector *v1, suNg_vector *v2)
{
	for(int i = 0; i < NG*NG; i++)
	{
		int row = i/NG;
		int col = i%NG;
		_complex_mul_star(u->c[i], v1->c[row], v2->c[col]);
	}
}

static void force_scalar_g(double dt, void *vpar)
{
	force_scalar_par *par = (force_scalar_par*)vpar;
	suNg_av_field *force = *par->g_momenta;

	_MASTER_FOR(&glattice,ix)
	{
		suNg s1, s2;
		suNg_algebra_vector f;

		for(int mu = 0; mu < 4; mu++)
		{
			outer_product(&s1, pu_scalar(iup(ix,mu)), pu_scalar(ix));
			_suNg_times_suNg(s2, *_4FIELD_AT(u_gauge,ix,mu), s1);
			_fund_algebra_project(f,s2);
			_algebra_vector_mul_add_assign_g(*_4FIELD_AT(force,ix,mu), -2.0*dt, f);
		}
	}

	apply_BCs_on_momentum_field(force);
}

void force_scalar(double dt, void* par)
{
	force_scalar_s(dt, par);
	force_scalar_g(dt, par);
}

#endif //WITH_QUATERNIONS