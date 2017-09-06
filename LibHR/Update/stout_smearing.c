/***************************************************************************\
* Copyright (c) 2017, Jarno Rantaharju, Martin Hansen                       *
* All rights reserved.                                                      * 
\***************************************************************************/

#include "global.h"
#include "representation.h"
#include "utils.h"
#include "spinor_field.h"
#include "suN_repr_func.h"
#include "random.h"
#include "update.h"
#include "random.h"
#include "io.h"
#include "utils.h"
#include "suN.h"
#include "memory.h"
#include "dirac.h"
#include "linear_algebra.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "logger.h"
#include "communications.h"

// Smearing parameters
static double rho_s = -0.2;
static double rho_t = -0.2;
static int taylor_order = 8;

// Global variables
static suNg_field *Sigma = NULL;
static suNg_field *L = NULL;
static suNg_field *Lt = NULL;
static int smear_init = 0;

// Set smearing parameters
void init_smearing(double s, double t)
{
	// Allocate smeared field
	u_gauge_s = alloc_gfield(&glattice);

	// Set parameters
	rho_s = -s;
	rho_t = -t;

	// Log info
	lprintf("STOUT", 10, "Smearing parameter (space,time) = (%1.2f,%1.2f)\n", s, t);
}

// Plaquette of the smeared field
static double plaq_s(int ix, int mu, int nu)
{
	int iy, iz;
	double p;
	suNg *v1, *v2, *v3, *v4;
	suNg w1, w2, w3;

	iy = iup(ix,mu);
	iz = iup(ix,nu);

	v1 = _4FIELD_AT(u_gauge_s,ix,mu);
	v2 = _4FIELD_AT(u_gauge_s,iy,nu);
	v3 = _4FIELD_AT(u_gauge_s,iz,mu);
	v4 = _4FIELD_AT(u_gauge_s,ix,nu);

	_suNg_times_suNg(w1,*v1,*v2);
	_suNg_times_suNg(w2,*v4,*v3);
	_suNg_times_suNg_dagger(w3,w1,w2);
	_suNg_trace_re(p,w3);

#ifdef PLAQ_WEIGHTS
	if(plaq_weight)
	{
		return plaq_weight[ix*16+mu*4+nu]*p;
	}
#endif

	return p;
}

double avr_smeared_plaquette()
{
	double pa = 0.;

	_MASTER_FOR(&glattice,ix)
	{
		pa += plaq_s(ix,1,0);
		pa += plaq_s(ix,2,0);
		pa += plaq_s(ix,2,1);
		pa += plaq_s(ix,3,0);
		pa += plaq_s(ix,3,1);
		pa += plaq_s(ix,3,2);
	}

	global_sum(&pa, 1);

#ifdef BC_T_OPEN
	pa /= 6.0*NG*GLB_VOLUME*(GLB_T-1)/GLB_T;
#else
	pa /= 6.0*NG*GLB_VOLUME;
#endif

	return pa;
}

// Staples used to smear the field
static void projected_stout_staples(int ix, int mu, suNg *v)
{
	int ixpmu, ixpnu, ixmnu, ixpmumnu;
	suNg staple, tr1, tr2, sum;
	suNg_algebra_vector av;
	double rho;

	ixpmu = iup(ix,mu);
	_suNg_zero(sum);

	for(int nu = 0; nu < 4; nu++)
	{
		if(nu == mu)
		{
			continue;
		}

		if(mu == 0 || nu == 0)
		{
			rho = rho_t;
		}
		else
		{
			rho = rho_s;
		}

		ixpnu = iup(ix,nu);
		ixmnu = idn(ix,nu);
		ixpmumnu = idn(ixpmu,nu);

		// Up Staple
		_suNg_times_suNg(tr2,*pu_gauge(ix,nu),*pu_gauge(ixpnu,mu));
		_suNg_dagger(tr1,*pu_gauge(ixpmu,nu));
		_suNg_times_suNg(staple,tr2,tr1);
		_suNg_mul(staple,rho,staple);
		_suNg_add_assign(sum,staple);
		
		// Down Staple
		_suNg_times_suNg(tr2,*pu_gauge(ixmnu,mu),*pu_gauge(ixpmumnu,nu));
		_suNg_dagger(tr1,*pu_gauge(ixmnu,nu));
		_suNg_times_suNg(staple,tr1,tr2);
		_suNg_mul(staple,rho,staple);
		_suNg_add_assign(sum,staple);
	}

	_suNg_times_suNg_dagger(staple,*pu_gauge(ix,mu),sum);
	_fund_algebra_project(av,staple);
	_fund_algebra_represent(*v,av);
}

// Determine order of taylor polynomial
static void estimate_taylor_order()
{
	suNg X, Xk, tmp;
	double error = 1.;
	int k = 1;

	if(rho_s < rho_t)
	{
		_suNg_unit(Xk);
		projected_stout_staples(0,1,&X);
	}
	else
	{
		_suNg_unit(Xk);
		projected_stout_staples(0,0,&X);
	}

	while(error > 1.e-24)
	{
		_suNg_times_suNg(tmp,Xk,X);
		_suNg_mul(Xk,1./k,tmp);
		_suNg_sqnorm(error,Xk);
		k++;
	}

	taylor_order = k;
}

// Taylor expand the exponential of an algebra vector
static void Exp_Taylor(suNg *m, suNg *v, int order)
{
	suNg t1, t2, t3;
	_suNg_unit(t3);
	_suNg_unit(t1);

 	for(int n = 1; n <= order; n++)
	{
		_suNg_mul(t2,1./n,t1);
		_suNg_times_suNg(t1,*v,t2);
		_suNg_add_assign(t3,t1);
	}

	*m = t3;
}

static void init_fields()
{
	if(smear_init == 0)
	{
		Sigma = alloc_gfield(&glattice);
		L = alloc_gfield(&glattice);
		Lt = alloc_gfield(&glattice);
		smear_init = 1;
	}
}

// Calculate smeared gauge field
void smear_gauge_field()
{
#ifdef WITH_SMEARING

	init_fields();
	estimate_taylor_order();

	_MASTER_FOR(&glattice,ix)
	{
		for(int mu = 0; mu < 4; mu++)
		{
			suNg s, t1;
			projected_stout_staples(ix,mu,&s);
			*_4FIELD_AT(Sigma,ix,mu) = s;
			Exp_Taylor(&t1, &s, taylor_order);
	  		_suNg_times_suNg(*_4FIELD_AT(u_gauge_s,ix,mu), t1, *pu_gauge(ix,mu));
		}
	}

	start_gf_sendrecv(u_gauge_s);
	complete_gf_sendrecv(u_gauge_s);

#endif
}

// Calculate force on the original gauge field
// The force on a smeared link is projected into an NxN matrix L
// and each term is reordered in all possible ways
void smeared_gauge_force(suNg_av_field *force, suNg_av_field *momenta)
{
#ifdef WITH_SMEARING

	_MASTER_FOR(&glattice,ix)
	{
		suNg t1, t2, t3, t4;
		suNg l, s;

		for(int mu = 0; mu < 4; mu++)
		{
			suNg s_power[taylor_order];
			suNg_algebra_vector av;
			double coeff = 1;
			
			// Represent force as matrix
			av = *_4FIELD_AT(force,ix,mu);
			_fund_algebra_represent(t1,av);
			_suNg_dagger_times_suNg(l,*_4FIELD_AT(u_gauge_s,ix,mu),t1);

			// U_s = U*exp(M)
			// dU_s/dU L = exp(M) L + U dM/dU exp(M) L
			s = *_4FIELD_AT(Sigma,ix,mu);
			Exp_Taylor(&t1, &s, taylor_order);

			// The first contribution, directly wrt the gauge field
			_suNg_times_suNg(t2,l,t1);
			*_4FIELD_AT(L,ix,mu) = t2;

			// dS/dQ = d exp(M) /dQ * ( U L ), t1 = U L
			_suNg_times_suNg(t1, *pu_gauge(ix,mu), l);

			// Calculate s^{n-1} for n = {1,..,n-1}
			_suNg_unit(s_power[0]);
			s_power[1] = s;
			for(int n = 2; n < taylor_order; n++)
			{
				_suNg_times_suNg(s_power[n], s_power[n-1], s);
			}

			// Accumulate in l
			l = t1;
			for(int n = 2; n <= taylor_order; n++)
			{
				coeff /= n;
				_suNg_zero(t4);
				for(int k = 0; k < n; k++)
				{
					_suNg_times_suNg(t2,s_power[k],t1);
					_suNg_times_suNg(t3,t2,s_power[n-k-1]);
					_suNg_add_assign(t4,t3);
				}
				_suNg_mul(t4,coeff,t4);
				_suNg_add_assign(l,t4);
			}

			// Now the derivative of the projection
			_fund_algebra_project(av,l);
			_fund_algebra_represent(l,av);

			// The derivative wrt to the plaquette
			*_4FIELD_AT(Lt,ix,mu) = l;
		}
	}

	// Communicate field
	start_gf_sendrecv(Lt);
	complete_gf_sendrecv(Lt);

	// Now the derivative of the plaquette
	_MASTER_FOR(&glattice,ix)
	{
		suNg_algebra_vector av;
		suNg w1, w2, sum;
		suNg v1, v2, v3;
		double rho;

		for(int mu = 0; mu < 4; mu++)
		{
			for(int nu = 0; nu < 4; nu++)
			{
				if(mu == nu)
				{
					continue;
				}

				if(mu == 0 || nu == 0)
				{
					rho = rho_t;
				}
				else
				{
					rho = rho_s;
				}

				int o1 = iup(ix,mu);
				int o2 = iup(ix,nu);
				int o3 = idn(ix,nu);
				int o4 = idn(o1,nu);

				_suNg_dagger_times_suNg(w1, *pu_gauge(o3,mu), *pu_gauge(o3,nu));
				_suNg_times_suNg_dagger(w2, *pu_gauge(o1,nu), *pu_gauge(o2,mu));
				_suNg_zero(sum);

				// Part 1
				v1 = *_4FIELD_AT(Lt,ix,mu);
				_suNg_sub_assign(v1, *_4FIELD_AT(Lt,ix,nu));
				_suNg_dagger_times_suNg(v2, *pu_gauge(ix,nu), v1);
				_suNg_times_suNg(v3, w2, v2);
				_suNg_add_assign(sum, v3);

				// Part 2
				v1 = *_4FIELD_AT(Lt,o3,nu);
				_suNg_sub_assign(v1, *_4FIELD_AT(Lt,o3,mu));
				_suNg_times_suNg(v2, v1, *pu_gauge(o3,nu));
				_suNg_times_suNg(v1, *pu_gauge(o3,mu), *pu_gauge(o4,nu));
				_suNg_dagger_times_suNg(v3, v1, v2);
				_suNg_add_assign(sum, v3);

				// Part 3
				_suNg_times_suNg(v1, w1, *_4FIELD_AT(Lt,ix,mu));
				_suNg_times_suNg(v2, *_4FIELD_AT(Lt,o4,nu), w1);
				_suNg_sub_assign(v1, v2);
				_suNg_dagger_times_suNg(v3, *pu_gauge(o4,nu), v1);
				_suNg_add_assign(sum, v3);

				// Part 4
				_suNg_times_suNg(v1, *_4FIELD_AT(Lt,o1,nu), w2);
				_suNg_times_suNg(v2, w2, *_4FIELD_AT(Lt,o2,mu));
				_suNg_sub_assign(v1, v2);
				_suNg_times_suNg_dagger(v3, v1, *pu_gauge(ix,nu));
				_suNg_add_assign(sum, v3);

				// Add contribution
				_suNg_mul(sum, rho, sum);
				_suNg_add_assign(*_4FIELD_AT(L,ix,mu), sum);
			}

			_suNg_times_suNg(v1, *pu_gauge(ix,mu),  *_4FIELD_AT(L,ix,mu));
			_fund_algebra_project(av, v1);
			_algebra_vector_add_assign_g(*_4FIELD_AT(momenta,ix,mu), av);
		}
	}

#endif
}
