/***************************************************************************\
* Copyright (c) 2008, Claudio Pica, Vincent Drach and Ari Hietanen          *
* Copyright (c) 2017, Martin Hansen                                         *
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
#include "clover_tools.h"
#include <string.h>

/* we need to compute  Tr  U(x,mu) g_5*(1-g_mu) chi2 # chi1^+
 * where # indicates the tensor product and Tr is the trace on Lorentz space.
 * the strategy is the following:
 * given the form of g_5(1-g_mu) one can compute only the first two lorentz
 * components of the spinor; so we first apply g_5(1-g_mu) to chi2 to find the first
 * two components; then we multiply these two vectors by U(x,mu) and
 * store the result in p.c[0], p.c[1]; when computing the trace we can factorize p.c[0] and p.c[1]
 * as they both multiply two components of chi1^+; we store these factors in p.c[2] and p.c[3].
 * the tensor product is performed by the macro 
 * _suNf_FMAT(u,p): u = p.c[0] # p.c[2]^+ + p.c[1] # p.c[3]^+
 */

/* these macros use the variables ptmp, p */
#ifdef BC_T_THETA
#define _T_theta_mulc(r) _vector_mulc_f(ptmp,eitheta[0],(r)); (r)=ptmp
#else
#define _T_theta_mulc(r)
#endif
#ifdef BC_X_THETA
#define _X_theta_mulc(r) _vector_mulc_f(ptmp,eitheta[1],(r)); (r)=ptmp
#else
#define _X_theta_mulc(r)
#endif
#ifdef BC_Y_THETA
#define _Y_theta_mulc(r) _vector_mulc_f(ptmp,eitheta[2],(r)); (r)=ptmp
#else
#define _Y_theta_mulc(r)
#endif
#ifdef BC_Z_THETA
#define _Z_theta_mulc(r) _vector_mulc_f(ptmp,eitheta[3],(r)); (r)=ptmp
#else
#define _Z_theta_mulc(r)
#endif

#define _F_DIR0(u,chi1,chi2)				      \
  _vector_add_f(ptmp,(chi2)->c[0],(chi2)->c[2]);		      \
  _suNf_multiply(p.c[0],*(pu_gauge_f(ix,0)),ptmp);		      \
  _T_theta_mulc(p.c[0]);                                      \
  _vector_add_f(ptmp,(chi2)->c[1],(chi2)->c[3]);		      \
  _suNf_multiply(p.c[1],*(pu_gauge_f(ix,0)),ptmp);		      \
  _T_theta_mulc(p.c[1]);                                      \
  _vector_sub_f(p.c[2],(chi1)->c[0],(chi1)->c[2]);	      \
  _vector_sub_f(p.c[3],(chi1)->c[1],(chi1)->c[3]);	      \
  _suNf_FMAT((u),p)

#define _F_DIR1(u,chi1,chi2)				      \
  _vector_i_add_f(ptmp,(chi2)->c[0],(chi2)->c[3]);		      \
  _suNf_multiply(p.c[0],*(pu_gauge_f(ix,1)),ptmp);		      \
  _X_theta_mulc(p.c[0]);                                      \
  _vector_i_add_f(ptmp,(chi2)->c[1],(chi2)->c[2]);		      \
  _suNf_multiply(p.c[1],*(pu_gauge_f(ix,1)),ptmp);		      \
  _X_theta_mulc(p.c[1]);                                      \
  _vector_i_sub_f(p.c[2],(chi1)->c[0],(chi1)->c[3]);	      \
  _vector_i_sub_f(p.c[3],(chi1)->c[1],(chi1)->c[2]);	      \
  _suNf_FMAT((u),p)

#define _F_DIR2(u,chi1,chi2)				      \
  _vector_add_f(ptmp,(chi2)->c[0],(chi2)->c[3]);		      \
  _suNf_multiply(p.c[0],*(pu_gauge_f(ix,2)),ptmp);		      \
  _Y_theta_mulc(p.c[0]);                                      \
  _vector_sub_f(ptmp,(chi2)->c[1],(chi2)->c[2]);		      \
  _suNf_multiply(p.c[1],*(pu_gauge_f(ix,2)),ptmp);		      \
  _Y_theta_mulc(p.c[1]);                                      \
  _vector_sub_f(p.c[2],(chi1)->c[0],(chi1)->c[3]);	      \
  _vector_add_f(p.c[3],(chi1)->c[1],(chi1)->c[2]);	      \
  _suNf_FMAT((u),p)

#define _F_DIR3(u,chi1,chi2)				      \
  _vector_i_add_f(ptmp,(chi2)->c[0],(chi2)->c[2]);		      \
  _suNf_multiply(p.c[0],*(pu_gauge_f(ix,3)),ptmp);		      \
  _Z_theta_mulc(p.c[0]);                                      \
  _vector_i_sub_f(ptmp,(chi2)->c[1],(chi2)->c[3]);		      \
  _suNf_multiply(p.c[1],*(pu_gauge_f(ix,3)),ptmp);		      \
  _Z_theta_mulc(p.c[1]);                                      \
  _vector_i_sub_f(p.c[2],(chi1)->c[0],(chi1)->c[2]);	      \
  _vector_i_add_f(p.c[3],(chi1)->c[1],(chi1)->c[3]);	      \
  _suNf_FMAT((u),p)

static suNg_av_field *force_sum = NULL;

#ifdef WITH_CLOVER

/* ----------------------------------- */
/* CALCULATE FORCE OF THE CLOVER TERM  */
/* ----------------------------------- */
static void g5_sigma(suNf_spinor *s, suNf_spinor *u, int mu, int nu)
{
	if(mu == 0 && nu == 1)
	{
		for(int i = 0; i < NF; i++)
		{
			s->c[0].c[i].re =  u->c[1].c[i].im;
			s->c[1].c[i].re =  u->c[0].c[i].im;
			s->c[2].c[i].re =  u->c[3].c[i].im;
			s->c[3].c[i].re =  u->c[2].c[i].im;
			s->c[0].c[i].im = -u->c[1].c[i].re;
			s->c[1].c[i].im = -u->c[0].c[i].re;
			s->c[2].c[i].im = -u->c[3].c[i].re;
			s->c[3].c[i].im = -u->c[2].c[i].re;
		}
	}
	else if(nu == 0 && mu == 1)
	{
		for(int i = 0; i < NF; i++)
		{
			s->c[0].c[i].re = -u->c[1].c[i].im;
			s->c[1].c[i].re = -u->c[0].c[i].im;
			s->c[2].c[i].re = -u->c[3].c[i].im;
			s->c[3].c[i].re = -u->c[2].c[i].im;
			s->c[0].c[i].im =  u->c[1].c[i].re;
			s->c[1].c[i].im =  u->c[0].c[i].re;
			s->c[2].c[i].im =  u->c[3].c[i].re;
			s->c[3].c[i].im =  u->c[2].c[i].re;
		}
	}
	else if(mu == 0 && nu == 2)
	{
		for(int i = 0; i < NF; i++)
		{
			s->c[0].c[i].re = -u->c[1].c[i].re;
			s->c[1].c[i].re =  u->c[0].c[i].re;
			s->c[2].c[i].re = -u->c[3].c[i].re;
			s->c[3].c[i].re =  u->c[2].c[i].re;
			s->c[0].c[i].im = -u->c[1].c[i].im;
			s->c[1].c[i].im =  u->c[0].c[i].im;
			s->c[2].c[i].im = -u->c[3].c[i].im;
			s->c[3].c[i].im =  u->c[2].c[i].im;
		}
	}
	else if(nu == 0 && mu == 2)
	{
		for(int i = 0; i < NF; i++)
		{
			s->c[0].c[i].re =  u->c[1].c[i].re;
			s->c[1].c[i].re = -u->c[0].c[i].re;
			s->c[2].c[i].re =  u->c[3].c[i].re;
			s->c[3].c[i].re = -u->c[2].c[i].re;
			s->c[0].c[i].im =  u->c[1].c[i].im;
			s->c[1].c[i].im = -u->c[0].c[i].im;
			s->c[2].c[i].im =  u->c[3].c[i].im;
			s->c[3].c[i].im = -u->c[2].c[i].im;
		}
	}
	else if(mu == 0 && nu == 3)
	{
		for(int i = 0; i < NF; i++)
		{
			s->c[0].c[i].re =  u->c[0].c[i].im;
			s->c[1].c[i].re = -u->c[1].c[i].im;
			s->c[2].c[i].re =  u->c[2].c[i].im;
			s->c[3].c[i].re = -u->c[3].c[i].im;
			s->c[0].c[i].im = -u->c[0].c[i].re;
			s->c[1].c[i].im =  u->c[1].c[i].re;
			s->c[2].c[i].im = -u->c[2].c[i].re;
			s->c[3].c[i].im =  u->c[3].c[i].re;
		}
	}
	else if(nu == 0 && mu == 3)
	{
		for(int i = 0; i < NF; i++)
		{
			s->c[0].c[i].re = -u->c[0].c[i].im;
			s->c[1].c[i].re =  u->c[1].c[i].im;
			s->c[2].c[i].re = -u->c[2].c[i].im;
			s->c[3].c[i].re =  u->c[3].c[i].im;
			s->c[0].c[i].im =  u->c[0].c[i].re;
			s->c[1].c[i].im = -u->c[1].c[i].re;
			s->c[2].c[i].im =  u->c[2].c[i].re;
			s->c[3].c[i].im = -u->c[3].c[i].re;
		}
	}
	else if(mu == 1 && nu == 2)
	{
		for(int i = 0; i < NF; i++)
		{
			s->c[0].c[i].re = -u->c[0].c[i].im;
			s->c[1].c[i].re =  u->c[1].c[i].im;
			s->c[2].c[i].re =  u->c[2].c[i].im;
			s->c[3].c[i].re = -u->c[3].c[i].im;
			s->c[0].c[i].im =  u->c[0].c[i].re;
			s->c[1].c[i].im = -u->c[1].c[i].re;
			s->c[2].c[i].im = -u->c[2].c[i].re;
			s->c[3].c[i].im =  u->c[3].c[i].re;
		}
	}
	else if(nu == 1 && mu == 2)
	{
		for(int i = 0; i < NF; i++)
		{
			s->c[0].c[i].re =  u->c[0].c[i].im;
			s->c[1].c[i].re = -u->c[1].c[i].im;
			s->c[2].c[i].re = -u->c[2].c[i].im;
			s->c[3].c[i].re =  u->c[3].c[i].im;
			s->c[0].c[i].im = -u->c[0].c[i].re;
			s->c[1].c[i].im =  u->c[1].c[i].re;
			s->c[2].c[i].im =  u->c[2].c[i].re;
			s->c[3].c[i].im = -u->c[3].c[i].re;
		}
	}
	else if(mu == 1 && nu == 3)
	{
		for(int i = 0; i < NF; i++)
		{
			s->c[0].c[i].re = -u->c[1].c[i].re;
			s->c[1].c[i].re =  u->c[0].c[i].re;
			s->c[2].c[i].re =  u->c[3].c[i].re;
			s->c[3].c[i].re = -u->c[2].c[i].re;
			s->c[0].c[i].im = -u->c[1].c[i].im;
			s->c[1].c[i].im =  u->c[0].c[i].im;
			s->c[2].c[i].im =  u->c[3].c[i].im;
			s->c[3].c[i].im = -u->c[2].c[i].im;
		}
	}
	else if(nu == 1 && mu == 3)
	{
		for(int i = 0; i < NF; i++)
		{
			s->c[0].c[i].re =  u->c[1].c[i].re;
			s->c[1].c[i].re = -u->c[0].c[i].re;
			s->c[2].c[i].re = -u->c[3].c[i].re;
			s->c[3].c[i].re =  u->c[2].c[i].re;
			s->c[0].c[i].im =  u->c[1].c[i].im;
			s->c[1].c[i].im = -u->c[0].c[i].im;
			s->c[2].c[i].im = -u->c[3].c[i].im;
			s->c[3].c[i].im =  u->c[2].c[i].im;
		}
	}
	else if(mu == 2 && nu == 3)
	{
		for(int i = 0; i < NF; i++)
		{
			s->c[0].c[i].re = -u->c[1].c[i].im;
			s->c[1].c[i].re = -u->c[0].c[i].im;
			s->c[2].c[i].re =  u->c[3].c[i].im;
			s->c[3].c[i].re =  u->c[2].c[i].im;
			s->c[0].c[i].im =  u->c[1].c[i].re;
			s->c[1].c[i].im =  u->c[0].c[i].re;
			s->c[2].c[i].im = -u->c[3].c[i].re;
			s->c[3].c[i].im = -u->c[2].c[i].re;
		}
	}
	else if(nu == 2 && mu == 3)
	{
		for(int i = 0; i < NF; i++)
		{
			s->c[0].c[i].re =  u->c[1].c[i].im;
			s->c[1].c[i].re =  u->c[0].c[i].im;
			s->c[2].c[i].re = -u->c[3].c[i].im;
			s->c[3].c[i].re = -u->c[2].c[i].im;
			s->c[0].c[i].im = -u->c[1].c[i].re;
			s->c[1].c[i].im = -u->c[0].c[i].re;
			s->c[2].c[i].im =  u->c[3].c[i].re;
			s->c[3].c[i].im =  u->c[2].c[i].re;
		}
	}
}

static suNf fmat_create(suNf_spinor *a_lhs, suNf_spinor *a_rhs, suNf_spinor *b_lhs, suNf_spinor *b_rhs)
{
	suNf fmat;
	_suNf_zero(fmat);
	for(int i = 0; i < NF; i++)
	{
		for(int j = 0; j < NF; j++)
		{
			for(int k = 0; k < 4; k++)
			{
#ifdef REPR_IS_REAL
				fmat.c[i*NF+j] += a_lhs->c[k].c[i].re*a_rhs->c[k].c[j].re + a_lhs->c[k].c[i].im*a_rhs->c[k].c[j].im;
				fmat.c[i*NF+j] += b_lhs->c[k].c[i].re*b_rhs->c[k].c[j].re + b_lhs->c[k].c[i].im*b_rhs->c[k].c[j].im;
#else
				fmat.c[i*NF+j].re += a_lhs->c[k].c[i].re*a_rhs->c[k].c[j].re + a_lhs->c[k].c[i].im*a_rhs->c[k].c[j].im;
				fmat.c[i*NF+j].im += a_lhs->c[k].c[i].im*a_rhs->c[k].c[j].re - a_lhs->c[k].c[i].re*a_rhs->c[k].c[j].im;
				fmat.c[i*NF+j].re += b_lhs->c[k].c[i].re*b_rhs->c[k].c[j].re + b_lhs->c[k].c[i].im*b_rhs->c[k].c[j].im;
				fmat.c[i*NF+j].im += b_lhs->c[k].c[i].im*b_rhs->c[k].c[j].re - b_lhs->c[k].c[i].re*b_rhs->c[k].c[j].im;
#endif
			}
		}
	}
	return fmat;
}

static void force_clover_core(double dt)
{
	double coeff = dt*(_REPR_NORM2/_FUND_NORM2)*(1./8.)*get_csw();

	// Communicate forces
	start_clover_force_sendrecv(cl_force);

	// Loop over lattice
	_PIECE_FOR(&glattice,xp)
	{
		if(xp == glattice.inner_master_pieces)
		{
			_OMP_PRAGMA(master)
			{
				complete_clover_force_sendrecv(cl_force);
			}
			_OMP_PRAGMA(barrier)
		}

		_SITE_FOR(&glattice,xp,ix)
		{
			suNf *Z[6], W[9];
			suNf s1, s2, s3, fmat;
			suNg_algebra_vector f;
			int num, sign;

			for(int mu = 0; mu < 4; mu++)
			{
				for(int nu = 0; nu < 4; nu++)
				{
					if(mu == nu) continue;
					
					// Coordinates
					int o1 = iup(ix,mu); // x + mu
					int o2 = iup(ix,nu); // x + nu
					int o3 = idn(ix,nu); // x - nu
					int o4 = iup(o3,mu); // x + mu - nu
					int o5 = iup(o2,mu); // x + mu + nu

					if(mu < nu)
					{
						num = nu*(nu-1)/2+mu;
						sign = +1;
					}
					else
					{
						num = mu*(mu-1)/2+nu;
						sign = -1;
					}

					// Force matrices
					Z[0] = _6FIELD_AT(cl_force,ix,num);
					Z[1] = _6FIELD_AT(cl_force,o1,num);
					Z[2] = _6FIELD_AT(cl_force,o3,num);
					Z[3] = _6FIELD_AT(cl_force,o4,num);
					Z[4] = _6FIELD_AT(cl_force,o5,num);
					Z[5] = _6FIELD_AT(cl_force,o2,num);

					// Construct links
					_suNf_dagger(W[0], *pu_gauge_f(o3,mu));
					W[1] = *pu_gauge_f(o3,nu);
					W[2] = *pu_gauge_f(o1,nu);
					_suNf_dagger(W[3], *pu_gauge_f(o2,mu));
					_suNf_dagger(W[4], *pu_gauge_f(ix,nu));
					_suNf_dagger(W[5], *pu_gauge_f(o4,nu));
					_suNf_times_suNf(W[6], W[0], W[1]);
					_suNf_times_suNf(W[7], W[2], W[3]);
					_suNf_times_suNf(s1, W[5], W[6]);
					_suNf_times_suNf(W[8], W[7], W[4]);
					_suNf_sub_assign(W[8], s1);

					// Calculate sum of forces
					_suNf_times_suNf(fmat, W[8], *Z[0]);
					_suNf_times_suNf(s1, *Z[1], W[8]);
					_suNf_add_assign(fmat, s1);
					_suNf_times_suNf(s1, W[0], *Z[2]);
					_suNf_times_suNf(s2, s1, W[1]);
					_suNf_times_suNf(s3, *Z[3], W[6]);
					_suNf_add_assign(s2, s3);
					_suNf_times_suNf(s1, W[5], s2);
					_suNf_sub_assign(fmat, s1);
					_suNf_times_suNf(s1, W[2], *Z[4]);
					_suNf_times_suNf(s2, s1, W[3]);
					_suNf_times_suNf(s3, W[7], *Z[5]);
					_suNf_add_assign(s2, s3);
					_suNf_times_suNf(s1, s2, W[4]);
					_suNf_add_assign(fmat, s1);
					_suNf_times_suNf(s1, *pu_gauge_f(ix,mu), fmat);

					// Project on force
					_algebra_project(f,s1);
					_algebra_vector_mul_add_assign_g(*_4FIELD_AT(force_sum,ix,mu), sign*coeff, f);
				} // nu
			} // mu
		} // sites
	} // pieces
}

void force_clover_fermion(spinor_field* Xs, spinor_field* Ys, double residue)
{
	// Construct force matrices
	_MASTER_FOR(&glattice,ix)
	{
		suNf_spinor tmp_lhs, tmp_rhs;
		suNf_spinor *rhs, *lhs;
		suNf *fm, fm_tmp;

		// (mu,nu) = (0,1)
		rhs = _FIELD_AT(Xs, ix);
		lhs = _FIELD_AT(Ys, ix);
		fm = _6FIELD_AT(cl_force, ix, 0);
		g5_sigma(&tmp_rhs, rhs, 0, 1);
		g5_sigma(&tmp_lhs, lhs, 0, 1);
		fm_tmp = fmat_create(&tmp_lhs, rhs, &tmp_rhs, lhs);
		_suNf_mul(fm_tmp, residue, fm_tmp);
		_suNf_add_assign(*fm, fm_tmp);

		// (mu,nu) = (0,2)
		rhs = _FIELD_AT(Xs, ix);
		lhs = _FIELD_AT(Ys, ix);
		fm = _6FIELD_AT(cl_force, ix, 1);
		g5_sigma(&tmp_rhs, rhs, 0, 2);
		g5_sigma(&tmp_lhs, lhs, 0, 2);
		fm_tmp = fmat_create(&tmp_lhs, rhs, &tmp_rhs, lhs);
		_suNf_mul(fm_tmp, residue, fm_tmp);
		_suNf_add_assign(*fm, fm_tmp);

		// (mu,nu) = (1,2)
		rhs = _FIELD_AT(Xs, ix);
		lhs = _FIELD_AT(Ys, ix);
		fm = _6FIELD_AT(cl_force, ix, 2);
		g5_sigma(&tmp_rhs, rhs, 1, 2);
		g5_sigma(&tmp_lhs, lhs, 1, 2);
		fm_tmp = fmat_create(&tmp_lhs, rhs, &tmp_rhs, lhs);
		_suNf_mul(fm_tmp, residue, fm_tmp);
		_suNf_add_assign(*fm, fm_tmp);
	
		// (mu,nu) = (0,3)
		rhs = _FIELD_AT(Xs, ix);
		lhs = _FIELD_AT(Ys, ix);
		fm = _6FIELD_AT(cl_force, ix, 3);
		g5_sigma(&tmp_rhs, rhs, 0, 3);
		g5_sigma(&tmp_lhs, lhs, 0, 3);
		fm_tmp = fmat_create(&tmp_lhs, rhs, &tmp_rhs, lhs);
		_suNf_mul(fm_tmp, residue, fm_tmp);
		_suNf_add_assign(*fm, fm_tmp);

		// (mu,nu) = (1,3)
		rhs = _FIELD_AT(Xs, ix);
		lhs = _FIELD_AT(Ys, ix);
		fm = _6FIELD_AT(cl_force, ix, 4);
		g5_sigma(&tmp_rhs, rhs, 1, 3);
		g5_sigma(&tmp_lhs, lhs, 1, 3);
		fm_tmp = fmat_create(&tmp_lhs, rhs, &tmp_rhs, lhs);
		_suNf_mul(fm_tmp, residue, fm_tmp);
		_suNf_add_assign(*fm, fm_tmp);

		// (mu,nu) = (2,3)
		rhs = _FIELD_AT(Xs, ix);
		lhs = _FIELD_AT(Ys, ix);
		fm = _6FIELD_AT(cl_force, ix, 5);
		g5_sigma(&tmp_rhs, rhs, 2, 3);
		g5_sigma(&tmp_lhs, lhs, 2, 3);
		fm_tmp = fmat_create(&tmp_lhs, rhs, &tmp_rhs, lhs);
		_suNf_mul(fm_tmp, residue, fm_tmp);
		_suNf_add_assign(*fm, fm_tmp);
	}
}

void force_clover_logdet(double mass, double residue)
{
	// Compute force matrices
	compute_force_logdet(mass, residue);
}

#endif //#ifdef WITH_CLOVER

/* ------------------------------------ */
/* CALCULATE FORCE OF THE HOPPING TERM  */
/* ------------------------------------ */
void force_fermion_core(spinor_field *Xs, spinor_field *Ys, int auto_fill_odd, double dt, double residue)
{
	double coeff;
	spinor_field Xtmp, Ytmp;

	coeff = residue * dt * (_REPR_NORM2/_FUND_NORM2);
	Xtmp = *Xs;
	Ytmp = *Ys;
	Xs->type = &glattice;
	Ys->type = &glattice;

#ifdef UPDATE_EO

	if(auto_fill_odd)
	{
		spinor_field Xe, Xo, Ye, Yo;

		Xe = *Xs;
		Xe.type = &glat_even;
		Xo = *Xs;
		Xo.ptr = Xs->ptr + glat_odd.master_shift;
		Xo.type = &glat_odd;

		Ye = *Ys;
		Ye.type = &glat_even;
		Yo = *Ys;
		Yo.type = &glat_odd;
		Yo.ptr = Ys->ptr + glat_odd.master_shift;

		Dphi_(&Xo, &Xe);
		Dphi_(&Yo, &Ye);
#ifdef WITH_CLOVER
		Cphi_diag_inv(get_dirac_mass(), &Xo, &Xo);
		Cphi_diag_inv(get_dirac_mass(), &Yo, &Yo);
#endif
	}

	coeff = -coeff;

#endif

	// Communicate spinor field
	start_sf_sendrecv(Xs);
	start_sf_sendrecv(Ys);

#ifdef WITH_CLOVER
	force_clover_fermion(Xs, Ys, residue);
#endif

	// Loop over lattice
	_PIECE_FOR(&glattice,xp)
	{
		suNg_algebra_vector f;
		suNf_vector ptmp;
		suNf_spinor p;
		suNf_FMAT s1;

		if(xp == glattice.inner_master_pieces)
		{
			_OMP_PRAGMA(master)
			{
				complete_sf_sendrecv(Xs);
				complete_sf_sendrecv(Ys);
			}
			_OMP_PRAGMA(barrier)
		}

		_SITE_FOR(&glattice,xp,ix)
		{
			int iy;
			suNf_spinor *chi1, *chi2;

			// Direction 0
			iy = iup(ix,0);
			_suNf_FMAT_zero(s1);
			chi1 = _FIELD_AT(Xs,ix);
			chi2 = _FIELD_AT(Ys,iy);
			_F_DIR0(s1,chi1,chi2);
			chi1 = _FIELD_AT(Ys,ix);
			chi2 = _FIELD_AT(Xs,iy);
			_F_DIR0(s1,chi1,chi2);

			_algebra_project_FMAT(f,s1);
			_algebra_vector_mul_add_assign_g(*_4FIELD_AT(force_sum,ix,0),coeff,f);

			// Direction 1
			iy = iup(ix,1);
			_suNf_FMAT_zero(s1);
			chi1 = _FIELD_AT(Xs,ix);
			chi2 = _FIELD_AT(Ys,iy);
			_F_DIR1(s1,chi1,chi2);
			chi1 = _FIELD_AT(Ys,ix);
			chi2 = _FIELD_AT(Xs,iy);
			_F_DIR1(s1,chi1,chi2);

			_algebra_project_FMAT(f,s1);
			_algebra_vector_mul_add_assign_g(*_4FIELD_AT(force_sum,ix,1),coeff,f);

			// Direction 2
			iy = iup(ix,2);
			_suNf_FMAT_zero(s1);
			chi1 = _FIELD_AT(Xs,ix);
			chi2 = _FIELD_AT(Ys,iy);
			_F_DIR2(s1,chi1,chi2);
			chi1 = _FIELD_AT(Ys,ix);
			chi2 = _FIELD_AT(Xs,iy);
			_F_DIR2(s1,chi1,chi2);

			_algebra_project_FMAT(f,s1);
			_algebra_vector_mul_add_assign_g(*_4FIELD_AT(force_sum,ix,2),coeff,f);

			// Direction 3
			iy = iup(ix,3);
			_suNf_FMAT_zero(s1);
			chi1 = _FIELD_AT(Xs,ix);
			chi2 = _FIELD_AT(Ys,iy);
			_F_DIR3(s1,chi1,chi2);
			chi1 = _FIELD_AT(Ys,ix);
			chi2 = _FIELD_AT(Xs,iy);
			_F_DIR3(s1,chi1,chi2);

			_algebra_project_FMAT(f,s1);
			_algebra_vector_mul_add_assign_g(*_4FIELD_AT(force_sum,ix,3),coeff,f);
		} // sites
	} // pieces

	// Reset spinor geometry
	Xs->type = Xtmp.type;
	Ys->type = Ytmp.type;
}

#undef _F_DIR0
#undef _F_DIR1
#undef _F_DIR2
#undef _F_DIR3

void fermion_force_begin()
{
	if(force_sum == NULL)
	{
		force_sum = alloc_avfield(&glattice);
	}

	// Clear force field
	_MASTER_FOR(&glattice,ix)
	{
		for(int mu = 0; mu < 4; mu++)
		{
			_algebra_vector_zero_g(*_4FIELD_AT(force_sum,ix,mu));
		}
	}

#ifdef WITH_CLOVER

	// Clear clover force field
	_MASTER_FOR(&glattice,ix)
	{
		for(int mu = 0; mu < 6; mu++)
		{
			_suNf_zero(*_6FIELD_AT(cl_force,ix,mu));
		}
	}

#endif
}

void fermion_force_end(double dt, suNg_av_field *force)
{
#ifdef WITH_CLOVER

	// Evaluate derivative of clover term
	force_clover_core(dt);

#endif

#ifdef WITH_SMEARING

	// Evaluate smeared force and add to global force field
	smeared_gauge_force(force_sum, force);

#else

	// Add force to global force field
	_MASTER_FOR(&glattice,ix)
	{
		for(int mu = 0; mu < 4; mu++)
		{
			_algebra_vector_add_assign_g(*_4FIELD_AT(force,ix,mu), *_4FIELD_AT(force_sum,ix,mu));
		}
	}

#endif

	// Boundary conditions
	apply_BCs_on_momentum_field(force);
}
