/***************************************************************************\
* Copyright (c) 2008, Claudio Pica, Vincent Drach and Ari Hietanen          *
* Copyright (c) 2016, Martin Hansen                                         *
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

#ifdef MEASURE_FORCE
static suNg_av_field *force_tmp = NULL;
static int force_init = 0;
#endif

void force_measure_begin()
{
#ifdef MEASURE_FORCE
	if(force_init == 0)
	{
		force_tmp = alloc_avfield(&glattice);
		force_init = 1;
	}
	_MASTER_FOR(&glattice,ix)
	{
		for(int mu = 0; mu < 4; mu++)
		{
			_algebra_vector_zero_g(*_4FIELD_AT(force_tmp,ix,mu));
		}
	}
#endif
}

void force_measure_end(int id, const char *name, double dt, int nit)
{
#ifdef MEASURE_FORCE
	double max = 0;
	double sum = 0;

	// This loop does not work with OpenMP
	_MASTER_FOR(&glattice,ix)
	{
		suNg_algebra_vector *f;
		double nsq;

		for(int mu = 0; mu < 4; mu++)
		{
			// Calculate |F|^2 = sum_{a,b} (F^a * F^b) * (T_R * delta^{ab})
			f = _4FIELD_AT(force_tmp, ix, mu);
			_algebra_vector_sqnorm_g(nsq, *f);
			nsq *= _REPR_NORM2;

			sum += nsq;
			if(nsq > max) max = nsq;
		}
	}

	global_max(&max, 1);
	global_sum(&sum, 1);

	force_ave[id] += sum;
	force_max[id] += max;
	n_inv_iter[id-1] += nit;

	lprintf("FORCE", 20, "%s: id = %d, sum|F|^2 = %1.8e, avg|F|^2 = %1.8e, max|F|^2 = %1.8e, dt = %1.8e\n", name, id, sum, sum/(4*GLB_VOLUME), max, dt);
#endif
}

void force_fermion_core(spinor_field *Xs, spinor_field *Ys, suNg_av_field *force, int auto_fill_odd, double dt, double residue)
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
	}

	coeff = -coeff;

#endif

	// Communicate spinor field
	start_sf_sendrecv(Xs);
	start_sf_sendrecv(Ys);

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

			_algebra_project(f,s1);
			_algebra_vector_mul_g(f,coeff,f);
			_algebra_vector_add_assign_g(*_4FIELD_AT(force,ix,0),f);
			#ifdef MEASURE_FORCE
			_algebra_vector_add_assign_g(*_4FIELD_AT(force_tmp,ix,0),f);
			#endif

			// Direction 1
			iy = iup(ix,1);
			_suNf_FMAT_zero(s1);
			chi1 = _FIELD_AT(Xs,ix);
			chi2 = _FIELD_AT(Ys,iy);
			_F_DIR1(s1,chi1,chi2);
			chi1 = _FIELD_AT(Ys,ix);
			chi2 = _FIELD_AT(Xs,iy);
			_F_DIR1(s1,chi1,chi2);

			_algebra_project(f,s1);
			_algebra_vector_mul_g(f,coeff,f);
			_algebra_vector_add_assign_g(*_4FIELD_AT(force,ix,1),f);
			#ifdef MEASURE_FORCE
			_algebra_vector_add_assign_g(*_4FIELD_AT(force_tmp,ix,1),f);
			#endif

			// Direction 2
			iy = iup(ix,2);
			_suNf_FMAT_zero(s1);
			chi1 = _FIELD_AT(Xs,ix);
			chi2 = _FIELD_AT(Ys,iy);
			_F_DIR2(s1,chi1,chi2);
			chi1 = _FIELD_AT(Ys,ix);
			chi2 = _FIELD_AT(Xs,iy);
			_F_DIR2(s1,chi1,chi2);

			_algebra_project(f,s1);
			_algebra_vector_mul_g(f,coeff,f);
			_algebra_vector_add_assign_g(*_4FIELD_AT(force,ix,2),f);
			#ifdef MEASURE_FORCE
			_algebra_vector_add_assign_g(*_4FIELD_AT(force_tmp,ix,2),f);
			#endif

			// Direction 3
			iy = iup(ix,3);
			_suNf_FMAT_zero(s1);
			chi1 = _FIELD_AT(Xs,ix);
			chi2 = _FIELD_AT(Ys,iy);
			_F_DIR3(s1,chi1,chi2);
			chi1 = _FIELD_AT(Ys,ix);
			chi2 = _FIELD_AT(Xs,iy);
			_F_DIR3(s1,chi1,chi2);

			_algebra_project(f,s1);
			_algebra_vector_mul_g(f,coeff,f);
			_algebra_vector_add_assign_g(*_4FIELD_AT(force,ix,3),f);
			#ifdef MEASURE_FORCE
			_algebra_vector_add_assign_g(*_4FIELD_AT(force_tmp,ix,3),f);
			#endif
		} // sites
	} // pieces

	// Boundary conditions
	apply_BCs_on_momentum_field(force);

	// Reset spinor geometry
	Xs->type = Xtmp.type;
	Ys->type = Ytmp.type;
}

#undef _F_DIR0
#undef _F_DIR1
#undef _F_DIR2
#undef _F_DIR3
