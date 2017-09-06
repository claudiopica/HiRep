/***************************************************************************\
 * Copyright (c) 2017, Agostino Patella, Martin Hansen                     *
 * All rights reserved.                                                    *
\***************************************************************************/

/*

Main functions:

===  double lw_action(double beta, double c0, double c1);
  Returns the LW action, already summed over the global volume.

===  void lw_force(suNg_av_field *force, double beta, double c0, double c1);
  It sets force equal to the LW force.


Core functions:

===  static void calculate_stfld(int comm);
  It calculates the staples in the local lattice and organizes them in the field stfld. If comm==true then it also
  retrieves the value of the staples in the communication buffers from neighbouring processes. Notice that not all
  communication buffers are actually filled but only the ones that are needed for the force calculation. IMPORTANT:
  This function assume that the communication buffers of the gauge field have been already filled.

===  static double lw_action_density(int ix, double beta, double c0, double c1);
  It calculates the LW action density. IMPORTANT: This function assumes that the staples have been already calculated.


Test functions:

=== void test_wilson_action_and_force(double beta);
  It tests the LW action and force with c0=1 and c1=0 against Wilson action and force already present in the code.

=== void test_ginv_lw_action(double beta, double c0, double c1);
  It tests the gauge invariance of the LW action.

=== void test_gcov_lw_force(double beta, double c0, double c1);
  It tests the gauge covariance of the LW force.

=== void test_lw_force(double beta, double c0, double c1) {
  It tests the LW force against the numerical derivative of the LW action, point by point in the global lattice.

===  static void random_g(suNg_field* g);
===  static void transform(suNg_field* gtransf, suNg_field* gfield);
===  static void transform_force(suNg_field* gtransf, suNg_av_field* force);
  Utilities for gauge transformations.

*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "global.h"
#include "communications.h"
#include "memory.h"
#include "update.h"
#include "suN.h"
#include "utils.h"
#include "representation.h"
#include "logger.h"
#include "observables.h"
#include "random.h"

#define COMM (1==1)
#define NOCOMM (1==0)

/*********************************************************************
mu = (nu+i+1)&0x3;
stfld[2*nu+0][3*ix+i] = ix -> ix-nu -> ix-nu+mu -> ix+mu
stfld[2*nu+1][3*ix+i] = ix -> ix+nu -> ix+nu+mu -> ix+mu
*********************************************************************/
static suNg* stfld[8] = {NULL};

static void calculate_stfld(int comm)
{
	suNg *staple, wu1;
	static int nb[8];
	static int source[8][32];
	static int *ib[8] = {NULL};

	if(stfld[0] == NULL)
	{
		stfld[0] = amalloc(glattice.gsize_gauge*sizeof(suNg)*8*3,ALIGN);
		for(int k = 1; k < 8; k++)
		{
			stfld[k] = stfld[k-1] + 3*glattice.gsize_gauge;
		}

		ib[0] = malloc(sizeof(int)*glattice.nbuffers_gauge*8);
		for(int k = 0; k < 8; k++)
		{
			nb[k] = glattice.nbuffers_gauge;
			ib[k] = ib[0] + k*glattice.nbuffers_gauge;
			for(int i = 0; i < glattice.nbuffers_gauge; i++)
			{
				ib[k][i] = i;
			}
		}

		for(int nu = 0; nu < 4; nu++)
		{
			source[2*nu+0][0] = proc_dn(CID,nu);
			for(int i = 0; i < 3; i++)
			{
				int mu = (nu+i+1)&0x3;
				source[2*nu+0][i+1] = proc_dn(source[2*nu+0][0],mu);
				source[2*nu+0][i+4] = proc_dn(CID,mu);
			}

			source[2*nu+1][0] = proc_up(CID,nu);
			for(int i = 0; i < 3; i++)
			{
				int mu = (nu+i+1)&0x3;
				source[2*nu+1][i+1] = proc_dn(source[2*nu+1][0],mu);
				source[2*nu+1][i+4] = proc_dn(CID,mu);
			}
		}

		for(int k = 0; k < 8; k++)
		{
			nb[k] = 0;
			for(int n = 0; n < glattice.nbuffers_gauge; n++)
			{
				for(int i = 0; i < 7; i++)
				{
					if(glattice.rbuf_from_proc[n] == source[k][i])
					{
						ib[k][nb[k]] = n;
						nb[k]++;
						break;
					}
				}
			}
		}

		lprintf("INIT", 0, "nbuffers_gauge = %d\n", glattice.nbuffers_gauge);
		for(int k = 0; k < 8; k++)
		{
			lprintf("INIT", 0, "nb[%d] = %d\n", k, nb[k]);
		}
	}

	memset(stfld[0], 0, glattice.gsize_gauge*sizeof(suNg)*8*3);
	//start_gf_sendrecv(u_gauge);
	//complete_gf_sendrecv(u_gauge);
	//apply_BCs_on_fundamental_gauge_field();

	_MASTER_FOR(&glattice,ix)
	{
		for(int nu = 0; nu < 4; nu++)
		{
			int ixpnu = iup(ix,nu);
			int ixmnu = idn(ix,nu);

			for(int i = 0; i < 3; i++)
			{
				int mu = (nu+i+1)&0x3;
				int ixpmu = iup(ix,mu);
				int ixpmumnu = idn(ixpmu,nu);

				// *---
				// |
				// *---
				staple = &stfld[2*nu+0][3*ix+i];
				_suNg_times_suNg(wu1, *pu_gauge(ixmnu,mu), *pu_gauge(ixpmumnu,nu));
				_suNg_dagger_times_suNg(*staple, *pu_gauge(ixmnu,nu), wu1);

				// ---*
				//    |
				// ---*
				staple = &stfld[2*nu+1][3*ix+i];
				_suNg_times_suNg(wu1, *pu_gauge(ix,nu), *pu_gauge(ixpnu,mu));
				_suNg_times_suNg_dagger(*staple, wu1, *pu_gauge(ixpmu,nu));
			}
		}
	}

	if(comm)
	{
		for(int k = 0; k < 8; k++)
		{
			for(int i = 0; i < glattice.ncopies_gauge; i++)
			{
				memcpy(stfld[k]+3*glattice.copy_to[i], stfld[k]+3*glattice.copy_from[i], 3*glattice.copy_len[i]*sizeof(suNg));
			}
		}

		for(int k = 0; k < 8; k++)
		{
			for(int i = 0; i < nb[k]; i++)
			{
				int n = ib[k][i];
				#ifdef WITH_MPI
				MPI_Sendrecv(stfld[k]+3*glattice.sbuf_start[n],
									glattice.sbuf_len[n]*sizeof(suNg)/sizeof(double)*3,
									MPI_DOUBLE,
									glattice.sbuf_to_proc[n],
									4*k+i,
									stfld[k]+3*glattice.rbuf_start[n],
									glattice.rbuf_len[n]*sizeof(suNg)/sizeof(double)*3,
									MPI_DOUBLE,
									glattice.rbuf_from_proc[n],
									4*k+i,
									cart_comm,
									MPI_STATUS_IGNORE);
				#else
				memcpy(stfld[k]+3*glattice.rbuf_start[n], stfld[k]+3*glattice.sbuf_start[n], glattice.sbuf_len[n]*sizeof(suNg)*3);
				#endif
			}
		}
	}
}

static double lw_action_density(int ix, double beta, double c0, double c1)
{
	double plaqs = 0;
	double rects = 0;
	double p;
	suNg w1;

	for(int nu = 0; nu < 3; nu++)
	{
		for(int mu = nu+1; mu < 4; mu++)
		{
			int i = mu-nu-1;
			_suNg_times_suNg_dagger(w1, stfld[2*nu+1][3*ix+i], *pu_gauge(ix,mu));
			_suNg_trace_re(p,w1);
#ifdef PLAQ_WEIGHTS
			p *= plaq_weight[16*ix+4*mu+nu];
#endif
			plaqs -= p;
		}
	}

	for(int nu = 0; nu < 4; nu++)
	{
		for(int i = 0; i < 3; i++)
		{
			int ixpnu = ix;
			_suNg_times_suNg_dagger(w1, stfld[2*nu+1][3*ixpnu+i], stfld[2*nu+0][3*ixpnu+i]);
			_suNg_trace_re(p,w1);
#ifdef PLAQ_WEIGHTS
			int mu = (nu+i+1)&0x3;
			ixpnu = idn(ix,nu);
			p *= rect_weight[16*ixpnu+4*mu+nu];
#endif
			rects -= p;
		}
	}

	return (beta/NG)*(c0*plaqs+c1*rects);
}

double lw_action(double beta, double c0, double c1)
{
	double s = 0;
	calculate_stfld(NOCOMM);

	_MASTER_FOR(&glattice,ix)
	{
		s += lw_action_density(ix,beta,c0,c1);
	}

	global_sum(&s,1);
	return s;
}

void lw_local_action(scalar_field *loc_action, double beta, double c0, double c1)
{
	calculate_stfld(COMM);
	int iy = ipt(2,0,0,0);
	_MASTER_FOR(&glattice,ix)
	{
		*_FIELD_AT(loc_action,iy) += lw_action_density(ix,beta,c0,c1);
	}
}

void lw_force(double dt, void *vpar)
{
	suNg ws[4], wu1, wu2;
	suNg_algebra_vector wf1;

	force_gauge_par *par = (force_gauge_par*)vpar;
	suNg_av_field *force = *par->momenta;
	double beta = par->beta;
	double c0 = par->c0;
	double c1 = par->c1;

	calculate_stfld(COMM);

	// Calculation of the force in (ix,mu).
	// In the drawings below, mu is the direction of the missing link.
	// The index nu labels the directions orthogonal to mu.
	_MASTER_FOR(&glattice,ix)
	{
		for(int mu = 0; mu < 4; mu++)
		{
			_suNg_zero(ws[mu]);
		}

		for(int nu = 0; nu < 4; nu++)
		{
			for(int i = 0; i < 3; i++)
			{
				int mu = (nu+i+1)&0x3;

				// *---
				// |
				// *---
#ifdef PLAQ_WEIGHTS
				int ixmnu = idn(ix,nu);
				_suNg_mul(wu1, plaq_weight[16*ixmnu+4*mu+nu], stfld[2*nu+0][3*ix+i]);
				_suNg_add_assign(ws[mu], wu1);
#else
				_suNg_add_assign(ws[mu], stfld[2*nu+0][3*ix+i]);
#endif

				// ---*
				//    |
				// ---*
#ifdef PLAQ_WEIGHTS
				_suNg_mul(wu1, plaq_weight[16*ix+4*mu+nu], stfld[2*nu+1][3*ix+i]);
				_suNg_add_assign(ws[mu], wu1);
#else
				_suNg_add_assign(ws[mu], stfld[2*nu+1][3*ix+i]);
#endif
			}
		}

		for(int mu = 0; mu < 4; mu++)
		{
			_suNg_times_suNg_dagger(wu1, *pu_gauge(ix,mu), ws[mu]);
			_fund_algebra_project(wf1, wu1);
			_algebra_vector_mul_add_assign_g(*_4FIELD_AT(force,ix,mu), dt*(-beta*c0/NG), wf1);
			_suNg_zero(ws[mu]);
		}

		for(int nu = 0; nu < 4; nu++)
		{
			int ixpnu = iup(ix,nu);
			int ixmnu = idn(ix,nu);

			for(int i = 0; i < 3; i++)
			{
				int mu = (nu+i+1)&0x3;
				int ixpmu = iup(ix,mu);
				int ixpmunnu = idn(ixpmu,nu);

				// *---*---
				// |
				// *---*---
				_suNg_dagger_times_suNg(wu1, *pu_gauge(ixmnu,nu), stfld[2*nu+0][3*ixmnu+i]);
				_suNg_times_suNg(wu2, wu1, *pu_gauge(ixpmunnu,nu));
#ifdef PLAQ_WEIGHTS
				int ixmnumnu = idn(ixmnu,nu);
				_suNg_mul(wu2, rect_weight[16*ixmnumnu+4*mu+nu], wu2);
#endif
				_suNg_add_assign(ws[mu], wu2);

				// ---*---*
				//        |
				// ---*---*
				_suNg_times_suNg(wu1, *pu_gauge(ix,nu), stfld[2*nu+1][3*ixpnu+i]);
				_suNg_times_suNg_dagger(wu2, wu1, *pu_gauge(ixpmu,nu));
#ifdef PLAQ_WEIGHTS
				_suNg_mul(wu2, rect_weight[16*ix+4*mu+nu], wu2);
#endif
				_suNg_add_assign(ws[mu], wu2);
			}
		}

		for(int mu = 0; mu < 4; mu++)
		{
			int ixpmu = iup(ix,mu);
			int ixmmu = idn(ix,mu);

			for(int i = 0; i < 3; i++)
			{
				int nu = (mu+i+1)&0x3;
				int ixpnu = iup(ix,nu);
				int ixmnu = idn(ix,nu);
				int ixmnupmu = iup(ixmnu,mu);
				int ixmnummu = idn(ixmnu,mu);

				// *---*
				// |   |
				// *   *
				// |
				// *---*
				_suNg_dagger_times_suNg(wu1, *pu_gauge(ixmnu,nu), *pu_gauge(ixmnu,mu));
				_suNg_times_suNg(wu2, wu1, stfld[2*mu+1][3*ixmnupmu+i]);
#ifdef PLAQ_WEIGHTS
				_suNg_mul(wu2, rect_weight[16*ixmnu+4*nu+mu], wu2);
#endif
				_suNg_add_assign(ws[mu], wu2);

				// *---*
				// |   |
				// *   *
				//     |
				// *---*
				_suNg_times_suNg(wu1, *pu_gauge(ix,nu), *pu_gauge(ixpnu,mu));
				_suNg_times_suNg_dagger(wu2, wu1, stfld[2*mu+1][3*ixpmu+i]);
#ifdef PLAQ_WEIGHTS
				_suNg_mul(wu2, rect_weight[16*ix+4*nu+mu], wu2);
#endif
				_suNg_add_assign(ws[mu], wu2);

				// *---*
				// |
				// *   *
				// |   |
				// *---*
				_suNg_dagger_times_suNg(wu1, stfld[2*mu+0][3*ixmnu+i], *pu_gauge(ixmnu,mu));
				_suNg_times_suNg(wu2, wu1, *pu_gauge(ixmnupmu,nu));
#ifdef PLAQ_WEIGHTS
				_suNg_mul(wu2,rect_weight[16*ixmnummu+4*nu+mu],wu2);
#endif
				_suNg_add_assign(ws[mu],wu2);

				// *---*
				//     |
				// *   *
				// |   |
				// *---*
				_suNg_times_suNg(wu1, stfld[2*mu+0][3*ix+i], *pu_gauge(ixpnu,mu));
				_suNg_times_suNg_dagger(wu2, wu1, *pu_gauge(ixpmu,nu));
#ifdef PLAQ_WEIGHTS
				_suNg_mul(wu2, rect_weight[16*ixmmu+4*nu+mu], wu2);
#endif
				_suNg_add_assign(ws[mu], wu2);
			}
		}

		for(int mu = 0; mu < 4; mu++)
		{
			_suNg_times_suNg_dagger(wu1, *pu_gauge(ix,mu), ws[mu]);
			_fund_algebra_project(wf1, wu1);
			_algebra_vector_mul_add_assign_g(*_4FIELD_AT(force,ix,mu), dt*(-beta*c1/NG), wf1);
		}
	}

	apply_BCs_on_momentum_field(force);
}

/*
void test_wilson_action_and_force(double beta) {
   double s1,s2,diff,err;
   suNg_av_field *f1,*f2;
   suNg_algebra_vector *v1, *v2;
	mon_lw_par par;

   random_u(u_gauge);
   start_gf_sendrecv(u_gauge);
   complete_gf_sendrecv(u_gauge);

   calculate_stfld(NOCOMM);

   err=0.;
   _MASTER_FOR(&glattice,ix) {
     s1=(beta/((double)NG))*(6.*NG-local_plaq(ix));
     s2=lw_action_density(ix,beta,1.0,0.0);
     diff=fabs(s1-s2);
     if(diff>err) err=diff;
   }
   global_max(&err,1);

   lprintf("TEST",0,"beta=%f :  Error Wilson action density = %e\n",beta,err);

   f1=alloc_avfield(&glattice);
   f2=alloc_avfield(&glattice);
   _MASTER_FOR(&glattice,ix) {
      for(int mu=0; mu<4; mu++) {
         _algebra_vector_zero_g(*_4FIELD_AT(f1,ix,mu));
         _algebra_vector_zero_g(*_4FIELD_AT(f2,ix,mu));
      }
   }
   force0(1.,f1,&beta);
	par.c0 = 1;
	par.c1 = 0;
	par.beta = beta;
   lw_force(1.,f2,&par);

   err=0.;
   _MASTER_FOR(&glattice,ix) {
      for(int mu=0; mu<4; mu++) {
         v1=_4FIELD_AT(f1,ix,mu);
         v2=_4FIELD_AT(f2,ix,mu);
         for(int i=0;i<NG*NG-1;i++) {
            diff=fabs(v1->c[i]-v2->c[i]);
            if(diff>err) err=diff;
         }
      }
   }
   global_max(&err,1);

   lprintf("TEST",0,"beta=%f :  Error Wilson force = %e\n",beta,err);

   free_avfield(f1);
   free_avfield(f2);
}


static void random_g(suNg_field* g)
{
  _MASTER_FOR(&glattice,ix) {
    random_suNg(_FIELD_AT(g,ix));
  }
  start_gt_sendrecv(g);
  complete_gt_sendrecv(g);
}

static void transform(suNg_field* gtransf, suNg_field* gfield)
{
  _MASTER_FOR(&glattice,ix) {
    for (int mu=0;mu<4;mu++) {
      int iy=iup(ix,mu);
      suNg * u=_4FIELD_AT(gfield,ix,mu);
      suNg v;
      _suNg_times_suNg_dagger(v,*u,*_FIELD_AT(gtransf,iy));
      _suNg_times_suNg(*u,*_FIELD_AT(gtransf,ix),v);
    }
  }
  start_gf_sendrecv(u_gauge);
  complete_gf_sendrecv(u_gauge);
}

static void transform_force(suNg_field* gtransf, suNg_av_field* force)
{
  _MASTER_FOR(&glattice,ix) {
    for (int mu=0;mu<4;mu++) {
      suNg v,m;
      _fund_algebra_represent(m,*_4FIELD_AT(force,ix,mu));
      _suNg_times_suNg_dagger(v,m,*_FIELD_AT(gtransf,ix));
      _suNg_times_suNg(m,*_FIELD_AT(gtransf,ix),v);
      _fund_algebra_project(*_4FIELD_AT(force,ix,mu),m);
    }
  }
 }


void test_ginv_lw_action(double beta, double c0, double c1) {
   suNg_field *g;
   double *s;
   double diff, err;

   g=alloc_gtransf(&glattice);
   s=malloc(glattice.gsize_gauge*sizeof(double));

   random_u(u_gauge);
   start_gf_sendrecv(u_gauge);
   complete_gf_sendrecv(u_gauge);

   calculate_stfld(NOCOMM);
   _MASTER_FOR(&glattice,ix) {
      s[ix]=lw_action_density(ix,beta,c0,c1);
   }

   random_g(g);
   transform(g,u_gauge);

   calculate_stfld(NOCOMM);
   err=0.;
   _MASTER_FOR(&glattice,ix) {
     diff=fabs(s[ix]-lw_action_density(ix,beta,c0,c1));
     if(diff>err) err=diff;
   }
   global_max(&err,1);

   lprintf("TEST",0,"pars=(%f,%f,%f) :  Gauge invariance LW action = %e\n",beta,c0,c1,err);

   free(s);
   free_gtransf(g);
}


void test_gcov_lw_force(double beta, double c0, double c1) {
   suNg_field *g;
   suNg_av_field *f1, *f2;
   suNg_algebra_vector *v1, *v2;
   double diff, err;

	mon_lw_par par;
	par.c0 = c0;
	par.c1 = c1;
	par.beta = beta;

   g=alloc_gtransf(&glattice);
   f1=alloc_avfield(&glattice);
   f2=alloc_avfield(&glattice);
	_MASTER_FOR(&glattice,ix) {
		for(int mu=0; mu<4; mu++) {
			_algebra_vector_zero_g(*_4FIELD_AT(f1,ix,mu));
			_algebra_vector_zero_g(*_4FIELD_AT(f2,ix,mu));
		}
	}
	random_u(u_gauge);
   start_gf_sendrecv(u_gauge);
   complete_gf_sendrecv(u_gauge);

   lw_force(1.,f1,&par);

   random_g(g);
   transform_force(g,f1);
   transform(g,u_gauge);

   lw_force(1.,f2,&par);

   err=0.;
   _MASTER_FOR(&glattice,ix) {
      for(int mu=0; mu<4; mu++) {
         v1=_4FIELD_AT(f1,ix,mu);
         v2=_4FIELD_AT(f2,ix,mu);
         double loc=0.;
         for(int i=0;i<NG*NG-1;i++) {
            diff=fabs(v1->c[i]-v2->c[i]);
            if(diff>err) loc=diff;
         }
         if(loc>err) err=loc;
         // lprintf("TEST",0,"ix=%d mu=%d   err = %e\n",ix,mu,loc);
      }
   }
   global_max(&err,1);

   lprintf("TEST",0,"pars=(%f,%f,%f) :  Gauge covariance LW force = %e\n",beta,c0,c1,err);

   free_avfield(f1);
   free_avfield(f2);
   free_gtransf(g);
}


void test_lw_force(double beta, double c0, double c1) {
   suNg_algebra_vector mom;
   suNg_field *u;
   suNg_av_field *f;
   double *s;
   double eps;
   double err, diff;
   int x[4];

	mon_lw_par par;
	par.c0 = c0;
	par.c1 = c1;
	par.beta = beta;

   u=alloc_gfield(&glattice);
   f=alloc_avfield(&glattice);
	_MASTER_FOR(&glattice,ix) {
		for(int mu=0; mu<4; mu++) {
			_algebra_vector_zero_g(*_4FIELD_AT(f,ix,mu));
		}
	}
	s=malloc(glattice.gsize_gauge*sizeof(double));

   random_u(u_gauge);
   start_gf_sendrecv(u_gauge);
   complete_gf_sendrecv(u_gauge);

   lw_force(1.,f,&par);

   calculate_stfld(NOCOMM);
   _MASTER_FOR(&glattice,ix) {
      s[ix]=lw_action_density(ix,beta,c0,c1);
   }

   memcpy(u->ptr,u_gauge->ptr,4*glattice.gsize_gauge*sizeof(suNg));

   eps=.01;
   for(int k=0;k<6;k++) {
      err=0.;
      for(x[0]=0;x[0]<GLB_T;x[0]++)
      for(x[1]=0;x[1]<GLB_X;x[1]++)
      for(x[2]=0;x[2]<GLB_Y;x[2]++)
      for(x[3]=0;x[3]<GLB_Z;x[3]++) {
         int local=(1==0);
         if((x[0] >= zerocoord[0] && x[0] < zerocoord[0]+T)
   		   && (x[1] >= zerocoord[1] && x[1] < zerocoord[1]+X)
   		   && (x[2] >= zerocoord[2] && x[2] < zerocoord[2]+Y)
   		   && (x[3] >= zerocoord[3] && x[3] < zerocoord[3]+Z)) local=(1==1);

         int ix=-1;
         if(local) ix=ipt(x[0]-zerocoord[0],x[1]-zerocoord[1],x[2]-zerocoord[2],x[3]-zerocoord[3]);

         for(int mu=0;mu<4;mu++) {
            if(local) {
               gauss((double*)(&mom),NG*NG-1);
               ExpX(eps,&mom,pu_gauge(ix,mu));
            }
            start_gf_sendrecv(u_gauge);
            complete_gf_sendrecv(u_gauge);

            double deltaS=0.;
            calculate_stfld(NOCOMM);
            _MASTER_FOR(&glattice,iy) {
              deltaS+=s[iy]-lw_action_density(iy,beta,c0,c1);
            }
            global_sum(&deltaS,1);

            double Xf=0.;
            if(local) {
               for(int i=0;i<NG*NG-1;i++) {
                  Xf += _FUND_NORM2 * mom.c[i] * _4FIELD_AT(f,ix,mu)->c[i];
               }
            }
            global_sum(&Xf,1);

            diff=fabs(Xf-deltaS/eps);
            if(diff>err) err=diff;
            // lprintf("TEST",0,"x=(%d,%d,%d,%d) mu=%d   X.force = %e   DeltaS/eps = %e\n",x[0],x[1],x[2],x[3],mu,Xf,deltaS/eps);


            memcpy(u_gauge->ptr,u->ptr,4*glattice.gsize_gauge*sizeof(suNg));
            start_gf_sendrecv(u_gauge);
            complete_gf_sendrecv(u_gauge);
         }
      }

      lprintf("TEST",0,"pars=(%f,%f,%f) :  Derivative of the action,  eps = %.3e     fabs(DeltaS - eps*X.force)/eps^2 = %e\n",beta,c0,c1,eps,err/eps);

      eps*=.1;
   }

   free_gfield(u);
   free(s);
}
*/
