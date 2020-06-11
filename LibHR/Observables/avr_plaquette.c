/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *
* All rights reserved.                                                      *
\***************************************************************************/

/*******************************************************************************
*
* File plaquette.c
*
* Routines for the average plaquette
*
*******************************************************************************/

#include "global.h"
#include "suN.h"
#include "utils.h"
#include "update.h"
#include "memory.h"
#include "random.h"
#include "dirac.h"
#include "representation.h"
#include "linear_algebra.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "logger.h"
#include "communications.h"

double plaq(int ix, int mu, int nu)
{
  int iy, iz;
  double p;
  suNg *v1, *v2, *v3, *v4, w1, w2, w3;

  iy = iup(ix, mu);
  iz = iup(ix, nu);

  v1 = pu_gauge(ix, mu);
  v2 = pu_gauge(iy, nu);
  v3 = pu_gauge(iz, mu);
  v4 = pu_gauge(ix, nu);

  _suNg_times_suNg(w1, (*v1), (*v2));
  _suNg_times_suNg(w2, (*v4), (*v3));
  _suNg_times_suNg_dagger(w3, w1, w2);

  _suNg_trace_re(p, w3);

#ifdef PLAQ_WEIGHTS
  if (plaq_weight == NULL)
    return p;
  return plaq_weight[ix * 16 + mu * 4 + nu] * p;
#else
  return p;
#endif
}

void cplaq(double complex *ret, int ix, int mu, int nu)
{
  int iy, iz;
  suNg *v1, *v2, *v3, *v4, w1, w2, w3;
  double tmpre = 0.;
  double tmpim = 0.;

  iy = iup(ix, mu);
  iz = iup(ix, nu);

  v1 = pu_gauge(ix, mu);
  v2 = pu_gauge(iy, nu);
  v3 = pu_gauge(iz, mu);
  v4 = pu_gauge(ix, nu);

  _suNg_times_suNg(w1, (*v1), (*v2));
  _suNg_times_suNg(w2, (*v4), (*v3));
  _suNg_times_suNg_dagger(w3, w1, w2);

  _suNg_trace_re(tmpre, w3);

#ifndef GAUGE_SON
  _suNg_trace_im(tmpim, w3);
#endif
  *ret = tmpre + I * tmpim;

#ifdef PLAQ_WEIGHTS
  if (plaq_weight != NULL)
  {
    *ret *= plaq_weight[ix * 16 + mu * 4 + nu];
  }
#endif
}

double avr_plaquette()
{
  static double pa = 0.;

  _OMP_PRAGMA(single)
  {
    pa = 0.;
  }

  _PIECE_FOR(&glattice, ixp)
  {
    if (ixp == glattice.inner_master_pieces)
    {
      _OMP_PRAGMA(master)
      /* wait for gauge field to be transfered */
      complete_gf_sendrecv(u_gauge);
      _OMP_PRAGMA(barrier)
    }
    _SITE_FOR_SUM(&glattice, ixp, ix, pa)
    {
      pa += plaq(ix, 1, 0);
      pa += plaq(ix, 2, 0);
      pa += plaq(ix, 2, 1);
      pa += plaq(ix, 3, 0);
      pa += plaq(ix, 3, 1);
      pa += plaq(ix, 3, 2);
    }
  }

  global_sum(&pa, 1);

#ifdef BC_T_OPEN
  pa /= 6.0 * NG * GLB_VOL3 * (GLB_T - 1);
#elif BASIC_SF
  pa /= 6.0 * NG * GLB_VOL3 * (GLB_T - 2);
#else
  pa /= 6.0 * NG * GLB_VOLUME;
#endif

  return pa;
}

void avr_plaquette_time(double *plaqt, double *plaqs)
{
  int ix;
  int tc;
  for (int nt = 0; nt < GLB_T; nt++)
    plaqt[nt] = plaqs[nt] = 0.0;

  for (int nt = 0; nt < T; nt++)
  {
    tc = (zerocoord[0] + nt + GLB_T) % GLB_T;
    for (int nx = 0; nx < X; nx++)
      for (int ny = 0; ny < Y; ny++)
        for (int nz = 0; nz < Z; nz++)
        {
          ix = ipt(nt, nx, ny, nz);
          plaqt[tc] += plaq(ix, 1, 0);
          plaqt[tc] += plaq(ix, 2, 0);
          plaqt[tc] += plaq(ix, 3, 0);
          plaqs[tc] += plaq(ix, 2, 1);
          plaqs[tc] += plaq(ix, 3, 1);
          plaqs[tc] += plaq(ix, 3, 2);
        }
    plaqt[tc] /= 3.0 * NG * GLB_VOL3;
    plaqs[tc] /= 3.0 * NG * GLB_VOL3;
  }

  for (int nt = 0; nt < GLB_T; nt++)
  {
    global_sum(&plaqt[nt], 1);
    global_sum(&plaqs[nt], 1);
  }
}

void full_plaquette()
{
  static double complex pa[6];
  static double complex r0;
  static double complex r1;
  static double complex r2;
  static double complex r3;
  static double complex r4;
  static double complex r5;

  _OMP_PRAGMA(single)
  {
    r0 = 0.;
    r1 = 0.;
    r2 = 0.;
    r3 = 0.;
    r4 = 0.;
    r5 = 0.;
  }

  _PIECE_FOR(&glattice, ixp)
  {
    if (ixp == glattice.inner_master_pieces)
    {
      _OMP_PRAGMA(master)
      /* wait for gauge field to be transfered */
      complete_gf_sendrecv(u_gauge);
      _OMP_PRAGMA(barrier)
    }

    _SITE_FOR_SUM(&glattice, ixp, ix, r0, r1, r2, r3, r4, r5)
    {
      double complex tmp;

      cplaq(&tmp, ix, 1, 0);
      r0 += tmp;

      cplaq(&tmp, ix, 2, 0);
      r1 += tmp;

      cplaq(&tmp, ix, 2, 1);
      r2 += tmp;

      cplaq(&tmp, ix, 3, 0);
      r3 += tmp;

      cplaq(&tmp, ix, 3, 1);
      r4 += tmp;

      cplaq(&tmp, ix, 3, 2);
      r5 += tmp;

      /*if(twbc_plaq[ix*16+2*4+1]==-1 &&
		    twbc_plaq[ix*16+3*4+1]==-1 &&
		    twbc_plaq[ix*16+3*4+2]==-1) {
		  cplaq(&tmp,ix,1,0);
		  lprintf("LOCPL",0,"Plaq( %d , %d , %d ) = ( %f , %f )\n",t,1,0,tmp.re,tmp.im);
		  cplaq(&tmp,ix,2,0);
		  lprintf("LOCPL",0,"Plaq( %d , %d , %d ) = ( %f , %f )\n",t,2,0,tmp.re,tmp.im);
		  cplaq(&tmp,ix,2,1);
		  lprintf("LOCPL",0,"Plaq( %d , %d , %d ) = ( %f , %f )\n",t,2,1,tmp.re,tmp.im);
		  cplaq(&tmp,ix,3,0);
		  lprintf("LOCPL",0,"Plaq( %d , %d , %d ) = ( %f , %f )\n",t,3,0,tmp.re,tmp.im);
		  cplaq(&tmp,ix,3,1);
		  lprintf("LOCPL",0,"Plaq( %d , %d , %d ) = ( %f , %f )\n",t,3,1,tmp.re,tmp.im);
		  cplaq(&tmp,ix,3,2);
		  lprintf("LOCPL",0,"Plaq( %d , %d , %d ) = ( %f , %f )\n",t,3,2,tmp.re,tmp.im);
		  t++;
		    } */
    }
  }

  _OMP_PRAGMA(single)
  {
    pa[0] = r0;
    pa[1] = r1;
    pa[2] = r2;
    pa[3] = r3;
    pa[4] = r4;
    pa[5] = r5;
  }

  global_sum((double *)pa, 12);
  _OMP_PRAGMA(single)
  for (int k = 0; k < 6; k++)
  {

#ifdef BC_T_OPEN
    pa[k] /= NG * GLB_VOLUME * (GLB_T - 1) / GLB_T;
#else
    pa[k] /= NG * GLB_VOLUME;
#endif
  }

  lprintf("PLAQ", 0, "Plaq(%d,%d) = ( %f , %f )\n", 1, 0, creal(pa[0]), cimag(pa[0]));
  lprintf("PLAQ", 0, "Plaq(%d,%d) = ( %f , %f )\n", 2, 0, creal(pa[1]), cimag(pa[1]));
  lprintf("PLAQ", 0, "Plaq(%d,%d) = ( %f , %f )\n", 2, 1, creal(pa[2]), cimag(pa[2]));
  lprintf("PLAQ", 0, "Plaq(%d,%d) = ( %f , %f )\n", 3, 0, creal(pa[3]), cimag(pa[3]));
  lprintf("PLAQ", 0, "Plaq(%d,%d) = ( %f , %f )\n", 3, 1, creal(pa[4]), cimag(pa[4]));
  lprintf("PLAQ", 0, "Plaq(%d,%d) = ( %f , %f )\n", 3, 2, creal(pa[5]), cimag(pa[5]));
}

/*void avr_ts_plaquette()
{
  int ix;
  static double pt = 0., ps = 0.;
  for (int nt = 0; nt < T; nt++)
  {
    for (int nx = 0; nx < X; nx++)
      for (int ny = 0; ny < Y; ny++)
        for (int nz = 0; nz < Z; nz++)
        {
          ix = ipt(nt, nx, ny, nz);
          pt += plaq(ix, 1, 0);
          pt += plaq(ix, 2, 0);
          ps += plaq(ix, 2, 1);
          pt += plaq(ix, 3, 0);
          ps += plaq(ix, 3, 1);
          ps += plaq(ix, 3, 2);
        }
    pt /= 3.0 * NG * GLB_VOL3;
    ps /= 3.0 * NG * GLB_VOL3;
    printf("%d %.10e %.10e\n",nt+zerocoord[0],pt,ps);

    pt=ps=0;
  }
}*/

double local_plaq(int ix)
{
  double pa;

  pa = plaq(ix, 1, 0);
  pa += plaq(ix, 2, 0);
  pa += plaq(ix, 2, 1);
  pa += plaq(ix, 3, 0);
  pa += plaq(ix, 3, 1);
  pa += plaq(ix, 3, 2);

  return pa;
}

void full_momenta(suNg_av_field *momenta)
{
  scalar_field *la = alloc_sfield(1, &glattice);
  _MASTER_FOR(&glattice, i)
  {
    double a = 0., tmp;
    /* Momenta */
    for (int j = 0; j < 4; ++j)
    {
      suNg_algebra_vector *cmom = momenta->ptr + coord_to_index(i, j);
      _algebra_vector_sqnorm_g(tmp, *cmom);
      a += tmp; /* this must be positive */
    }
    a *= 0.5 * _FUND_NORM2;
    *_FIELD_AT(la, i) = a;
  }
  static double mom;
  _OMP_PRAGMA(single)
  {
    mom = 0.;
  }
  _MASTER_FOR_SUM(la->type, i, mom)
  {
    mom += *_FIELD_AT(la, i);
  }
  lprintf("MOMENTA", 0, "%1.8g\n", mom);
  free_sfield(la);
}

void cplaq_wrk(double complex *ret, int ix, int mu, int nu)
{
  int iy, iz;
  suNg *v1, *v2, *v3, *v4, w1, w2, w3;
#ifdef GAUGE_SON
  double tmpre;
#endif
  iy = iup_wrk(ix, mu);
  iz = iup_wrk(ix, nu);

  v1 = pu_gauge_wrk(ix, mu);
  v2 = pu_gauge_wrk(iy, nu);
  v3 = pu_gauge_wrk(iz, mu);
  v4 = pu_gauge_wrk(ix, nu);

  _suNg_times_suNg(w1, (*v1), (*v2));
  _suNg_times_suNg(w2, (*v4), (*v3));
  _suNg_times_suNg_dagger(w3, w1, w2);

#ifndef GAUGE_SON
  _suNg_trace(*ret, w3);
#else
  _suNg_trace_re(tmpre, w3);
  *ret = tmpre;
#endif

#ifdef PLAQ_WEIGHTS
  if (plaq_weight != NULL)
  {
    *ret *= plaq_weight[ix * 16 + mu * 4 + nu];
  }
#endif
}

double complex avr_plaquette_wrk()
{
  static double complex pa, tmp;
  suNg_field *_u = u_gauge_wrk();
  start_gf_sendrecv(_u);

  _OMP_PRAGMA(single)
  {
    pa = tmp = 0.;
  }

  _PIECE_FOR(&glattice, ixp)
  {
    if (ixp == glattice.inner_master_pieces)
    {
      _OMP_PRAGMA(master)
      /* wait for gauge field to be transfered */
      complete_gf_sendrecv(_u);
      _OMP_PRAGMA(barrier)
    }
    _SITE_FOR_SUM(&glattice, ixp, ix, pa)
    {
      cplaq_wrk(&tmp, ix, 1, 0);
      pa += tmp;
      cplaq_wrk(&tmp, ix, 2, 0);
      pa += tmp;
      cplaq_wrk(&tmp, ix, 2, 1);
      pa += tmp;
      cplaq_wrk(&tmp, ix, 3, 0);
      pa += tmp;
      cplaq_wrk(&tmp, ix, 3, 1);
      pa += tmp;
      cplaq_wrk(&tmp, ix, 3, 2);
      pa += tmp;
    }
  }

  global_sum((double *)(&pa), 2);

#ifdef BC_T_OPEN
  pa /= 6.0 * NG * GLB_VOLUME * (GLB_T - 1) / GLB_T;
#else
  pa /= 6.0 * NG * GLB_VOLUME;
#endif

  return pa;
}
