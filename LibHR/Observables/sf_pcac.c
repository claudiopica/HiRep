#include "global.h"
#include "communications.h"
#include "linear_algebra.h"
#include "inverters.h"
#include "suN.h"
#include "observables.h"
#include "dirac.h"
#include "utils.h"
#include "memory.h"
#include "error.h"
#include "geometry.h"
#include "spinor_field.h"
#include "gamma_spinor.h"
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "logger.h"
#include <stdio.h>
#include "data_storage.h"

#ifdef ROTATED_SF
#error This code has never been tested
#endif

#ifdef ROTATED_SF
#include "update.h"

extern rhmc_par _update_par; /* Update/update_rhmc.c */

#endif /* ROTATED_SF */

#ifdef BASIC_SF
static double hmass;
static void H_sf(spinor_field *out, spinor_field *in)
{
#if defined(WITH_CLOVER) || defined(WITH_EXPCLOVER)
  g5Cphi(hmass, out, in);
#else
  g5Dphi(hmass, out, in);
#endif
}
#endif
#ifdef ROTATED_SF
static double hmass;
static void H2_sf(spinor_field *out, spinor_field *in)
{
  g5Dphi_sq(hmass, out, in);
}
#endif

/*
  this computes Hu^{-1}
*/
int SF_quark_propagator(spinor_field *in, double mass, spinor_field *out, double acc)
{
  int cgiter = 0;
#ifdef BASIC_SF

  static MINRES_par MINRESpar;
  hmass = mass;

  MINRESpar.err2 = acc;
  MINRESpar.max_iter = 0;
  cgiter = MINRES(&MINRESpar, &H_sf, in, out, 0);
#elif defined(ROTATED_SF)

  static mshift_par inv_par;
  static spinor_field *chi = NULL;
  if (chi == NULL)
    chi = alloc_spinor_field_f(1, &glattice);
  hmass = mass;

  inv_par.n = 1;
  inv_par.shift = malloc(sizeof(double));
  inv_par.shift[0] = 0.;
  inv_par.err2 = acc;
  inv_par.max_iter = 0;
  /*this change of sign is equivalent to evaluate H^\dagger*/
  _update_par.SF_sign = -_update_par.SF_sign;
  g5Dphi(mass, chi, in);
  _update_par.SF_sign = -_update_par.SF_sign;

  cgiter = cg_mshift(&inv_par, &H2_sf, chi, out);
  free(inv_par.shift);

#endif
  return cgiter;
}

/*
  f1 can be computed inserting the sources in the correlator in the lower or higher border. The two procedures give identical results. We keep only the computation with sources on the lower border.
*/

data_storage_array *SF_PCAC_wall_corr(double mass, double acc, storage_switch swh)
{
  data_storage_array *ret = NULL;

#ifdef BASIC_SF

  int i, j, ix0, ix1, ix2, ix3;
  double f_P[GLB_T], f_A[GLB_T], g_P[GLB_T], g_A[GLB_T], f_1 = 0, temp;
  spinor_field *prop;
  spinor_field *source;
  suNf_spinor *sptr;
  suNf_spinor *stmp;
  stmp = malloc(2 * sizeof(suNf_spinor));
  suNf_spinor *sbord;
  sbord = malloc(4 * NF * sizeof(suNf_spinor));
  suNf *uptr;
  prop = alloc_spinor_field_f(4 * NF, &glattice);
  source = alloc_spinor_field_f(4 * NF, &glattice);

  if (swh == STORE)
  {
    ret = allocate_data_storage_array(3);
    int idx[2] = {2, GLB_T};
    allocate_data_storage_element(ret, 0, 2, idx);
    allocate_data_storage_element(ret, 1, 2, idx);
    idx[0] = 1;
    allocate_data_storage_element(ret, 2, 1, idx);
  }
  for (ix0 = 0; ix0 < GLB_T; ix0++)
  {
    f_A[ix0] = 0;
    g_A[ix0] = 0;
    f_P[ix0] = 0;
    g_P[ix0] = 0;
  }

  /*U' and P+ on source (actually P- since there is a g5 that needs to be commuted through)*/
  for (int s = 0; s < 4 * NF; s++)
  {
    spinor_field_zero_f(&source[s]);
    if (COORD[0] == 0)
    {
      _spinor_zero_f(stmp[0]);
      stmp[0].c[s % 4].c[s / 4] = 1.;
      _spinor_g5_assign_f(stmp[0]);

      for (ix1 = 0; ix1 < X; ix1++)
        for (ix2 = 0; ix2 < Y; ix2++)
          for (ix3 = 0; ix3 < Z; ix3++)
          {
            i = ipt(2, ix1, ix2, ix3);
            uptr = pu_gauge_f(idn(i, 0), 0);
            for (j = 0; j < 4; j++)
            {
              _suNf_inverse_multiply(stmp[1].c[j], *uptr, stmp[0].c[j]);
            }
            sptr = _FIELD_AT(&source[s], i);
            _spinor_pminus_f(*sptr, stmp[1]);
          }
    }
  }

  /*get propagator to all points*/
  for (int s = 0; s < 4 * NF; s++)
  {
    spinor_field_zero_f(&prop[s]);
    SF_quark_propagator(&source[s], mass, &prop[s], acc);
  }

  /*get time averaged correlators for each timeslice*/
  /*f_P*/
  for (int s = 0; s < 4 * NF; s++)
  {
    for (ix0 = 0; ix0 < T; ix0++)
      for (ix1 = 0; ix1 < X; ix1++)
        for (ix2 = 0; ix2 < Y; ix2++)
          for (ix3 = 0; ix3 < Z; ix3++)
          {
            i = ipt(ix0, ix1, ix2, ix3);
            sptr = _FIELD_AT(&prop[s], i);
            /*f_P*/
            _spinor_prod_re_f(temp, *sptr, *sptr);
            f_P[(zerocoord[0] + ix0 - 1 + GLB_T) % GLB_T] += temp;
            /*f_A*/
            /*gamma_0*/
            _spinor_g0_f(stmp[0], (*sptr));
            _spinor_prod_re_f(temp, *sptr, stmp[0]);
            f_A[(zerocoord[0] + ix0 - 1 + GLB_T) % GLB_T] += temp;
          }
  }

  /* f_P = prop^dag prop */
  /* f_A = prop^dag gamma0 prop */

  global_sum((double *)f_P, GLB_T);

  for (ix0 = 0; ix0 < GLB_T - 1; ix0++)
  {
    f_P[ix0] /= 2.0 * (double)GLB_VOL3;
    lprintf("PC_wall_AC", 10, "f_P( %d ) = %.16e \n", ix0, f_P[ix0]);
  }
  global_sum((double *)f_A, GLB_T);
  for (ix0 = 0; ix0 < GLB_T - 1; ix0++)
  {
    f_A[ix0] /= -2.0 * (double)GLB_VOL3;
    lprintf("PC_wall_AC", 10, "f_A( %d )= %.16e\n", ix0, f_A[ix0]);
  }

  if (swh == STORE)
  {
    int idx[2];
    for (ix0 = 0; ix0 < GLB_T - 1; ix0++)
    {
      idx[1] = ix0;
      idx[0] = 0;
      *data_storage_element(ret, 0, idx) = f_P[ix0];
      idx[0] = 1;
      *data_storage_element(ret, 0, idx) = f_A[ix0];
    }
  }

  /*f_1 - NEED TO DO EACH color/dirac component separately, then combine at the end*/
  /*U' and P- on prop at T-2 (actually P+ since there is a g5 that needs to be commuted through)*/

  for (int s = 0; s < 4 * NF; s++)
  {
    _spinor_zero_f(sbord[s]);
    if (COORD[0] == NP_T - 1)
    {
      for (ix1 = 0; ix1 < X; ix1++)
        for (ix2 = 0; ix2 < Y; ix2++)
          for (ix3 = 0; ix3 < Z; ix3++)
          {
            i = ipt(T - 2, ix1, ix2, ix3);
            uptr = pu_gauge_f(i, 0);
            sptr = _FIELD_AT(&prop[s], i);
            for (j = 0; j < 4; j++)
            {
              _suNf_inverse_multiply(stmp[0].c[j], *uptr, sptr->c[j]);
            }
            _spinor_pplus_f(stmp[1], stmp[0]);
            _spinor_add_assign_f(sbord[s], stmp[1]);
          }
    }
  }

  global_sum((double *)sbord, sizeof(suNf_spinor) / sizeof(double) * 4 * NF);

  f_1 = 0;
  for (int s = 0; s < 4 * NF; s++)
  {
    _spinor_prod_re_f(temp, sbord[s], sbord[s]);
    f_1 += temp;
  }

  f_1 /= 2.0 * ((double)GLB_VOL3) * ((double)GLB_VOL3);

  lprintf("PC_wall_AC", 0, "f1 = %.16e\n", f_1);
  if (swh == STORE)
  {
    int idx = 0;
    *data_storage_element(ret, 2, &idx) = f_1;
  }

  /*Create wall source with g5 factor at t=T-2*/
  /*U and P- on source (again actually use P+ to account for commuting with g5 in source)*/
  for (int s = 0; s < 4 * NF; s++)
  {
    spinor_field_zero_f(&source[s]);
    if (COORD[0] == NP_T - 1)
    {
      _spinor_zero_f(stmp[0]);
      stmp[0].c[s % 4].c[s / 4] = 1.;
      _spinor_g5_assign_f(stmp[0]);

      for (ix1 = 0; ix1 < X; ix1++)
        for (ix2 = 0; ix2 < Y; ix2++)
          for (ix3 = 0; ix3 < Z; ix3++)
          {
            i = ipt(T - 2, ix1, ix2, ix3);
            uptr = pu_gauge_f(i, 0);
            for (j = 0; j < 4; j++)
            {
              _suNf_multiply(stmp[1].c[j], *uptr, stmp[0].c[j]);
            }
            sptr = _FIELD_AT(&source[s], i);
            _spinor_pplus_f(*sptr, stmp[1]);
          }
    }
  }

  /*get propagator to all points*/
  for (int s = 0; s < 4 * NF; s++)
  {
    spinor_field_zero_f(&prop[s]);
    SF_quark_propagator(&source[s], mass, &prop[s], acc);
  }

  /*get time averaged correlators for each timeslice (going back from T in time)*/
  for (int s = 0; s < 4 * NF; s++)
  {
    for (ix1 = 0; ix1 < X; ix1++)
      for (ix2 = 0; ix2 < Y; ix2++)
        for (ix3 = 0; ix3 < Z; ix3++)
          for (ix0 = 0; ix0 < T; ix0++)
          {
            i = ipt(ix0, ix1, ix2, ix3);
            sptr = _FIELD_AT(&prop[s], i);
            /*f_P*/
            _spinor_prod_re_f(temp, *sptr, *sptr);
            g_P[((GLB_T - 1) - (zerocoord[0] + ix0)) % GLB_T] += temp;
            /*f_A*/
            /*gamma_0*/
            _spinor_g0_f(stmp[0], (*sptr));
            _spinor_prod_re_f(temp, *sptr, stmp[0]);
            g_A[((GLB_T - 1) - (zerocoord[0] + ix0)) % GLB_T] += temp;
          }
  }

  global_sum((double *)g_P, GLB_T);
  for (ix0 = 0; ix0 < GLB_T - 1; ix0++)
  {
    g_P[ix0] /= 2.0 * (double)GLB_VOL3;
    lprintf("PC_wall_AC", 10, "g_P( %d )= %.16e \n", ix0, g_P[ix0]);
  }
  global_sum((double *)g_A, GLB_T);
  for (ix0 = 0; ix0 < GLB_T - 1; ix0++)
  {
    g_A[ix0] /= 2.0 * (double)GLB_VOL3;
    lprintf("PC_wall_AC", 10, "g_A( %d )= %.16e \n", ix0, g_A[ix0]);
  }

  if (swh == STORE)
  {
    int idx[2];
    for (ix0 = 0; ix0 < GLB_T - 1; ix0++)
    {
      idx[1] = ix0;
      idx[0] = 0;
      *data_storage_element(ret, 1, idx) = g_P[ix0];
      idx[0] = 1;
      *data_storage_element(ret, 1, idx) = g_A[ix0];
    }
  }

  free_spinor_field_f(source);
  free_spinor_field_f(prop);

  free(stmp);
  free(sbord);

  /*   return (double)(f_A[(int)(GLB_T/2)] - f_A[(int)(GLB_T/2)-2])/(4*f_P[(int)((GLB_T/2)-1)]); */

#elif defined(ROTATED_SF)

  static spinor_field *prop_uu = NULL;
  static spinor_field *source = NULL;
  static spinor_field *prop_dd = NULL;
  suNf_spinor chi[4 * NF + 1];
  suNf *uptr;
  int i, j, ix1, ix2, ix3, s1;

  if (source == NULL)
  {
    source = alloc_spinor_field_f(4 * NF, &glattice);
    prop_uu = alloc_spinor_field_f(4 * NF, &glattice);
    prop_dd = alloc_spinor_field_f(4 * NF, &glattice);
  }

  static chisf_mem *corr_mem = NULL;

  if (corr_mem == NULL)
    corr_mem = init_rotated_corr_mem();

  /****** We construct the 2 propagators

	  prop[s] = Huu^{-1} source[s]
	  prop[s] = Hdd^{-1} source[s]

	  once and for all *************/

  /******** Construction of the 2 sources**************/

  for (s1 = 0; s1 < 4 * NF; s1++)
  {
    spinor_field_zero_f(&source[s1]);
    _spinor_zero_f(chi[s1]);
    chi[s1].c[s1 % 4].c[s1 / 4] = 1.;
    if (COORD[0] == 0)
    {

      for (ix1 = 0; ix1 < X; ix1++)
        for (ix2 = 0; ix2 < Y; ix2++)
          for (ix3 = 0; ix3 < Z; ix3++)
          {
            i = ipt(2, ix1, ix2, ix3);
            uptr = pu_gauge_f(idn(i, 0), 0);

            for (j = 0; j < 4; j++)
            {
              _suNf_inverse_multiply(_FIELD_AT(&source[s1], i)->c[j], *uptr, chi[s1].c[j]);
            }
          }
    }
  }

  /****** Construction of the 2 propagators **********/
  /*
    prop[s] = Huu^{-1} source[s]
    prop[s] = Hdd^{-1} source[s]

  */
  for (s1 = 0; s1 < 4 * NF; s1++)
  {
    spinor_field_zero_f(&prop_uu[s1]);
    spinor_field_zero_f(&prop_dd[s1]);

    SF_quark_propagator(&source[s1], mass, &prop_uu[s1], acc);

    _update_par.SF_sign = -_update_par.SF_sign;
    SF_quark_propagator(&source[s1], mass, &prop_dd[s1], acc);
    _update_par.SF_sign = -_update_par.SF_sign;
  }

  rotated_gXuup(corr_mem, chi, prop_uu, prop_dd);

  rotated_gXddp(corr_mem, chi, prop_uu, prop_dd);

  rotated_gXudp(corr_mem, chi, prop_uu, prop_dd);

  rotated_gXdup(corr_mem, chi, prop_uu, prop_dd);

  rotated_gvtuup(corr_mem, chi, prop_uu, prop_dd);

  rotated_gvtddp(corr_mem, chi, prop_uu, prop_dd);

  rotated_gvtdup(corr_mem, chi, prop_uu, prop_dd);

  rotated_gvtudp(corr_mem, chi, prop_uu, prop_dd);

  rotated_g1uup(corr_mem, chi, prop_uu, prop_dd);

  rotated_g1ddp(corr_mem, chi, prop_uu, prop_dd);

  rotated_g1dup(corr_mem, chi, prop_uu, prop_dd);

  rotated_g1udp(corr_mem, chi, prop_uu, prop_dd);

  rotated_lXuup(corr_mem, chi, prop_uu, prop_dd);

  rotated_lXddp(corr_mem, chi, prop_uu, prop_dd);

  rotated_lXudp(corr_mem, chi, prop_uu, prop_dd);

  rotated_lXdup(corr_mem, chi, prop_uu, prop_dd);

  rotated_gXuum(corr_mem, chi, prop_uu, prop_dd);

  rotated_gXddm(corr_mem, chi, prop_uu, prop_dd);

  rotated_gXudm(corr_mem, chi, prop_uu, prop_dd);

  rotated_gXdum(corr_mem, chi, prop_uu, prop_dd);

  rotated_gvtuum(corr_mem, chi, prop_uu, prop_dd);

  rotated_gvtddm(corr_mem, chi, prop_uu, prop_dd);

  rotated_gvtdum(corr_mem, chi, prop_uu, prop_dd);

  rotated_gvtudm(corr_mem, chi, prop_uu, prop_dd);

  rotated_g1uum(corr_mem, chi, prop_uu, prop_dd);

  rotated_g1ddm(corr_mem, chi, prop_uu, prop_dd);

  rotated_g1dum(corr_mem, chi, prop_uu, prop_dd);

  rotated_g1udm(corr_mem, chi, prop_uu, prop_dd);

  rotated_lXuum(corr_mem, chi, prop_uu, prop_dd);

  rotated_lXddm(corr_mem, chi, prop_uu, prop_dd);

  rotated_lXudm(corr_mem, chi, prop_uu, prop_dd);

  rotated_lXdum(corr_mem, chi, prop_uu, prop_dd);

#endif
  return ret;
}
