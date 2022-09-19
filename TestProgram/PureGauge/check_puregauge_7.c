/*******************************************************************************
 *
 * Gauge and N-ality invariance of the torellons operators
 *
 *******************************************************************************/

#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "global.h"
#include "io.h"
#include "update.h"
#include "geometry.h"
#include "logger.h"
#include "random.h"
#include "memory.h"
#include "dirac.h"
#include "linear_algebra.h"
#include "representation.h"
#include "communications.h"
#include "setup.h"
#include "utils.h"
#include "observables.h"
#include "glueballs.h"
#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288419716939937510
#endif

static suNg_field *g;

static void random_g(void)
{
  _MASTER_FOR(&glattice, ix)
  {
    random_suNg(_FIELD_AT(g, ix));
  }

  start_gt_sendrecv(g);
  complete_gt_sendrecv(g);
}

static void transform_u(void)
{
  _MASTER_FOR(&glattice, ix)
  {
    suNg v;
    for (int mu = 0; mu < 4; mu++)
    {
      int iy = iup(ix, mu);
      suNg *u = pu_gauge(ix, mu);
      _suNg_times_suNg_dagger(v, *u, *_FIELD_AT(g, iy));
      _suNg_times_suNg(*u, *_FIELD_AT(g, ix), v);
    }
  }

  start_gf_sendrecv(u_gauge);
  represent_gauge_field();
}

static void n_ality_transform(int dir)
{
  int n[4];
  double complex centre = cexp(2. * M_PI * I / NG);
  suNg unew;
  suNg *uold;
  for (n[0] = 0; n[0] < T; n[0]++)
    for (n[1] = 0; n[1] < X; n[1]++)
      for (n[2] = 0; n[2] < Y; n[2]++)
        for (n[3] = 0; n[3] < Z; n[3]++)
        {
          if (n[dir] + zerocoord[dir] == 1)
          {
            int ix = ipt(n[0], n[1], n[2], n[3]);
            uold = pu_gauge(ix, dir);
            _suNg_mul(unew, centre, *uold);
            *uold = unew;
          }
        }
  start_gf_sendrecv(u_gauge);
  represent_gauge_field();
}

int main(int argc, char *argv[])
{

  int return_value = 0, n, nt;

  double complex *dop, *dop1;
  double max_diff[2];
  double min_size;
  setup_process(&argc, &argv);

  setup_gauge_fields();

  report_tor_group_setup();

  initialize_spatial_active_slices(NULL);
  cor_list corrs;

  corrs.n_entries = corrs.n_corrs = GLB_T;

  corrs.list = malloc(sizeof(cor_points) * GLB_T);
  for (n = 0; n < GLB_T / 2; n++)
  {
    corrs.list[n].t1 = GLB_T / 2 - 1;
    corrs.list[n].t2 = GLB_T / 2 + n;
    corrs.list[n].n_pairs = 1;
  }
  for (n = 0; n < GLB_T / 2; n++)
  {
    corrs.list[n + GLB_T / 2].t1 = 0;
    corrs.list[n + GLB_T / 2].t2 = GLB_T / 2 + n;
    corrs.list[n + GLB_T / 2].n_pairs = 1;
  }
  /* allocate additional memory */
  g = alloc_gtransf(&glattice);

  lprintf("MAIN", 0, "Generating a random gauge field... ");
  random_u(u_gauge);
  start_gf_sendrecv(u_gauge);
  represent_gauge_field();
  lprintf("MAIN", 0, "done.\n\n");
  double complex **polyf;

  dop = malloc(T * total_n_tor_op * sizeof(double complex));
  dop1 = malloc(T * total_n_tor_op * sizeof(double complex));
  polyf = malloc(sizeof(double complex *) * 3);
  polyf[0] = amalloc(sizeof(double complex) * Y * Z * T, ALIGN);
  polyf[1] = amalloc(sizeof(double complex) * X * Z * T, ALIGN);
  polyf[2] = amalloc(sizeof(double complex) * X * Y * T, ALIGN);
  memset(polyf[0], 0, sizeof(double complex) * Y * Z * T);
  memset(polyf[1], 0, sizeof(double complex) * X * Z * T);
  memset(polyf[2], 0, sizeof(double complex) * X * Y * T);

  for (n = 0; n < T * total_n_tor_op; n++)
  {
    dop[n] = 0.;
    dop1[n] = 0.;
  }

  for (nt = 0; nt < T; nt++)
    eval_all_torellon_ops(nt, dop + nt * total_n_tor_op, polyf);

  collect_1pt_torellon_functions(&corrs, dop, polyf);

  for (n = 0; n < T * total_n_tor_op; n++)
    dop[n] /= NG * GLB_VOLUME;

  lprintf("MAIN", 0, "Generating and applying a random gauge transf... ");
  random_g();
  transform_u();

  lprintf("MAIN", 0, "done.\n");

  memset(polyf[0], 0, sizeof(double complex) * Y * Z * T);
  memset(polyf[1], 0, sizeof(double complex) * X * Z * T);
  memset(polyf[2], 0, sizeof(double complex) * X * Y * T);

  for (nt = 0; nt < T; nt++)
    eval_all_torellon_ops(nt, dop1 + nt * total_n_tor_op, polyf);

  collect_1pt_torellon_functions(&corrs, dop1, polyf);

  for (n = 0; n < T * total_n_tor_op; n++)
    dop1[n] /= NG * GLB_VOLUME;

  for (n = 0; n < T * total_n_tor_op; n++)
    dop1[n] -= dop[n];

  lprintf("MAIN", 0, "Checking gauge invariance of the %d torellon operators on each timeslice.\n", total_n_tor_op);

  max_diff[0] = max_diff[1] = -1.0;
  min_size = 10;

  for (n = 0; n < T * total_n_tor_op; n++)
  {
    if (fabs(creal(dop1[n])) > max_diff[0])
      max_diff[0] = fabs(creal(dop1[n]));

    if (fabs(cimag(dop1[n])) > max_diff[1])
      max_diff[1] = fabs(cimag(dop1[n]));

    if (fabs(creal(dop1[n])) > 1.e-12)
      return_value++;
    if (fabs(cimag(dop1[n])) > 1.e-12)
      return_value++;

    if (sqrt(creal(dop[n]) * creal(dop[n]) + cimag(dop[n]) * cimag(dop[n])) < min_size)
      min_size = sqrt(creal(dop[n]) * creal(dop[n]) + cimag(dop[n]) * cimag(dop[n]));

    if (sqrt(creal(dop[n]) * creal(dop[n]) + cimag(dop[n]) * cimag(dop[n])) < 10.e-10)
    {
      lprintf("MAIN", 0, "Operator %d on timeslice %d seems to be numerically zero\n", n % total_n_tor_op, n / total_n_tor_op);
      return_value++;
    }
  }
  global_sum_int(&return_value, 1);
  global_max(max_diff, 2);
  global_min(&min_size, 1);

  lprintf("MAIN", 0, "Maximal normalized real difference = %.16e\n", max_diff[0]);
  lprintf("MAIN", 0, "(should be around 1*10^(-15) or so)\n");
  lprintf("MAIN", 0, "Maximal normalized imaginary difference = %.16e\n", max_diff[1]);
  lprintf("MAIN", 0, "(should be around 1*10^(-15) or so)\n");
  lprintf("MAIN", 0, "Minimal amplitude of the operators = %.4e\n", min_size);
  lprintf("MAIN", 0, "(should be greater 1*10^(-15) or so)\n\n");

  free_gtransf(g);
  int inx = ipt(3, 3, 3, 3);

  wilson_lines *pol = polyleg(inx, 1);
  double complex cb = pol->tr;

  lprintf("MAIN", 0, "Pa: single polyakov line:  %.10e +I*(%.10e).\n", creal(cb), cimag(cb));
  lprintf("MAIN", 0, "Applying a N-ality transf to a X hyperplane... ");
  n_ality_transform(1);
  lprintf("MAIN", 0, "done.\n");

  inx = ipt(3, 2, 3, 3);
  pol = polyleg(inx, 1);
  double complex ca = pol->tr;

  lprintf("MAIN", 0, "\nPb: single polyakov line after the transf. :   %.10e +I*(%.10e).\n", creal(ca), cimag(ca));

  ca = (ca / cb) / cexp(2. * M_PI * I / NG) - 1.0;

  lprintf("MAIN", 0, "Checking |(Pa/Pb)/cexp(2. * M_PI * I / NG) - 1|= %.10e\n", sqrt(creal(ca) * creal(ca) + cimag(ca) * cimag(ca)));
  lprintf("MAIN", 0, "(should be around 1*10^(-15) or so)\n\n");
  if (sqrt(creal(ca) * creal(ca) + cimag(ca) * cimag(ca)) > 1.e-12)
    return_value++;

  for (n = 0; n < T * total_n_tor_op; n++)
    dop1[n] = 0.;

  for (nt = 0; nt < T; nt++)
    eval_all_torellon_ops(nt, dop1 + nt * total_n_tor_op, polyf);

  for (n = 0; n < T * total_n_tor_op; n++)
    dop1[n] /= NG * GLB_VOLUME;

  for (n = 0; n < T * total_n_tor_op; n++)
    dop1[n] -= dop[n];

  lprintf("MAIN", 0, "Checking N-ality invariance of the %d torellon operators on each timeslice.\n", total_n_tor_op);

  for (n = 0; n < T * total_n_tor_op; n++)
  {

    if (fabs(creal(dop1[n])) > max_diff[0])
      max_diff[0] = fabs(creal(dop1[n]));

    if (fabs(cimag(dop1[n])) > max_diff[1])
      max_diff[1] = fabs(cimag(dop1[n]));

    if (fabs(creal(dop1[n])) > 1.e-12)
      return_value++;
    if (fabs(cimag(dop1[n])) > 1.e-12)
      return_value++;
  }
  global_max(max_diff, 2);

  lprintf("MAIN", 0, "Maximal normalized real difference = %.16e\n", max_diff[0]);
  lprintf("MAIN", 0, "(should be around 1*10^(-15) or so)\n");
  lprintf("MAIN", 0, "Maximal normalized imaginary difference = %.16e\n", max_diff[1]);
  lprintf("MAIN", 0, "(should be around 1*10^(-15) or so)\n\n");

  lprintf("MAIN", 0, "Applying a N-ality transf to a Y hyperplane... ");
  n_ality_transform(2);
  lprintf("MAIN", 0, "done.\n");
  for (n = 0; n < T * total_n_tor_op; n++)
    dop1[n] = 0.;

  for (nt = 0; nt < T; nt++)
    eval_all_torellon_ops(nt, dop1 + nt * total_n_tor_op, polyf);

  for (n = 0; n < T * total_n_tor_op; n++)
    dop1[n] /= NG * GLB_VOLUME;

  for (n = 0; n < T * total_n_tor_op; n++)
    dop1[n] -= dop[n];

  lprintf("MAIN", 0, "Checking N-ality invariance of the %d torellon operators on each timeslice.\n", total_n_tor_op);

  for (n = 0; n < T * total_n_tor_op; n++)
  {

    if (fabs(creal(dop1[n])) > max_diff[0])
      max_diff[0] = fabs(creal(dop1[n]));

    if (fabs(cimag(dop1[n])) > max_diff[1])
      max_diff[1] = fabs(cimag(dop1[n]));

    if (fabs(creal(dop1[n])) > 1.e-12)
      return_value++;
    if (fabs(cimag(dop1[n])) > 1.e-12)
      return_value++;
  }
  global_max(max_diff, 2);

  lprintf("MAIN", 0, "Maximal normalized real difference = %.16e\n", max_diff[0]);
  lprintf("MAIN", 0, "(should be around 1*10^(-15) or so)\n");
  lprintf("MAIN", 0, "Maximal normalized imaginary difference = %.16e\n", max_diff[1]);
  lprintf("MAIN", 0, "(should be around 1*10^(-15) or so)\n\n");

  lprintf("MAIN", 0, "Applying a N-ality transf to a Z hyperplane... ");
  n_ality_transform(3);
  lprintf("MAIN", 0, "done.\n");
  for (n = 0; n < T * total_n_tor_op; n++)
    dop1[n] = 0.;

  for (nt = 0; nt < T; nt++)
    eval_all_torellon_ops(nt, dop1 + nt * total_n_tor_op, polyf);

  for (n = 0; n < T * total_n_tor_op; n++)
    dop1[n] /= NG * GLB_VOLUME;

  for (n = 0; n < T * total_n_tor_op; n++)
    dop1[n] -= dop[n];

  lprintf("MAIN", 0, "Checking N-ality invariance of the %d torellon operators on each timeslice.\n", total_n_tor_op);

  for (n = 0; n < T * total_n_tor_op; n++)
  {

    if (fabs(creal(dop1[n])) > max_diff[0])
      max_diff[0] = fabs(creal(dop1[n]));

    if (fabs(cimag(dop1[n])) > max_diff[1])
      max_diff[1] = fabs(cimag(dop1[n]));

    if (fabs(creal(dop1[n])) > 1.e-12)
      return_value++;
    if (fabs(cimag(dop1[n])) > 1.e-12)
      return_value++;
  }
  global_max(max_diff, 2);

  lprintf("MAIN", 0, "Maximal normalized real difference = %.16e\n", max_diff[0]);
  lprintf("MAIN", 0, "(should be around 1*10^(-15) or so)\n");
  lprintf("MAIN", 0, "Maximal normalized imaginary difference = %.16e\n", max_diff[1]);
  lprintf("MAIN", 0, "(should be around 1*10^(-15) or so)\n");

  global_sum_int(&return_value, 1);

  // write_gauge_field("myconf.dat");

  finalize_process();
  return return_value;
}
