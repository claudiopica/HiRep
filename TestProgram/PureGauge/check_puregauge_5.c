/*******************************************************************************
*
* Gauge invariance of the glueball operators
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

int main(int argc, char *argv[])
{

  int return_value = 0, n, nt;

  double complex *dop, *dop1;

  setup_process(&argc, &argv);

  setup_gauge_fields();

  /* allocate additional memory */
  g = alloc_gtransf(&glattice);

  lprintf("MAIN", 0, "Generating a random gauge field... ");
  random_u(u_gauge);
  start_gf_sendrecv(u_gauge);
  represent_gauge_field();
  lprintf("MAIN", 0, "done.\n\n");

  dop = malloc(total_n_glue_op * sizeof(double complex));
  dop1 = malloc(total_n_glue_op * sizeof(double complex));

  for (n = 0; n < total_n_glue_op; n++)
  {
    dop[n] = 0.;
    dop1[n] = 0.;
  }

  for (nt = 0; nt < T; nt++)
    eval_all_glueball_ops(nt, dop);

  global_sum((double *)dop, 2 * total_n_glue_op);
  for (n = 0; n < total_n_glue_op; n++)
    dop[n] /= NG * GLB_VOLUME;

  lprintf("MAIN", 0, "Generating and applying a random gauge transf... ");
  random_g();
  transform_u();

  lprintf("MAIN", 0, "done.\n");

  for (nt = 0; nt < T; nt++)
    eval_all_glueball_ops(nt, dop1);

  global_sum((double *)dop1, 2 * total_n_glue_op);
  for (n = 0; n < total_n_glue_op; n++)
    dop1[n] /= NG * GLB_VOLUME;

  for (n = 0; n < total_n_glue_op; n++)
    dop[n] -= dop1[n];

  for (n = 0; n < total_n_glue_op; n++)
  {
    lprintf("MAIN", 0, "Checking gauge invariance of the operator %d.\n ", n);
    lprintf("MAIN", 0, "Maximal normalized real difference = %.4e\n", fabs(creal(dop[n])));
    lprintf("MAIN", 0, "(should be around 1*10^(-15) or so)\n");
    lprintf("MAIN", 0, "Maximal normalized imaginary difference = %.4e\n", fabs(cimag(dop[n])));
    lprintf("MAIN", 0, "(should be around 1*10^(-15) or so)\n\n");
    if (fabs(creal(dop[n])) > 10.e-14)
      return_value++;
    if (fabs(cimag(dop[n])) > 10.e-14)
      return_value++;
  }

  free_gtransf(g);

  finalize_process();
  return return_value;
}
