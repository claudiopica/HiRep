/*******************************************************************************
*
* Gauge invariance of the spatial smearing operators
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
#define blk_level 7

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

double doubleop(int in, int mu, int nu)
{
  suNg *w1, *w2;
  suNg res, res1;
  int pit;
  double p;

  w1 = pu_gauge(in, mu);

  pit = iup(in, mu);
  w2 = pu_gauge(pit, mu);
  _suNg_times_suNg(res, *w1, *w2);

  pit = iup(pit, mu);
  w1 = pu_gauge(pit, nu);
  _suNg_times_suNg(res1, res, *w1);

  pit = iup(pit, nu);
  w1 = pu_gauge(pit, nu);
  _suNg_times_suNg(res, res1, *w1);

  pit = iup(pit, nu);
  pit = idn(pit, mu);
  w1 = pu_gauge(pit, mu);
  _suNg_times_suNg_dagger(res1, res, *w1);

  pit = idn(pit, mu);
  w1 = pu_gauge(pit, mu);
  _suNg_times_suNg_dagger(res, res1, *w1);

  pit = idn(pit, nu);
  w1 = pu_gauge(pit, nu);
  _suNg_times_suNg_dagger(res1, res, *w1);

  pit = idn(pit, nu);
  w1 = pu_gauge(pit, nu);
  _suNg_times_suNg_dagger(res, res1, *w1);

  _suNg_trace_re(p, res);
  return p;
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
  // smear_gauge_field();
}

int main(int argc, char *argv[])
{

  int return_value = 0;

  double complex splaq, splaq2;

  logger_map("DEBUG", "debug");

  setup_process(&argc, &argv);

  setup_gauge_fields();

  /* allocate additional memory */
  g = alloc_gtransf(&glattice);

  lprintf("MAIN", 0, "Generating a random gauge field... ");
  random_u(u_gauge);
  start_gf_sendrecv(u_gauge);
  represent_gauge_field();
  lprintf("MAIN", 0, "done.\n");
  lprintf("MAIN", 0, "\n\n");

  initialize_spatial_blocking(NULL);

  double dop = 0.;
  int id, ix, iy, iz, it;

  for (ix = 0; ix < X; ix++)
    for (iy = 0; iy < Y; iy++)
      for (iz = 0; iz < Z; iz++)
        for (it = 0; it < T; it++)
        {
          id = ipt(it, ix, iy, iz);
          dop += doubleop(id, 1, 2);
          dop += doubleop(id, 1, 3);
          dop += doubleop(id, 2, 3);
        }

  global_sum((double *)&dop, 2);
  dop /= NG * GLB_VOLUME;

  spatial_blocking(1, 1);

  splaq = spatial_plaquette_sblk();

  dop -= creal(splaq);
  dop = sqrt(dop * dop);


  lprintf("MAIN", 0, "Comparison of the double plaquette operator vs single plaquette on level of blocking 1.\n ");
  lprintf("MAIN", 0, "Maximal normalized difference = %.4e\n", dop);
  lprintf("MAIN", 0, "(should be around 1*10^(-15) or so)\n");
  if (dop > 10.e-14)
    return_value++;

  lprintf("MAIN", 0, "\n\n");
  spatial_blocking(blk_level, 1);

  splaq = spatial_plaquette_sblk();

  lprintf("MAIN", 0, "Generating and applying a random gauge transf... ");
  random_g();
  transform_u();

  lprintf("MAIN", 0, "done.\n");

  spatial_blocking(blk_level, 1);

  splaq2 = splaq - spatial_plaquette_sblk();

  dop = sqrt(_complex_prod_re(splaq2, splaq2));

  lprintf("MAIN", 0, "Checking gauge invariance of the spatial plaq on blocking level %d.\n ", blk_level);
  lprintf("MAIN", 0, "Maximal normalized difference = %.4e\n", dop);
  lprintf("MAIN", 0, "(should be around 1*10^(-15) or so)\n");
  if (dop > 10.e-14)
    return_value++;

  lprintf("MAIN", 0, "\n\n");

  spatial_blocking(blk_level - 2, 1);
  spatial_blocking(blk_level - 1, 0);
  spatial_blocking(blk_level, 0);

  splaq2 = splaq - spatial_plaquette_sblk();

  dop = sqrt(_complex_prod_re(splaq2, splaq2));

  lprintf("MAIN", 0, "Checking gauge invariance of the spatial plaq on blocking level %d (composing different blocking levels).\n ", blk_level);
  lprintf("MAIN", 0, "Maximal normalized difference = %.4e\n", dop);
  lprintf("MAIN", 0, "(should be around 1*10^(-15) or so)\n");
  if (dop > 10.e-14)
    return_value++;

  lprintf("MAIN", 0, "\n\n");

  spatial_blocking(blk_level + 2, 0);
  spatial_blocking(blk_level , 0);

  splaq2 = splaq - spatial_plaquette_sblk();

  dop = sqrt(_complex_prod_re(splaq2, splaq2));

  lprintf("MAIN", 0, "Checking gauge invariance of the spatial plaq on blocking level %d (without forcing the re-evaluation).\n ", blk_level );
  lprintf("MAIN", 0, "Maximal normalized difference = %.4e\n", dop);
  lprintf("MAIN", 0, "(should be around 1*10^(-15) or so)\n");
  if (dop > 10.e-14)
    return_value++;


lprintf("MAIN", 0, "\n\n");

  spatial_blocking(blk_level , 0);

  splaq2 = splaq - spatial_plaquette_sblk();

  dop = sqrt(_complex_prod_re(splaq2, splaq2));

  lprintf("MAIN", 0, "Checking gauge invariance of the spatial plaq on blocking level %d (on a level already evaluated).\n ", blk_level );
  lprintf("MAIN", 0, "Maximal normalized difference = %.4e\n", dop);
  lprintf("MAIN", 0, "(should be around 1*10^(-15) or so)\n");
  if (dop > 10.e-14)
    return_value++;





  free_gtransf(g);

  free_spatial_blocking();

  finalize_process();
  return return_value;

}
