/******************************************************************************
*
* Test of modules
*
******************************************************************************/

#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "io.h"
#include "random.h"
#include "error.h"
#include "geometry.h"
#include "memory.h"
#include "statistics.h"
#include "update.h"
#include "global.h"
#include "observables.h"
#include "suN.h"
#include "suN_types.h"
#include "dirac.h"
#include "linear_algebra.h"
#include "inverters.h"
#include "representation.h"
#include "utils.h"
#include "logger.h"
#include "communications.h"
#include "setup.h"

int nhb, nor, nit, nth, nms, level, seed;
double beta;

static double hmass = 1.15;

void M(spinor_field *out, spinor_field *in)
{
  g5Dphi_sq(-hmass, out, in);
}

int main(int argc, char *argv[])
{
  int i;
  int return_value = 0;
  double tau;
  spinor_field *s1, *s2;
  spinor_field *res, *res2;

  mshift_par par;

  int cgiters;

  eva_prec e_par;

  logger_map("DEBUG", "debug");

  setup_process(&argc, &argv);

  setup_gauge_fields();

  lprintf("MAIN", 0, "Generating a random gauge field... ");
  fflush(stdout);
  random_u(u_gauge);
  lprintf("MAIN", 0, "done.\n");
  represent_gauge_field();

  par.n = 6;
  par.shift = (double *)malloc(sizeof(double) * (par.n));
  par.err2 = 1.e-20;
  par.max_iter = 0;
  res = alloc_spinor_field_f(2 * par.n + 2, &glattice);
  s1 = res + par.n;
  s2 = s1 + 1;
  res2 = s2 + 1;

  par.shift[0] = +0.;
  par.shift[1] = -0.21;
  par.shift[2] = -0.05;
  par.shift[3] = -0.01;
  par.shift[4] = -0.15;
  par.shift[5] = -0.07;

  gaussian_spinor_field(s1);

  /* EVA prec parameters */
  e_par.nevt = 14;
  e_par.nev = 8;
  e_par.kmax = 400;
  e_par.maxiter = 100;
  e_par.omega1 = 1.e-10;
  e_par.omega2 = 5.e-7;

  set_def_matrix(&e_par, &M, s1->type);

  /* TEST CG_M */

  lprintf("CG TEST", 0, "Testing CG multishift PRE\n");
  lprintf("CG TEST", 0, "---------------------\n");

  cgiters = cg_mshift_def(&par, &M, &eva_def, &eva_def_inv, s1, res);
  lprintf("CG TEST", 0, "Converged in %d iterations\n", cgiters);
  for (i = 0; i < par.n; ++i)
  {
    M(s2, &res[i]);
    spinor_field_mul_add_assign_f(s2, -par.shift[i], &res[i]);
    spinor_field_sub_assign_f(s2, s1);
    tau = spinor_field_sqnorm_f(s2) / spinor_field_sqnorm_f(s1);
    lprintf("CG TEST", 0, "test cg[%d] = %e (req. %e)\n", i, tau, par.err2);
    if (tau > par.err2)
      return_value += 1;
  }

  lprintf("cg test", 0, "testing CG multishift DOUBLE\n");
  lprintf("cg test", 0, "---------------------\n");

  cgiters = cg_mshift(&par, &M, s1, res2);
  lprintf("CG TEST", 0, "Converged in %d iterations\n", cgiters);
  for (i = 0; i < par.n; ++i)
  {
    M(s2, &res2[i]);
    spinor_field_mul_add_assign_f(s2, -par.shift[i], &res2[i]);
    spinor_field_sub_assign_f(s2, s1);
    tau = spinor_field_sqnorm_f(s2) / spinor_field_sqnorm_f(s1);
    lprintf("CG TEST", 0, "test cg[%d] = %e (req. %e)\n", i, tau, par.err2);
    if (tau > par.err2)
      return_value += 1;
  }

  lprintf("CG TEST", 0, "Consistency FLOAT <-> DOUBLE\n");
  lprintf("CG TEST", 0, "---------------------\n");

  for (i = 0; i < par.n; ++i)
  {
    spinor_field_sub_assign_f(res + i, res2 + i);
    M(s2, res + i);
    spinor_field_mul_add_assign_f(s2, -par.shift[i], &res[i]);
    tau = spinor_field_sqnorm_f(s2) / spinor_field_sqnorm_f(s1);
    lprintf("CG TEST", 0, "test double-single[%d] = %e (req. %e)\n", i, tau, par.err2);
    if (tau > 10 * par.err2)
      return_value += 1;
  }

  free_spinor_field_f(res);

  free(par.shift);

  finalize_process();

  return return_value;
}
