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

int main(int argc, char *argv[])
{
  int return_value = 0;
  int i;
  double tau1, tau2;
  spinor_field *s1, *s2, *s3;
  spinor_field *res, *res_trunc;

  mshift_par par;

  int cgiters;
  logger_map("DEBUG", "debug");

  setup_process(&argc, &argv);

  setup_gauge_fields();

  lprintf("MAIN", 0, "Generating a random gauge field... ");
  lprintf("MAIN", 0, "done.\n");

  random_u(u_gauge);

  start_gf_sendrecv(u_gauge);
  complete_gf_sendrecv(u_gauge);

  represent_gauge_field();

  par.n = 6;
  par.shift = (double *)malloc(sizeof(double) * (par.n));
  par.err2 = 1.e-28;
  par.max_iter = 0;
  res = alloc_spinor_field_f(par.n * 2 + 3, &glattice);
  res_trunc = res + par.n;
  s1 = res_trunc + par.n;
  s2 = s1 + 1;
  s3 = s2 + 1;

  par.shift[0] = +0.1;
  par.shift[1] = -0.21;
  par.shift[2] = +0.05;
  par.shift[3] = -0.01;
  par.shift[4] = -0.15;
  par.shift[5] = -0.05;

  gaussian_spinor_field(s1);
  tau1 = spinor_field_sqnorm_f(s1);
  lprintf("QMR TEST", 0, "Norma iniziale: %e\n", tau1);

  /* TEST g5QMR_M */

  lprintf("QMR TEST", 0, "\n");
  lprintf("QMR TEST", 0, "Testing g5QMR multishift\n");
  lprintf("QMR TEST", 0, "------------------------\n");

  cgiters = g5QMR_mshift_trunc(&par, 20, &D, s1, res_trunc, res);
  lprintf("QMR TEST", 0, "Converged in %d iterations\n", cgiters);

  for (i = 0; i < par.n; ++i)
  {
    D(s2, &res[i]);
    D(s3, &res_trunc[i]);
    spinor_field_mul_add_assign_f(s2, -par.shift[i], &res[i]);
    spinor_field_mul_add_assign_f(s3, -par.shift[i], &res_trunc[i]);
    spinor_field_sub_assign_f(s2, s1);
    spinor_field_sub_assign_f(s3, s1);
    tau1 = spinor_field_sqnorm_f(s2) / spinor_field_sqnorm_f(s1);
    tau2 = spinor_field_sqnorm_f(s3) / spinor_field_sqnorm_f(s1);
    lprintf("QMR TEST", 0, " g5QMR[%d] = %e, trunc = %e (req. %e)\n", i, tau1, tau2, par.err2);
    if (tau1 > par.err2)
      return_value += 1;
  }

  free_spinor_field_f(res);
  free(par.shift);
  finalize_process();

  return return_value;
}
