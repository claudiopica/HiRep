/******************************************************************************
 *
 * Checks the invesrion of of the Dphi eopre tm
 * Author: Antonio Rago & Fernando Romero Lopez
 * NOCOMPILE= BC_X_ANTIPERIODIC
 * NOCOMPILE= BC_Y_ANTIPERIODIC
 * NOCOMPILE= BC_Z_ANTIPERIODIC
 * NOCOMPILE= BASIC_SF
 * NOCOMPILE= ROTATED_SF
 * NOCOMPILE= FERMION_THETA
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

static double hmass = 0.1;

static spinor_field *tmp;

void M(spinor_field *out, spinor_field *in)
{
  g5Dphi_eopre_tw_sq(hmass,0.0, out, in);
}

void Morig(spinor_field *out, spinor_field *in)
{
  g5Dphi_eopre_sq(hmass, out, in);
}


int main(int argc, char *argv[])
{

  int i;
  double tau;
  spinor_field *s1, *s2;
  spinor_field *res, *res2;
  struct timeval start, end, etime;
  int return_value = 0;
  mshift_par par;

  int cgiters;

  logger_map("DEBUG", "debug");

  setup_process(&argc, &argv);

  setup_gauge_fields();

  lprintf("MAIN", 0, "Generating a random gauge field... ");
  random_u(u_gauge);

  represent_gauge_field();

  par.n = 6;
  par.shift = (double *)malloc(sizeof(double) * (par.n));
  par.err2 = 1.e-20;
  par.max_iter = 0;
  res = alloc_spinor_field_f(2 * par.n + 3, &glat_even);
  s1 = res + par.n;
  s2 = s1 + 1;
  tmp = s2 + 1;
  res2 = tmp + 1;

  par.shift[0] = -0.0;
  par.shift[1] = -0.1;
  par.shift[2] = -0.2;
  par.shift[3] = +0.1;
  par.shift[4] = +0.2;
  par.shift[5] = +0.3;

  gaussian_spinor_field(s1);

  /* TEST CG_M */


  lprintf("cg test", 0, "testing CG multishift DOUBLE\n");
  lprintf("cg test", 0, "---------------------\n");

  gettimeofday(&start, 0);
  cgiters = cg_mshift(&par, &Morig, s1, res2);
  gettimeofday(&end, 0);
  timeval_subtract(&etime, &end, &start);
  lprintf("CG TEST", 0, "Inversion took: %ld sec %ld usec\n", etime.tv_sec, etime.tv_usec);
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

  free_spinor_field_f(res);
  free(par.shift);

  finalize_process();

  return return_value;
}
