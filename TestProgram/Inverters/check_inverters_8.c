/******************************************************************************
*
* NOCOMPILE= UPDATE_EO
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
  double tau;
  spinor_field *s1, *s2, *s3;

  g5QMR_fltacc_par par;
  mshift_par mpar;

  int cgiters;

  logger_map("DEBUG", "debug");

  setup_process(&argc, &argv);

  setup_gauge_fields();
  u_gauge_f_flt = alloc_gfield_f_flt(&glattice);

  lprintf("MAIN", 0, "Generating a random gauge field... ");
  lprintf("MAIN", 0, "done.\n");

  random_u(u_gauge);

  start_gf_sendrecv(u_gauge);
  complete_gf_sendrecv(u_gauge);

  represent_gauge_field();
  assign_ud2u_f();

  s1 = alloc_spinor_field_f(3, &glattice);
  s2 = s1 + 1;
  s3 = s2 + 1;

  gaussian_spinor_field(s1);

  /* TEST g5QMR_M */

  par.max_iter = 0;
  par.err2 = 1e-14;
  par.max_iter_flt = 0;
  par.err2_flt = 1e-6;

  set_dirac_mass(0.);

  lprintf("QMR TEST", 0, "\n");
  lprintf("QMR TEST", 0, "Testing g5QMR with single-precision acceleration\n");
  lprintf("QMR TEST", 0, "------------------------\n");

  cgiters = g5QMR_fltacc(&par, &D, &D_flt, s1, s3);
  lprintf("QMR TEST", 0, "Converged in %d iterations\n", cgiters);

  D(s2, s3);
  spinor_field_sub_assign_f(s2, s1);
  tau = spinor_field_sqnorm_f(s2) / spinor_field_sqnorm_f(s1);
  lprintf("QMR TEST", 0, "Res = %e\n", tau);

  if (tau > par.err2)
    return_value += 1;

  mpar.max_iter = 0;
  mpar.err2 = 1e-14;
  mpar.n = 1;
  double shift = 0;
  mpar.shift = &shift;

  lprintf("QMR TEST", 0, "\n");
  lprintf("QMR TEST", 0, "Testing g5QMR multishift\n");
  lprintf("QMR TEST", 0, "------------------------\n");

  spinor_field_zero_f(s3);
  cgiters = g5QMR_mshift(&mpar, &D, s1, s3);
  lprintf("QMR TEST", 0, "Converged in %d iterations\n", cgiters);

  D(s2, s3);
  spinor_field_sub_assign_f(s2, s1);
  tau = spinor_field_sqnorm_f(s2) / spinor_field_sqnorm_f(s1);
  lprintf("QMR TEST", 0, "Res = %e\n", tau);

  if (tau > par.err2)
    return_value += 1;

  free_spinor_field_f(s1);
  finalize_process();

  return return_value;
}
