/******************************************************************************
*
* Test of hermiticity
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
#include "setup.h"

#include "communications.h"

static double hmass = 0.1;

void DD(spinor_field *out, spinor_field *in)
{
  Dphi(hmass, out, in);
}

void HH(spinor_field *out, spinor_field *in)
{
  g5Dphi(-hmass, out, in);
}

void MM(spinor_field *out, spinor_field *in)
{
#ifdef UPDATE_EO
  g5Dphi_eopre_sq(-hmass, out, in);
#else
  g5Dphi_sq(-hmass, out, in);
#endif
}

int test_herm(spinor_operator S, char *name)
{
  spinor_field *s1, *s2, *s3, *s4;
  double tau;
  int return_val = 0;
#ifdef UPDATE_EO
  s1 = alloc_spinor_field_f(4, &glat_even);
#else
  s1 = alloc_spinor_field_f(4, &glattice);
#endif
  s2 = s1 + 1;
  s3 = s2 + 1;
  s4 = s3 + 1;

  lprintf("RESULT", 0, "Test if %s is hermitean: ", name);

  gaussian_spinor_field(s1);
  gaussian_spinor_field(s2);
  S(s3, s1);
  S(s4, s2);

  tau = spinor_field_prod_re_f(s2, s3);
  tau -= spinor_field_prod_re_f(s4, s1);
  tau += spinor_field_prod_im_f(s2, s3);
  tau -= spinor_field_prod_im_f(s4, s1);
  tau /= sqrt(spinor_field_sqnorm_f(s1));
  tau /= sqrt(spinor_field_sqnorm_f(s2));
  if (fabs(tau) > 1.e-14)
  {
    lprintf("RESULT", 0, "FAILED ");
    return_val = 1;
  }
  else
    lprintf("RESULT", 0, "OK ");
  lprintf("RESULT", 0, "[norm = %e]\n", tau);

  free_spinor_field_f(s1);
  return return_val;
}

int main(int argc, char *argv[])
{
  int return_value;
  /* setup process id and communications */
  logger_map("DEBUG", "debug");

  setup_process(&argc, &argv);

  setup_gauge_fields();

  lprintf("MAIN", 0, "Generating a random gauge field... ");
  random_u(u_gauge);
  lprintf("MAIN", 0, "done.\n");

  start_gf_sendrecv(u_gauge);

  represent_gauge_field();

  return_value=test_herm(&MM, "M");

  finalize_process();

  return return_value;
}
