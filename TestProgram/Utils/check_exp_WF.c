/*******************************************************************************
*
* Check for the implementation of the exponential using the Horner scheme 
* by Fernando Romero-Lopez
*
*******************************************************************************/

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
#include "random.h"
#include "wilsonflow.h"
#include "setup.h"

int main(int argc, char *argv[])
{
  int i = 3, j = 0, evaluations = 1;
  suNg test, exptest, exptest2;
  double complex tr;

  setup_process(&argc, &argv);
  setup_gauge_fields();

  random_suNg(&test);

  // Traceless and real in the diagonal
  for (i = 0; i < NG; i++)
    test.c[i * (NG + 1)] = 0.5 * (test.c[i * (NG + 1)] - conj(test.c[i * (NG + 1)]));
  _suNg_trace(tr, test);
  test.c[NG * NG - 1] = test.c[NG * NG - 1] - tr;

  //Change offdiagonal elements to hermitian matrix
  for (i = 0; i < NG; i++)
    for (j = 0; j < i; j++)
    {
      test.c[NG * (i) + j] = - conj(test.c[NG * (j) + i]);
    }

  _suNg_trace(tr, test);

  if (creal(conj(tr)*tr) > 1.e-30)
  {
    lprintf("ERROR", 0, "random matrix not traceless!! Trace = %f\n", tr);
    return 1;
  }

  for (i = 0; i < evaluations; i++)
    WF_Exp_Taylor(&exptest, &test);

  for (i = 0; i < evaluations; i++)
    WF_Exp(&exptest2, &test);

  _suNg_sub_assign(exptest, exptest2);

  double norm = 0.;

  for (i = 0; i < NG * NG; i++)
    norm += conj(exptest.c[i]) * exptest.c[i];

  norm = sqrt(norm);

  int check = 0;

  lprintf("MAIN", 0, "Maximal normalized difference = %.2e\n", norm / (NG));
  lprintf("MAIN", 0, "(should be around 1*10^(-15) or so)\n");
  if (norm / (NG) > 1.e-14)
    check++;

  finalize_process();
  return check;
}
