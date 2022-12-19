/*******************************************************************************
*
* NOCOMPILE= !NG==3
* NOCOMPILE= !REPR_ANTISYMMETRIC
* NOCOMPILE= WITH_MPI
*
* Test of modules
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
#include "global.h"
#include "suN.h"
#include "suN_types.h"
#include "linear_algebra.h"
#include "inverters.h"
#include "representation.h"
#include "utils.h"
#include "update.h"
#include "setup.h"
#include "logger.h"

int main(int argc, char *argv[])
{

  int return_value = 0;
  double test_val;
  /* setup process id and communications */
  logger_map("DEBUG", "debug");

  setup_process(&argc, &argv);

  setup_gauge_fields();

  suNg A;
  suNf a, s, tmp;
  int i, j;
  hr_complex zp, zm;

  lprintf("MAIN", 0, "Gauge group: SU(%d)\n", NG);
  lprintf("MAIN", 0, "Fermion representation: dim = %d\n", NF);
  lprintf("MAIN", 0, "Check 2AS = fund* for SU(3)\n");
  lprintf("MAIN", 0, "\n");

  random_suNg(&A);
  _group_represent2(&a, &A);

  _complex_1(zp);
  _complex_mulr(zm, -1.0, zp);

  _suNf_zero(s);
  (s).c[2] = zp;
  (s).c[4] = zm;
  (s).c[6] = zp;
  _suNf_times_suNf(tmp, a, s);
  _suNf_times_suNf(a, s, tmp);

  lprintf("MAIN", 0, "fundamental representation:\n");
  for (i = 0; i < NF; i++)
  {
    for (j = 0; j < NF; j++)
    {
      lprintf("MAIN", 0, "%.4e + i %.4e   ", creal(A.c[i * NF + j]), cimag(A.c[i * NF + j]));
    }
    lprintf("MAIN", 0, "\n");
  }

  lprintf("MAIN", 0, "2AS representation:\n");
  for (i = 0; i < NF; i++)
  {
    for (j = 0; j < NF; j++)
    {
      lprintf("MAIN", 0, "%.4e + i %.4e   ", creal(a.c[i * NF + j]), cimag(a.c[i * NF + j]));
    }
    lprintf("MAIN", 0, "\n");
  }

  for (i = 0; i < NF * NF; i++)
    A.c[i] -= conj(a.c[i]);

  _suNg_sqnorm(test_val, A);

  if (sqrt(test_val) > 10.e-14)
    return_value = 1;

  finalize_process();
  return return_value;
}
