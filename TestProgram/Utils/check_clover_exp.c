/*******************************************************************************
*
* NOCOMPILE= !WITH_EXPCLOVER
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
#include "clover_exp.h"
#include "clover_tools.h"

static void random_mat(suNfc *u)
{
  double gr[2 * NF * NF];
  int i;
  gauss(gr, 2 * NF * NF);
  for (i = 0; i < NF * NF; i++)
  {
    u->c[i] = gr[2 * i] + I * gr[2 * i + 1];
  }
}

int main(int argc, char *argv[])
{

  int i, j, return_value = 0;
  suNfc test[4], test2[4];
  suNfc exptest[4], exptest2[4];
  double tr=0;
  double mass=-3.6;

  setup_process(&argc, &argv);
  setup_gauge_fields();

  random_mat(&test[0]);
  random_mat(&test[3]);
  random_mat(&test[1]);

  _suNfc_dagger(test[2], test[1]);
  evaluate_sw_order(&mass);

  // Traceless and real in the diagonal
  for (i = 0; i < NF; i++)
    test[0].c[i * (NF + 1)] = 0.5 * (test[0].c[i * (NF + 1)] + conj(test[0].c[i * (NF + 1)]));
  _suNfc_trace(tr, test[0]);
  test[0].c[NF * NF - 1] = test[0].c[NF * NF - 1] - tr;
  //Change offdiagonal elements to hermitian matrix
  for (i = 0; i < NF; i++)
    for (j = 0; j < i; j++)
    {
      test[0].c[NF * (i) + j] = conj(test[0].c[NF * (j) + i]);
    }

  for (i = 0; i < NF; i++)
    test[3].c[i * (NF + 1)] = 0.5 * (test[3].c[i * (NF + 1)] + conj(test[3].c[i * (NF + 1)]));
  _suNfc_trace(tr, test[3]);
  test[3].c[NF * NF - 1] = test[3].c[NF * NF - 1] - tr;
  //Change offdiagonal elements to hermitian matrix
  for (i = 0; i < NF; i++)
    for (j = 0; j < i; j++)
    {
      test[3].c[NF * (i) + j] = conj(test[3].c[NF * (j) + i]);
    }

  _suNfc_trace(tr, test[0]);
  if (tr > 1.e-20)
  {
    lprintf("ERROR", 0, "random matrix not traceless!! Trace = %f\n", tr);
    return_value++;
  }
  _suNfc_trace(tr, test[3]);
  if (tr > 1.e-20)
  {
    lprintf("ERROR", 0, "random matrix not traceless!! Trace = %f\n", tr);
    return_value++;
  }
  for (i = 0; i < 4; i++)
  {
    test2[i] = test[i];
  }

  lprintf("MAIN", 0, "Evaluating the clover exp via Taylor expansion\n");

  clover_exp_taylor(test2, exptest2);

  lprintf("MAIN", 0, "Evaluating the clover exp via Cayley Hamilton repr\n");
  clover_exp(test, exptest);

  lprintf("MAIN", 0, "Comparing the resutls\n");

  for (i = 0; i < 4; i++)
  {
    _suNf_sub_assign(exptest[i], exptest2[i]);
  }

  double norm = 0.;

  for (j = 0; j < 4; j++)
  {
    norm = 0.;
    for (i = 0; i < NF * NF; i++)
      norm += _complex_prod(exptest[j].c[i], exptest[j].c[i]);
    norm = sqrt(norm);
    lprintf("MAIN", 0, "NORM of the difference %2.8e [should be around 1e-15] \n", norm);

    if (norm > NF * 1.0e-13)
      return_value++;
  }

  finalize_process();
  
  return return_value;
}
