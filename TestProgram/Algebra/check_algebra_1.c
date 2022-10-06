/*******************************************************************************
*
*  NOCOMPILE= WITH_MPI
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
#include "setup.h"
#include "logger.h"

#ifdef REPR_FUNDAMENTAL
static float C2 = (float)(NG * NG - 1) / (float)(2 * NG);
static float Tr = 0.5;
#endif

#ifdef REPR_ADJOINT
static float C2 = (float)NG;
static float Tr = (float)NG;
#endif

#ifdef REPR_ANTISYMMETRIC
static float C2 = (float)(NG - 2) * (NG + 1) / (float)NG;
static float Tr = (float)(NG - 2) / 2;
#endif

#ifdef REPR_SYMMETRIC
static float C2 = (float)(NG + 2) * (NG - 1) / (float)NG;
static float Tr = (float)(NG + 2) / 2;
#endif

static int dAdj = NG * NG - 1;
static float fund = (float)(NG * NG - 1) / (2 * (float)(NG));

int main(int argc, char *argv[])
{

      int return_value = 0;
      /* setup process id and communications */
      logger_map("DEBUG", "debug");

      setup_process(&argc, &argv);

      setup_gauge_fields();

      suNg_algebra_vector f[dAdj];
      suNg A, B, TMP, CAS;
      suNf a, b, tmp, cas;
      double tau, trace,test_val;
      int i, j;

      for (i = 0; i < dAdj; i++)
      {
            _algebra_vector_zero_g(f[i]);
            f[i].c[i] = 1.;
      }

      for (i = 0; i < dAdj; i++)
      {
            for (j = 0; j < dAdj; j++)
            {
                  _algebra_represent(a, f[i]);
                  _algebra_represent(b, f[j]);

                  _suNf_times_suNf(tmp, a, b);
                  _suNf_trace_re(trace, tmp);
                  lprintf("MAIN", 0, "tr_R (T[%d] T[%d]): %.4e ", i, j, trace);
                  if (i == j)
                  {
                        lprintf("MAIN", 0, "  [should be: %.4e]\n", -Tr);
                        test_val=-Tr;
                  }
                  else
                  {
                        lprintf("MAIN", 0, "  [should be: 0.00]\n");
                        test_val=0.;
                  }
                  if (fabs(test_val-trace) > 10.e-14)
                        return_value += 1;

                  _fund_algebra_represent(A, f[i]);
                  _fund_algebra_represent(B, f[j]);

                  _suNg_times_suNg(TMP, A, B);
                  _suNg_trace_re(trace, TMP);
                  lprintf("MAIN", 0, "tr_f (T[%d] T[%d]): %.4e ", i, j, trace);
                  if (i == j)
                  {
                        lprintf("MAIN", 0, "  [should be: %.4e]\n", -0.5);
                        test_val=-0.5;
                  }
                  else
                  {
                        lprintf("MAIN", 0, "  [should be: 0.00]\n");
                        test_val=0.;
                  }
                  if (fabs(test_val-trace) > 10.e-14)
                        return_value += 1;

            }
      }

      _algebra_represent(a, f[0]);
      _fund_algebra_represent(A, f[0]);
      _suNf_times_suNf(cas, a, a);
      _suNg_times_suNg(CAS, A, A);

      for (i = 1; i < dAdj; i++)
      {
            _algebra_represent(a, f[i]);
            _fund_algebra_represent(A, f[i]);

            _suNf_times_suNf(tmp, a, a);
            _suNf_add_assign(cas, tmp);
            _suNg_times_suNg(TMP, A, A);
            _suNg_add_assign(CAS, TMP);
      }

      _suNf_unit(tmp);
      _suNf_mul(tmp, C2, tmp);
      _suNf_add_assign(cas, tmp);
      _suNf_sqnorm(tau, cas);
      lprintf("MAIN", 0, "casimir check: %.4e\n", tau);
      lprintf("MAIN", 0, "(should be around 1*10^(-15) or so)\n\n");
      if (fabs(tau) > 10.e-14)
            return_value += 1;


      _suNg_unit(TMP);
      _suNg_mul(TMP, fund, TMP);
      _suNg_add_assign(CAS, TMP);
      _suNg_sqnorm(tau, CAS);
      lprintf("MAIN", 0, "casimir check: %.4e\n", tau);
      lprintf("MAIN", 0, "(should be around 1*10^(-15) or so)\n\n");
      if (fabs(tau) > 10.e-14)
            return_value += 1;

      finalize_process();
      return return_value;
}
