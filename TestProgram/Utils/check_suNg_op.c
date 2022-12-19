/*
 * NOCOMPILE = !GAUGE_SUN
 */
#define MAIN_PROGRAM

#include <stdio.h>
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
#include "communications.h"

void print_cmplx(hr_complex c)
{
  lprintf("MAIN", 0, "%.14e", creal(c));
  if (cimag(c) > 0)
  {
    lprintf("MAIN", 0, "+I*%.14e", cimag(c));
  }
  else
  {
    lprintf("MAIN", 0, "-I*%.14e", fabs(cimag(c)));
  }
}

double norm_suNg_minus_id(suNg *a)
{
  int i, j;
  double r = 0.;
  for (i = 0; i < NG; i++)
  {
    for (j = 1; j < NG; j++)
    {
#ifdef WITH_QUATERNIONS
      r += fabs(a->c[i * NG + j]);
#else
      r += cabs(a->c[i * NG + j]);
#endif
      if (i == j)
        r -= 1.0;
    }
  }
  return r;
}
void printML(suNg *a)
{
  int i, j;
  lprintf("MAIN", 0, "{{");
  for (i = 0; i < NG; i++)
  {
    print_cmplx(a->c[i * NG]);
    for (j = 1; j < NG; j++)
    {
      lprintf("MAIN", 0, ",");
      print_cmplx(a->c[i * NG + j]);
    }
    if (i < NG - 1.0)
      lprintf("MAIN", 0, "},{");
  }
  lprintf("MAIN", 0, "}}\n\n");
}
#define niter 20

int main(int argc, char *argv[])
{
  int return_value = 0;
  hr_complex det, detinv;
  double test;
  suNg a, b, c;
  setup_process(&argc, &argv);
 
  for (int i = 0; i < niter; i++)
  {
    random_suNg(&a);
    random_suNg(&b);

    det_Cmplx_Ng(&det, &a);
    if (fabs(creal(det - 1.0)) > 1e-14 || fabs(cimag(det)) > 1e-14)
    {
      return_value += 1;
      lprintf("ERROR", 0, "fail\n");
    }
    if (return_value != 0 || i == niter - 1)
    {
      lprintf("MAIN", 0, "A=Random suNg\n");
      printML(&a);
      lprintf("MAIN", 0, "Det(A):");
      print_cmplx(det);
      lprintf("MAIN", 0, " should be 1.0\n");
    }

    _suNg_add_assign(a, b); /* a = a+ b */
    det_Cmplx_Ng(&det, &a);
    if (fabs(creal(det - 1.0)) < 1e-10 || fabs(cimag(det)) > 1e-14)
    {
      return_value += 1;
      lprintf("ERROR", 0, "fail\n");
    }
    if (return_value != 0 || i == niter - 1)
    {
      lprintf("MAIN", 0, "-----------------------------------\n");
      lprintf("MAIN", 0, "C=Sum of random suNg\n");
      printML(&a);
      lprintf("MAIN", 0, "Det: ");
      print_cmplx(det);
      lprintf("MAIN", 0, " should not be 1.0\n");
    }

    ranlxd((double *)(a.c), 2 * NG * NG);
    det_Cmplx_Ng(&det, &a);
    if ((fabs(creal(det - 1.0)) < 1e-10) && (fabs(cimag(det)) < 1e-10))
    {
      return_value += 1;
      lprintf("ERROR", 0, "fail1\n");
    }
    b = a;
    inv_Cmplx_Ng(&b);
    det_Cmplx_Ng(&detinv, &b);
    if (fabs(creal(detinv - 1.0)) < 1e-10 && fabs(cimag(detinv)) < 1e-10)
    {
      return_value += 1;
      lprintf("ERROR", 0, "fail2\n");
    }
    if (fabs(creal(det * detinv - 1.0)) > 1e-14 || fabs(cimag(det * detinv)) > 1e-14)
    {
      return_value += 1;
      lprintf("ERROR", 0, "fail3\n");
    }
    _suNg_times_suNg(c, a, b);
    test = norm_suNg_minus_id(&c);
    if (test > 1e-12)
    {
      return_value += 1;
      lprintf("ERROR", 0, "fail4 %.14e\n", test);
    }
    if (return_value != 0 || i == niter - 1)
    {
      lprintf("MAIN", 0, "-----------------------------------\n");
      lprintf("MAIN", 0, "A=Random Matrix (Not Herm, Not SU(N))\n");
      printML(&a);
      lprintf("MAIN", 0, "Det A: ");
      print_cmplx(det);
      lprintf("MAIN", 0, " should not be 1.0\n");
      lprintf("MAIN", 0, "B=A^-1 Inverse of the random matrix\n");
      printML(&b);
      lprintf("MAIN", 0, "Det B: ");
      print_cmplx(detinv);
      lprintf("MAIN", 0, " should not be 1.0\n");
      lprintf("MAIN", 0, "Det A * Det B: ");
      print_cmplx(det * detinv);
      lprintf("MAIN", 0, " should  be 1.0\n\n");
      lprintf("MAIN", 0, "C= A^-1 * A Should be the identity matrix:\n");
      printML(&c);
    }
  }

  global_sum_int(&return_value,1);

  finalize_process();

  return return_value;
}
