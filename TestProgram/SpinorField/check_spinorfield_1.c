
/******************************************************************************
 *
 * File check11.c
 *
 * Consistency checks on the programs in the module linalg
 *
 * Author: luigi del debbio <luigi.del.debbio@ed.ac.uk>
 *
 ******************************************************************************/

#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "suN.h"
#include "random.h"
#include "error.h"
#include "geometry.h"
#include "memory.h"
#include "update.h"
#include "linear_algebra.h"
#include "representation.h"
#include "global.h"
#include "logger.h"
#include "utils.h"
#include "io.h"
#include "setup.h"

#include "communications.h"

#define MAX_ROTATE 50

static double complex v[25];
static double EPSILON = 1.e-12;
static spinor_field *ppk[5];

static int initr = 0;

static suNf_spinor *psi;

static void alloc_ws_rotate(void)
{
  psi = calloc(MAX_ROTATE, sizeof(suNf_spinor));

  error((psi == NULL), 1, "alloc_ws_rotate [linalg.c]",
        "Unable to allocate workspace");

  initr = 1;
}

static void rotate_ptr(int n, spinor_field *pkk[], double complex vl[])
{
  if (initr == 0)
    alloc_ws_rotate();

  error((n < 1) || (n > MAX_ROTATE), 1, "rotate [eva.c]",
        "Parameter n is out of range");

  for (int i = 0; i < n; i++)
    error((*pkk)->type != (*(pkk + i))->type, 1, "not available", "Spinors don't match!");

#undef _OMP_PRAGMA //This doesn't works with multiple threads
#define _OMP_PRAGMA(s)

  _MASTER_FOR(pkk[0]->type, ix)
  {
    for (int k = 0; k < n; k++)
    {
      suNf_spinor *pk = &(psi[k]);
      suNf_spinor *pj = _FIELD_AT(pkk[0], ix);
      double complex *z = &vl[k];

      _spinor_mulc_f(*pk, *z, *pj);

      for (int j = 1; j < n; j++)
      {
        pj = _FIELD_AT(pkk[j], ix);
        z += n;

        _spinor_mulc_add_assign_f(*pk, *z, *pj);
      }
    }

    for (int k = 0; k < n; k++)
      *_FIELD_AT(pkk[k], ix) = psi[k];
  }
}

static void project(spinor_field *pk, spinor_field *pl)
{
  double complex sp;

  sp = -spinor_field_prod_f(pl, pk);

  spinor_field_mulc_add_assign_f(pk, sp, pl);
}

static double normalize(spinor_field *ps)
{
  double r, ri;

  r = spinor_field_sqnorm_f(ps);
  r = sqrt(r);
  error(r < EPSILON, 1, "normalize [eva.c]", "vector has vanishing norm");

  ri = 1.0 / r;
  spinor_field_mul_f(ps, ri, ps);

  return (double)(r);
}

static double complex sp(spinor_field *pk, spinor_field *pl)
{

  double complex z = 0.0;

  _TWO_SPINORS_FOR_SUM(pk, pl, x, y)
  {
    for (int i = 0; i < (4 * NF); i++)
    {
      double complex *rpk = (double complex *)_SPINOR_PTR(pk) + i;
      double complex *rpl = (double complex *)_SPINOR_PTR(pl) + i;
      z += conj(*rpk) * (*rpl);
      /* x+=(double)((*rpk).re*(*rpl).re+(*rpk).im*(*rpl).im); */
      /* y+=(double)((*rpk).re*(*rpl).im-(*rpk).im*(*rpl).re); */
      //rpk+=1; //?? why these increment
      //rpl+=1; //??
    }
  }

#ifdef WITH_MPI
  global_sum((double *)&z, 2);
#endif

  return z;
}

int main(int argc, char *argv[])
{
  int i, j;
  double r;
  double rd, zsqd;
  double d, dmax;
  double complex w;
  double complex zd, wd;
  spinor_field *ws;
  spinor_field *pk, *pl;
  spinor_field *tmp;
  int return_value = 0;

  logger_map("DEBUG", "debug");

  /* setup process id and communications */
  setup_process(&argc, &argv);

  lprintf("CPTEST", 0, "spinor gsize=%d\n", glattice.gsize_spinor);
  lprintf("CPTEST", 0, "spinor nbuffers=%d\n", glattice.nbuffers_spinor);
  lprintf("CPTEST", 0, "spinor ncopies=%d\n", glattice.ncopies_spinor);
  lprintf("CPTEST", 0, "gauge gsize=%d\n", glattice.gsize_gauge);
  lprintf("CPTEST", 0, "gauge nbuffers=%d\n", glattice.nbuffers_gauge);
  lprintf("CPTEST", 0, "gauge ncopies=%d\n", glattice.ncopies_gauge);
  lprintf("CPTEST", 0, "lmp=%d\n", glattice.local_master_pieces);

  lprintf("LA TEST", 0, "Consistency of the programs in the module linalg\n");
  lprintf("LA TEST", 0, "------------------------------------------------\n");

  tmp = alloc_spinor_field_f(1, &glattice);
  ws = alloc_spinor_field_f(10, &glattice);

  for (i = 0; i < 10; i++)
    gaussian_spinor_field(&ws[i]);

  dmax = 0.0;

  for (i = 0; i < 10; i++)
  {
    pk = &ws[i];
    pl = &ws[9 - i];
    w = sp(pk, pl);

    zd = spinor_field_prod_f(pk, pl);
    rd = spinor_field_sqnorm_f(pk) * spinor_field_sqnorm_f(pl);
    d = (zd - w) * conj(zd - w);
    /* d=((zd.re-(double)w.re)*(zd.re-(double)w.re)+ */
    /*    (zd.im-(double)w.im)*(zd.im-(double)w.im)); */
    d = sqrt(d / rd);
    if (d > dmax)
      dmax = d;

    rd = spinor_field_prod_re_f(pk, pl);
    d = fabs(creal(zd) / rd - 1.0);
    if (d > dmax)
      dmax = d;

    zd = spinor_field_prod_f(pk, pk);
    rd = spinor_field_sqnorm_f(pk);

    d = fabs(cimag(zd) / rd);
    if (d > dmax)
      dmax = d;

    d = fabs(creal(zd) / rd - 1.0f);
    if (d > dmax)
      dmax = d;
  }
  lprintf("LA TEST", 0, "Check of spinor_field_prod, spinor_field_prod_re\n");
  lprintf("LA TEST", 0, "and spinor_field_sqnorm: %.2e\n\n", dmax);
  if (dmax > 1e-14)
  {
    lprintf("LA TEST", 0, "Test failed ?\n");
    return_value += 1;
  }

  dmax = 0.0;
  zd = 0.345 - I * 0.876;
  zsqd = zd * conj(zd);

  for (i = 0; i < 9; i++)
  {
    pk = &ws[i];
    pl = &ws[i + 1];

    wd = spinor_field_prod_f(pk, pl);
    rd = spinor_field_sqnorm_f(pk) + zsqd * spinor_field_sqnorm_f(pl) + 2.0 * (creal(zd * wd));

    spinor_field_mulc_add_assign_f(pk, zd, pl);

    d = fabs(rd / spinor_field_sqnorm_f(pk) - 1.0);
    if (d > dmax)
      dmax = d;
  }
  lprintf("LA TEST", 0, "Consistency of spinor_prod, norm_square\n");
  lprintf("LA TEST", 0, "and mulc_spinor_add: %.2e\n\n", dmax);
  if (dmax > 1e-14)
  {
    lprintf("LA TEST", 0, "Test failed ?\n");
    return_value += 1;
  }
  for (i = 0; i < 10; i++)
    gaussian_spinor_field(&ws[i]);

  dmax = 0.0;

  for (i = 0; i < 10; i++)
  {
    pk = &ws[i];

    if (i > 0)
    {
      pl = &ws[i - 1];
      project(pk, pl);
      zd = spinor_field_prod_f(pk, pl);

      d = (fabs(creal(zd)) +
           fabs(cimag(zd))) /
          sqrt(spinor_field_sqnorm_f(pk));

      if (d > dmax)
        dmax = d;
    }

    normalize(pk);
    rd = spinor_field_sqnorm_f(pk);

    d = fabs(rd - 1.0f);
    if (d > dmax)
      dmax = d;
  }

  lprintf("LA TEST", 0, "Consistency of spinor_prod, norm_square,\n");
  lprintf("LA TEST", 0, "normalize and project: %.2e\n\n", dmax);
  if (dmax > 1e-14)
  {
    lprintf("LA TEST", 0, "Test failed ?\n");
    return_value += 1;
  }

  for (i = 0; i < 5; i++)
  {
    pk = &ws[i];
    pl = &ws[i + 5];

    gaussian_spinor_field(pk);
    spinor_field_copy_f(pl, pk);

    for (j = 0; j < 5; j++)
    {
      v[5 * i + j] = 0.1234f * (double)(i ^ 2) - 0.8976f * (double)(j) + I * (0.2231f * (double)(i) + 0.9922f * (double)(j ^ 2));
    }

    ppk[i] = pl;
  }

  rotate_ptr(5, ppk, v);
  dmax = 0.0;

  for (i = 5; i < 10; i++)
  {
    pk = &ws[i];

    for (j = 0; j < 5; j++)
    {
      zd = -v[5 * j + (i - 5)];

      pl = &ws[j];
      spinor_field_mulc_add_assign_f(pk, zd, pl);
    }

    rd = spinor_field_sqnorm_f(pk);

    d = fabs(rd);
    if (d > dmax)
      dmax = d;
  }

  dmax /= spinor_field_sqnorm_f(&ws[0]);
  dmax = sqrt(dmax);

  lprintf("LA TEST", 0, "Consistency of mulc_spinor_add\n");
  lprintf("LA TEST", 0, "and rotate: %.2e\n\n", dmax);
  if (dmax > 1e-14)
  {
    lprintf("LA TEST", 0, "Test failed ?\n");
    return_value += 1;
  }

  dmax = 0.0;

  for (i = 0; i < 5; i++)
  {
    pk = &ws[i];
    pl = &ws[9 - i];
    gaussian_spinor_field(pk);
    spinor_field_copy_f(pl, pk);
    spinor_field_g5_f(tmp, pk);
    spinor_field_g5_f(pk, tmp);

    zd = -1.0;

    spinor_field_mulc_add_assign_f(pl, zd, pk);
    r = spinor_field_sqnorm_f(pl) / spinor_field_sqnorm_f(pk);
    d = sqrt(r);
    if (d > dmax)
      dmax = d;

    gaussian_spinor_field(pl);
    zd = spinor_field_prod_f(pk, pl);
    spinor_field_g5_f(pk, pk);
    spinor_field_g5_f(pl, pl);
    wd = spinor_field_prod_f(pk, pl);

    d = (fabs(creal(zd - wd)) + fabs(cimag(zd - wd))) /
        (fabs(creal(zd)) + fabs(cimag(zd)));
    if (d > dmax)
      dmax = d;
  }

  lprintf("LA TEST", 0, "Check of spinor_field_g5_f: %.2e\n\n", dmax);
  if (dmax > 1e-30)
  {
    lprintf("LA TEST", 0, "Test failed ?\n");
    return_value += 1;
  }

  dmax = 0.0;

  for (i = 0; i < 5; i++)
  {
    pk = &ws[i];
    pl = &ws[9 - i];
    gaussian_spinor_field(pk);
    spinor_field_copy_f(pl, pk);
    d = -2.5;
    spinor_field_lc1_f(d, pk, pl);

    zd = 1.5;
    spinor_field_mulc_add_assign_f(pk, zd, pl);
    d = spinor_field_sqnorm_f(pk) / spinor_field_sqnorm_f(pl);

    if (d > dmax)
      dmax = d;
  }

  lprintf("LA TEST", 0, "Check of lc1: %.2e\n\n", dmax);
  if (dmax > 1e-31)
  {
    lprintf("LA TEST", 0, "Test failed ?\n");
    return_value += 1;
  }
  dmax = 0.0;

  for (i = 0; i < 5; i++)
  {
    pk = &ws[i];
    pl = &ws[9 - i];
    gaussian_spinor_field(pk);
    spinor_field_copy_f(pl, pk);
    d = 1.0;
    r = 2.5;
    spinor_field_lc2_f(d, r, pk, pl);

    zd = -3.5;
    spinor_field_mulc_add_assign_f(pk, zd, pl);
    d = spinor_field_sqnorm_f(pk) / spinor_field_sqnorm_f(pl);

    if (d > dmax)
      dmax = d;
  }

  lprintf("LA TEST", 0, "Check of lc2: %.2e\n\n", dmax);
  if (dmax > 1e-31)
  {
    lprintf("LA TEST", 0, "Test failed ?\n");
    return_value += 1;
  }
  dmax = 0.0;

  for (i = 0; i < 5; i++)
  {
    pk = &ws[i];
    pl = &ws[9 - i];
    gaussian_spinor_field(pk);
    spinor_field_copy_f(pl, pk);
    d = 3.5;
    r = -1.5;
    spinor_field_lc3_f(d, r, pk, pl, pk);

    zd = -1.0;

    spinor_field_mulc_add_assign_f(pk, zd, pl);
    d = spinor_field_sqnorm_f(pk) / spinor_field_sqnorm_f(pl);

    if (d > dmax)
      dmax = d;
  }

  lprintf("LA TEST", 0, "Check of lc3: %.2e\n\n", dmax);
  if (dmax > 1e-31)
  {
    lprintf("LA TEST", 0, "Test failed ?\n");
    return_value += 1;
  }
  finalize_process();

  return return_value;
}
