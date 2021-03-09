/*******************************************************************************
*
* File  inv_hermNg.c
*
* Function to calculate determinant of an hermitian-matrix using LU-decomposition
* LU-decomposition. Modified from the real number version of Numerical Recipes.
*
* Ari Hietanen
*
*******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "utils.h"
#include "suN.h"
#include "representation.h"

#ifndef GAUGE_SON

#ifdef WITH_QUATERNIONS
void inv_hermNg(suNg *a)
{
  suNg b;
  double norm;

  b = *a;
  _suNg_quat_det(norm, b);
  norm = 1 / (norm);
  _suNg_dagger(b, b);
  _suNg_mul_assign(b, norm);
  *a = b;
}

#else
void inv_hermNg(suNg *a)
{
  suNg b;
  double complex col[NG];
  int indx[NG];
  double d;
  int i, j;
  b = *a;
  ludcmp(b.c, indx, &d, NG);
  for (j = 0; j < NG; j++)
  {
    for (i = 0; i < NG; i++)
    {
      _complex_0(col[i]);
    }
    _complex_1(col[j]);
    lubksb(b.c, indx, col, NG);
    for (i = 0; i < NG; i++)
      a->c[i * NG + j] = col[i];
  }
}

#endif

void ludcmp(double complex *a, int *indx, double *d, int N)
{
  const double tiny = 1.0e-20;
  int i, j, k, imax;
  double big, tmp, dum;
  double vv[N];
  double complex ctmp, csum;
  *d = 1;
  for (j = 0; j < N; ++j)
  {
    big = 0;
    for (i = 0; i < N; ++i)
    {
      tmp = _complex_prod_re(a[j * N + i], a[j * N + i]);
      if (tmp > big)
        big = tmp;
    }
    error(big == 0.0, 1, "ludcmp", "Singular matrix");
    vv[j] = 1 / sqrt(big);
  }
  imax = 0;
  for (j = 0; j < N; j++)
  { //Loop ower columns Crout's method
    //Calculate upper triangular matrix
    for (i = 0; i < j; i++)
    {
      csum = a[i * N + j];
      for (k = 0; k < i; k++)
      {
        csum -= a[i * N + k] * a[k * N + j];
      }
      a[i * N + j] = csum;
    }
    //Calculate lower triangular matrix and
    // find the largest pivot element
    big = 0;
    for (i = j; i < N; ++i)
    {
      csum = a[i * N + j];
      for (k = 0; k < j; ++k)
      {
        csum -= a[i * N + k] * a[k * N + j];
      }
      a[i * N + j] = csum;
      dum = sqrt(_complex_prod_re(csum, csum));
      dum *= vv[i];
      if (dum >= big)
      {
        big = dum;
        imax = i;
      }
    }
    if (j != imax)
    { //Do we need to interchage rows
      for (k = 0; k < N; ++k)
      {
        ctmp = a[imax * N + k];
        a[imax * N + k] = a[j * N + k];
        a[j * N + k] = ctmp;
      }
      *d = -*d;
      vv[imax] = vv[j];
    }
    indx[j] = imax;

    tmp = sqrt(_complex_prod_re(a[j * N + j], a[j * N + j]));
    if (tmp == 0)
      a[j * N + j] = tiny;

    //Divide by pivot element
    if (j != N - 1)
    {
      _complex_inv(csum, a[j * N + j]);
      for (int i1 = j + 1; i1 < N; ++i1)
      {
        a[i1 * N + j] = csum * a[i1 * N + j];
      }
    }
  }
}

void lubksb(double complex *a, int *indx, double complex *b, int N)
{
  int i, ii, ip, j;
  double complex csum, ctmp;
  double tmp;
  ii = 0;
  for (i = 0; i < N; ++i)
  {
    ip = indx[i];
    csum = b[ip];
    b[ip] = b[i];
    if (ii != 0)
    {
      for (j = ii - 1; j < i; ++j)
      {
        _complex_mul(ctmp, a[i * N + j], b[j]);
        _complex_sub_assign(csum, ctmp);
      }
    }
    else
    {
      tmp = _complex_prod_re(csum, csum);
      if (tmp != 0.0)
        ii = i + 1;
    }
    b[i] = csum;
  }
  for (i = N - 1; i >= 0; i--)
  {
    csum = b[i];
    for (j = i + 1; j < N; ++j)
    {
      _complex_mul(ctmp, a[i * N + j], b[j]);
      _complex_sub_assign(csum, ctmp);
    }
    _complex_div(b[i], csum, a[i * N + i]);
  }
}

#else

#endif
