/* Exponential of a SU(N) matrix using the Caley Hamilton or Taylor representation */
/* arXiv:1006.4518 [hep-lat] */

#include "global.h"
#include "geometry.h"
#include "suN.h"
#include "suN_repr_func.h"
#include "memory.h"
#include "global.h"
#include "logger.h"
#include "update.h"
#include "utils.h"
#include "communications.h"
#include <math.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327
#endif

// Here I start modifying FRL
/*                                                                                                                                                                                    
 *  u = exp(X), with X traceless 
 */
#if (NG >= 3) && (NG <= 6)

static double *inverse_fact = NULL;
static double factorial(int N)
{

  int i;
  double fact = 1.;

  for (i = 1; i <= N; ++i)
  {
    fact *= i;
  }

  return fact;
}
#endif

#if (NG == 3)
static void suNg_Exp_NG3(suNg *u, suNg *Xin)
{
  int NN = 30, i = 0, j = 0;

  double complex p[NG - 1];
  suNg X0, X2, X3;

  if (inverse_fact == NULL)
  {
    inverse_fact = malloc(sizeof(double) * (NN + 1));
    for (i = 0; i < NN + 1; i++)
      inverse_fact[i] = 1. / factorial(i);
  }

  _suNg_times_suNg(X2, *Xin, *Xin);
  _suNg_times_suNg(X3, X2, *Xin);
  _suNg_unit(X0);

  _suNg_trace(p[0], X3);
  _suNg_trace(p[1], X2);

  p[0] = -p[0] / 3;
  p[1] = -p[1] / 2;

  double complex q[NG];
  for (i = 0; i < NG; i++)
    q[i] = 0.;

  double complex qlast;
  q[0] = inverse_fact[NN];

  for (i = NN - 1; i >= 0; i--)
  {
    qlast = q[NG - 1];
    q[NG - 1] = q[NG - 2];
    for (j = NG - 2; j > 0; j--)
      q[j] = q[j - 1] - p[j] * qlast;
    q[0] = inverse_fact[i] - p[0] * qlast;
  }

  _suNg_mul_add(*u, q[0], X0, q[1], *Xin);
  _suNg_mulc(X3, q[2], X2);
  _suNg_add_assign(*u, X3);
}
#endif

#if (NG == 4)
static void suNg_Exp_NG4(suNg *u, suNg *Xin)
{

  int NN = 30, i = 0, j = 0;

  double complex p[NG - 1];

  if (inverse_fact == NULL)
  {
    inverse_fact = malloc(sizeof(double) * (NN + 1));
    for (i = 0; i < NN + 1; i++)
      inverse_fact[i] = 1. / factorial(i);
  }
  suNg X0, X2, X3, X4;

  _suNg_times_suNg(X2, *Xin, *Xin);
  _suNg_times_suNg(X3, X2, *Xin);
  _suNg_times_suNg(X4, X3, *Xin);
  _suNg_unit(X0);

  _suNg_trace(p[0], X4);
  _suNg_trace(p[1], X3);
  _suNg_trace(p[2], X2);

  p[0] = -p[0] / 4 + p[2] * p[2] / 8;
  p[1] = -p[1] / 3;
  p[2] = -p[2] / 2;

  double complex q[NG];
  for (i = 0; i < NG; i++)
    q[i] = 0.;

  double complex qlast;
  q[0] = inverse_fact[NN];

  for (i = NN - 1; i >= 0; i--)
  {
    qlast = q[NG - 1];
    q[NG - 1] = q[NG - 2];
    for (j = NG - 2; j > 0; j--)
      q[j] = q[j - 1] - p[j] * qlast;
    q[0] = inverse_fact[i] - p[0] * qlast;
  }

  _suNg_mul_add(*u, q[0], X0, q[1], *Xin);
  _suNg_mulc(X4, q[2], X2);
  _suNg_add_assign(*u, X4);
  _suNg_mulc(X4, q[3], X3);
  _suNg_add_assign(*u, X4);
}
#endif

#if (NG == 5)
static void suNg_Exp_NG5(suNg *u, suNg *Xin)
{

  int NN = 30, i = 0, j = 0;
  double complex p[NG - 1];

  if (inverse_fact == NULL)
  {
    inverse_fact = malloc(sizeof(double) * (NN + 1));
    for (i = 0; i < NN + 1; i++)
      inverse_fact[i] = 1. / factorial(i);
  }

  suNg X0, X2, X3, X4, X5;

  _suNg_times_suNg(X2, *Xin, *Xin);
  _suNg_times_suNg(X3, X2, *Xin);
  _suNg_times_suNg(X4, X3, *Xin);
  _suNg_times_suNg(X5, X4, *Xin);
  _suNg_unit(X0);

  _suNg_trace(p[0], X5);
  _suNg_trace(p[1], X4);
  _suNg_trace(p[2], X3);
  _suNg_trace(p[3], X2);

  p[0] = -p[0] / 5 + p[3] * p[2] / 6;
  p[1] = -p[1] / 4 + p[3] * p[3] / 8;
  p[2] = -p[2] / 3;
  p[3] = -p[3] / 2;

  double complex q[NG];
  for (i = 0; i < NG; i++)
    q[i] = 0.;

  double complex qlast;
  q[0] = inverse_fact[NN];

  for (i = NN - 1; i >= 0; i--)
  {
    qlast = q[NG - 1];
    q[NG - 1] = q[NG - 2];
    for (j = NG - 2; j > 0; j--)
      q[j] = q[j - 1] - p[j] * qlast;
    q[0] = inverse_fact[i] - p[0] * qlast;
  }

  _suNg_mul_add(*u, q[0], X0, q[1], *Xin);
  _suNg_mulc(X5, q[2], X2);
  _suNg_add_assign(*u, X5);
  _suNg_mulc(X5, q[3], X3);
  _suNg_add_assign(*u, X5);
  _suNg_mulc(X5, q[4], X4);
  _suNg_add_assign(*u, X5);
}
#endif

#if (NG == 6)
static void suNg_Exp_NG6(suNg *u, suNg *Xin)
{

  int NN = 30, i = 0, j = 0;

  double complex p[NG - 1];

  if (inverse_fact == NULL)
  {
    inverse_fact = malloc(sizeof(double) * (NN + 1));
    for (i = 0; i < NN + 1; i++)
      inverse_fact[i] = 1. / factorial(i);
  }

  suNg X0, X2, X3, X4, X5, X6;

  _suNg_times_suNg(X2, *Xin, *Xin);
  _suNg_times_suNg(X3, X2, *Xin);
  _suNg_times_suNg(X4, X3, *Xin);
  _suNg_times_suNg(X5, X4, *Xin);
  _suNg_times_suNg(X6, X5, *Xin);
  _suNg_unit(X0);

  _suNg_trace(p[0], X6);
  _suNg_trace(p[1], X5);
  _suNg_trace(p[2], X4);
  _suNg_trace(p[3], X3);
  _suNg_trace(p[4], X2);

  p[0] = -p[0] / 6 + p[4] * p[2] / 8 + p[3] * p[3] / 18 - p[4] * p[4] * p[4] / 48;
  p[1] = -p[1] / 5 + p[4] * p[3] / 6;
  p[2] = -p[2] / 4 + p[4] * p[4] / 8;
  p[3] = -p[3] / 3;
  p[4] = -p[4] / 2;

  double complex q[NG];
  for (i = 0; i < NG; i++)
    q[i] = 0.;

  double complex qlast;
  q[0] = inverse_fact[NN];

  for (i = NN - 1; i >= 0; i--)
  {
    qlast = q[NG - 1];
    q[NG - 1] = q[NG - 2];
    for (j = NG - 2; j > 0; j--)
      q[j] = q[j - 1] - p[j] * qlast;
    q[0] = inverse_fact[i] - p[0] * qlast;
  }

  _suNg_mul_add(*u, q[0], X0, q[1], *Xin);
  _suNg_mulc(X6, q[2], X2);
  _suNg_add_assign(*u, X6);
  _suNg_mulc(X6, q[3], X3);
  _suNg_add_assign(*u, X6);
  _suNg_mulc(X6, q[4], X4);
  _suNg_add_assign(*u, X6);
  _suNg_mulc(X6, q[5], X5);
  _suNg_add_assign(*u, X6);
}
#endif

#if (NG == 2)
/*
 *  u = exp(X)
 *
 * I AM ASSUMING
 * X^dag = -X
 * tr X = 0
 */
static void suNg_Exp_NG2(suNg *u, suNg *Xin)
{
  suNg_algebra_vector h, v;

  h.c[0] = cimag(Xin->c[1]);
  h.c[1] = creal(Xin->c[1]);
  h.c[2] = cimag(Xin->c[0]);

  double z = sqrt(h.c[0] * h.c[0] + h.c[1] * h.c[1] + h.c[2] * h.c[2]);
  double s = 1.;
  if (z > 1e-16)
    s = sin(z) / z;
  double c = cos(z);
  v.c[0] = h.c[0] * s;
  v.c[1] = h.c[1] * s;
  v.c[2] = h.c[2] * s;

  u->c[0] = c + I * v.c[2];
  u->c[1] = v.c[1] + I * v.c[0];
  u->c[2] = -v.c[1] + I * v.c[0];
  u->c[3] = c + I * -v.c[2];
}
#endif

void suNg_Exp_Taylor(suNg *u, suNg *Xin)
{
  suNg Xk, tmp;
  _suNg_unit(*u);
  _suNg_unit(Xk);

  int k = 1;
  double error;
  while (1)
  {
    _suNg_times_suNg(tmp, Xk, *Xin);
    _suNg_mul(Xk, 1. / k, tmp);
    k++;
    _suNg_add_assign(*u, Xk);

    _suNg_sqnorm(error, Xk);
    if (error < 1e-28)
      break;
  }
}

inline void suNg_Exp(suNg *u, suNg *Xin)
{
#if (NG == 2)
  suNg_Exp_NG2(u, Xin);
#elif (NG == 3)
  suNg_Exp_NG3(u, Xin);
#elif (NG == 4)
  suNg_Exp_NG4(u, Xin);
#elif (NG == 5)
  suNg_Exp_NG5(u, Xin);
#elif (NG == 6)
  suNg_Exp_NG6(u, Xin);
#else
  suNg_Exp_Taylor(u, Xin);
#endif
}


#ifdef GAUGE_SON
void ExpX(double dt, suNg_algebra_vector *h, suNg *r)
{
	error(0 == 0, 1, "ExpX [suN_epx.c]", "This function has yet not been implementd for SON");
}
#else
void ExpX(double dt, suNg_algebra_vector *h, suNg *u)
{
#ifdef WITH_QUATERNIONS
	suNg v_tmp, u_tmp;

	u_tmp = *u;
	_suNg_exp(dt, *h, v_tmp);
	_suNg_times_suNg(*u, v_tmp, u_tmp);
#else //WITH_QUATERNIONS
	suNg tmp1, tmp2;

	_fund_algebra_represent(tmp1, *h);
	_suNg_mul(tmp1, dt, tmp1);

	suNg_Exp(&tmp2, &tmp1);
	tmp1 = *u;

	_suNg_times_suNg(*u, tmp2, tmp1);

#endif //WITH_QUATERNIONS
}
#endif