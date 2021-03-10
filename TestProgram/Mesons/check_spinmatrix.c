
/******************************************************************************
 *
 * File check_spinmatrix.c
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
#include "spin_matrix.h"
#include "communications.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288419716939937510
#endif

//#error "Old version of Mesons, it should be updated"
static void print_suNf_spinor(suNf_spinor r, const char name[])
{
  int i, j;
  lprintf("MAIN", 0, "%s = \n", name);
  for (i = 0; i < 4; i++)
  {
    lprintf("MAIN", 0, "[ ");
    for (j = 0; j < NG; j++)
    {
      lprintf("MAIN", 0, "(%.2f,%.2f) ", creal(r.c[i].c[j]), cimag(r.c[i].c[j]));
    }
    lprintf("MAIN", 0, "]\n");
  }
}
static void print_spin_matrix(suNf_spin_matrix S, const char name[])
{
  int i, j, k;
  lprintf("MAIN", 0, "%s = \n", name);
  for (k = 0; k < NF; k++)
  {
    lprintf("MAIN", 0, "color %i \n", k);
    for (i = 0; i < 4; i++)
    {
      lprintf("MAIN", 0, "[ ");
      for (j = 0; j < 4; j++)
      {
        lprintf("MAIN", 0, "(%.2f,%.2f) ", creal(S.c[i].c[j].c[k]), cimag(S.c[i].c[j].c[k]));
      }
      lprintf("MAIN", 0, "]\n");
    }
  }
}
static double complex spin_matrix_norm_diff(suNf_spin_matrix sm1, suNf_spin_matrix sm2)
{
  int k;
  double complex z = 0.;
  suNf_spin_matrix sm_diff;

  _spinmatrix_sub(sm_diff, sm1, sm2);

  for (k = 0; k < 4; k++)
  {
    _spinor_prod_assign_f(z, sm_diff.c[k], sm_diff.c[k]);
  }
  return z;
}
static void transpose_mat(double complex At[4][4], double complex A[4][4])
{
  int i, j;
  for (i = 0; i < 4; i++)
  {
    for (j = 0; j < 4; j++)
    {
      At[i][j] = A[j][i];
    }
  }
}

#define assign_spin_matrix(sm, s)        \
  {                                      \
    _spinmatrix_assign_row(sm, s[0], 0); \
    _spinmatrix_assign_row(sm, s[1], 1); \
    _spinmatrix_assign_row(sm, s[2], 2); \
    _spinmatrix_assign_row(sm, s[3], 3); \
  }
#define set_zero_mat(A)          \
  {                              \
    int _i, _j;                  \
    for (_i = 0; _i < 4; _i++)   \
      for (_j = 0; _j < 4; _j++) \
      {                          \
        A[_i][_j] = 0.;          \
      }                          \
  }
#define mult_mat(r, A, B)                      \
  {                                            \
    int _i, _j, _k;                            \
    double complex wm[4][4];                   \
    for (_i = 0; _i < 4; _i++)                 \
      for (_j = 0; _j < 4; _j++)               \
      {                                        \
        _complex_0(wm[_i][_j]);                \
        for (_k = 0; _k < 4; _k++)             \
        {                                      \
          wm[_i][_j] += A[_i][_k] * B[_k][_j]; \
        }                                      \
      }                                        \
    for (_i = 0; _i < 4; _i++)                 \
      for (_j = 0; _j < 4; _j++)               \
      {                                        \
        r[_i][_j] = wm[_i][_j];                \
      }                                        \
  }
#define mult_mat_spinor(out, A, in)                      \
  {                                                      \
    int _a, _i, _j;                                      \
    for (_a = 0; _a < NF; _a++)                          \
      for (_i = 0; _i < 4; _i++)                         \
      {                                                  \
        _complex_0(out.c[_i].c[_a]);                     \
        for (_j = 0; _j < 4; _j++)                       \
        {                                                \
          out.c[_i].c[_a] += A[_i][_j] * in.c[_j].c[_a]; \
        }                                                \
      }                                                  \
  }

#define mult_mat_4spinor(sm, A, in)    \
  {                                    \
    suNf_spinor _out[4];                \
    mult_mat_spinor(_out[0], A, in[0]); \
    mult_mat_spinor(_out[1], A, in[1]); \
    mult_mat_spinor(_out[2], A, in[2]); \
    mult_mat_spinor(_out[3], A, in[3]); \
    assign_spin_matrix(sm, _out);       \
  }

#define _spinmatrix_assign_col(r, s, i) \
  (r).c[0].c[i] = (s).c[0];             \
  (r).c[1].c[i] = (s).c[1];             \
  (r).c[2].c[i] = (s).c[2];             \
  (r).c[3].c[i] = (s).c[3];

#define transpose_spin_matrix(out, in)      \
  {                                         \
    suNf_spinor tmp[4];                     \
    tmp[0] = in.c[0];                       \
    tmp[1] = in.c[1];                       \
    tmp[2] = in.c[2];                       \
    tmp[3] = in.c[3];                       \
    _spinmatrix_zero(out);                  \
    _spinmatrix_assign_col(out, tmp[0], 0); \
    _spinmatrix_assign_col(out, tmp[1], 1); \
    _spinmatrix_assign_col(out, tmp[2], 2); \
    _spinmatrix_assign_col(out, tmp[3], 3); \
  }

// compute the transpose of (A^T * sm^T)^T = sm A
#define mult_4spinor_mat(sm, A, in)          \
  {                                          \
    double complex At[4][4];                 \
    suNf_spin_matrix smtmp1, smtmp2, out;    \
    suNf_spinor tmp_spinor[4];               \
    transpose_mat(At, A);                    \
    assign_spin_matrix(smtmp1, in);          \
    transpose_spin_matrix(out, smtmp1);      \
    tmp_spinor[0] = out.c[0];                \
    tmp_spinor[1] = out.c[1];                \
    tmp_spinor[2] = out.c[2];                \
    tmp_spinor[3] = out.c[3];                \
    mult_mat_4spinor(smtmp2, A, tmp_spinor); \
    transpose_spin_matrix(sm, smtmp2);       \
  }

static void print_mat2(double complex mat[4][4], const char name[])
{
  int i, j;
  lprintf("MAIN", 0, "%s = \n", name);
  for (i = 0; i < 4; i++)
  {
    lprintf("MAIN", 0, "[ ");
    for (j = 0; j < 4; j++)
    {
      lprintf("MAIN", 0, "(%.2f,%.2f) ", creal(mat[i][j]), cimag(mat[i][j]));
    }
    lprintf("MAIN", 0, "]\n");
  }
}

int main(int argc, char *argv[])
{
  int sign, i;
  double complex ctest = 0.;
  suNf_spinor in[4];
  suNf_spin_matrix sma, smc;
  suNf_spin_matrix smb[16];
  double complex g[16][4][4];
  char *list_g[16] = {"g0", "g1", "g2", "g3", "g5", "g0g5", "g5g0", "g5g1", "g5g2", "g5g3", "g0g1", "g0g2", "g0g3", "g5g0g1", "g5g0g2", "g5g0g3"};

  int return_value = 0;

  /* setup process id and communications */
  setup_process(&argc, &argv);

  logger_map("DEBUG", "debug");

  set_zero_mat(g[6]);
  set_zero_mat(g[13]);
  set_zero_mat(g[14]);
  set_zero_mat(g[15]);

  g0_debug(g[0], &sign);
  g1_debug(g[1], &sign);
  g2_debug(g[2], &sign);
  g3_debug(g[3], &sign);
  g5_debug(g[4], &sign);
  g0g5_debug(g[5], &sign);
  mult_mat(g[6], g[4], g[0]) //  g5g0
      g5g1_debug(g[7], &sign);
  g5g2_debug(g[8], &sign);
  g5g3_debug(g[9], &sign);
  g0g1_debug(g[10], &sign);
  g0g2_debug(g[11], &sign);
  g0g3_debug(g[12], &sign);
  mult_mat(g[13], g[4], g[10])     //  g5g0g1
      mult_mat(g[14], g[4], g[11]) //  g5g0g2
      mult_mat(g[15], g[4], g[12]) //  g5g0g3

      for (i = 0; i < 5; i++)
  {
    print_mat2(g[i], list_g[i]);
  }

  // initialise random spinmatrix: sma

  ranlxd((double *)(&in[0]), 4 * sizeof(suNf_spinor) / sizeof(double));

  print_suNf_spinor(in[0], "test");
  assign_spin_matrix(sma, in);

  print_spin_matrix(sma, "test spinmatrix ");

  // compute Gamma * sma

  _spinmatrix_g0(smb[0], sma);
  _spinmatrix_g1(smb[1], sma);
  _spinmatrix_g2(smb[2], sma);
  _spinmatrix_g3(smb[3], sma);
  _spinmatrix_g5(smb[4], sma);
  _spinmatrix_g0g5(smb[5], sma);
  _spinmatrix_g5g0(smb[6], sma);
  _spinmatrix_g5g1(smb[7], sma);
  _spinmatrix_g5g2(smb[8], sma);
  _spinmatrix_g5g3(smb[9], sma);
  _spinmatrix_g0g1(smb[10], sma);
  _spinmatrix_g0g2(smb[11], sma);
  _spinmatrix_g0g3(smb[12], sma);
  _spinmatrix_g5g0g1(smb[13], sma);
  _spinmatrix_g5g0g2(smb[14], sma);
  _spinmatrix_g5g0g3(smb[15], sma);

  for (i = 0; i < 16; i++)
  {
    mult_mat_4spinor(smc, g[i], in);
    ctest = spin_matrix_norm_diff(smb[i], smc);
    if (creal(ctest) > 1e-15 || cimag(ctest) > 1e-15)
    {
      lprintf("MAIN", 0, "ERROR! _spinmatrix_%s mismatch, ctest=(%e,%e) should be 0 !\n", list_g[i], creal(ctest), cimag(ctest));
      return_value += 1;
    }
  }

  _g0_spinmatrix(smb[0], sma);
  _g1_spinmatrix(smb[1], sma);
  _g2_spinmatrix(smb[2], sma);
  _g3_spinmatrix(smb[3], sma);
  _g5_spinmatrix(smb[4], sma);
  _g0g5_spinmatrix(smb[5], sma);
  _g5g0_spinmatrix(smb[6], sma);
  _g5g1_spinmatrix(smb[7], sma);
  _g5g2_spinmatrix(smb[8], sma);
  _g5g3_spinmatrix(smb[9], sma);
  _g0g1_spinmatrix(smb[10], sma);
  _g0g2_spinmatrix(smb[11], sma);
  _g0g3_spinmatrix(smb[12], sma);
  _g5g0g1_spinmatrix(smb[13], sma);
  _g5g0g2_spinmatrix(smb[14], sma);
  _g5g0g3_spinmatrix(smb[15], sma);

  for (i = 0; i < 16; i++)
  {

    mult_4spinor_mat(smc, g[i], in);
    ctest = spin_matrix_norm_diff(smb[i], smc);
    if (creal(ctest) > 1e-15 || cimag(ctest) > 1e-15)
    {
      lprintf("MAIN", 0, "ERROR! _%s_spinmatrix mismatch, ctest=(%e,%e) should be 0 !\n", list_g[i], creal(ctest), cimag(ctest));
      return_value += 1;
    }
  }

  global_sum_int(&return_value, 1);
  lprintf("MAIN", 0, "return_value= %d\n ", return_value);
  finalize_process();

  return return_value;
}
