
/******************************************************************************
 *
 *
 * File check_triplets_4.c
 *
 * Check of the triplet mesons (free case): semwall source case
 *
 * Author: Vincent Drach
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
#include "spectrum.h"
#include "clover_tools.h"
#include "communications.h"
#include "clover_tools.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288419716939937510
#endif

static double gid[16] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}; /* gid = tr \bar{Gamma} Gamma / 4 */
static double g0[16] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};  /* g0 = tr gamma_0 Gamma gamma_0 \bar{Gamma} / 4  Because our matrices Gamma are not necessary hermitiean !*/
static double g1[16] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};  /* g1 = tr gamma_1 \bar{Gamma} gamma_1 Gamma / 4 */
static double g2[16] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};  /* g2 = tr gamma_2 \bar{Gamma} gamma_2 Gamma / 4 */
static double g3[16] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};  /* g3 = tr gamma_3 \bar{Gamma} gamma_3 Gamma / 4 */
static double gb[16] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
char *mes_channel_names[16] = {"g5", "id", "g0", "g1", "g2", "g3", "g0g1", "g0g2", "g0g3", "g0g5", "g5g1", "g5g2", "g5g3", "g0g5g1", "g0g5g2", "g0g5g3"};

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

#define set_zero_mat(A)          \
  {                              \
    int _i, _j;                  \
    for (_i = 0; _i < 4; _i++)   \
      for (_j = 0; _j < 4; _j++) \
      {                          \
        A[_i][_j] = 0.;          \
      }                          \
  }
#define trace_mat(r, A)        \
  {                            \
    int _i;                    \
    _complex_0(r);             \
    for (_i = 0; _i < 4; _i++) \
    {                          \
      r += A[_i][_i];          \
    }                          \
  }

static void adj_mat(double complex At[4][4], double complex A[4][4])
{
  int i, j;
  for (i = 0; i < 4; i++)
  {
    for (j = 0; j < 4; j++)
    {
      At[i][j] = conj(A[j][i]);
    }
  }
}

/*VD: This is kept to to some debugging.
  static void print_mat(double complex mat[4][4], const char name[]) {
    int i,j;
    lprintf("MAIN",0,"%s = \n", name);
    for(i=0; i<4; i++) {
        lprintf("MAIN",0,"[ ");
          for(j=0; j<4; j++) {
            lprintf("MAIN",0,"(%.2f,%.2f) ",creal(mat[i][j]),cimag(mat[i][j]));
          }
          lprintf("MAIN",0,"]\n");
        }
      }*/
/* Mesons parameters */
typedef struct _input_mesons
{
  char mstring[256];

  /* for the reading function */
  input_record_t read[4];
  double csw;
  int nhits_2pt;

} input_mesons;

#define init_input_mesons(varname)                                                                        \
  {                                                                                                       \
    .read = {                                                                                             \
      {"fermion mass", "mes:mass = %s", STRING_T, (varname).mstring},                          \
      {"csw coefficient", "mes:csw = %lg", DOUBLE_T, &(varname).csw},                                     \
      {"number of noisy sources per cnfg for 2pt fn", "mes:nhits_2pt = %d", INT_T, &(varname).nhits_2pt}, \
      {NULL, NULL, INT_T, NULL}                                                                           \
    }                                                                                                     \
  }

static double mass;
void free_correlators(double **triplets);

input_glb glb_ip = init_input_glb(glb_ip);
input_mesons mes_ip = init_input_mesons(mes_ip);

/*\bar{Gamma} Gamma / 4 */
double complex get_gid(double complex Gamma[4][4])
{
  double complex tmp[4][4], tmp2[4][4];
  double complex Gammabar[4][4];
  int sign;
  double complex r;
  double complex g0[4][4];
  g0_debug(g0, &sign);

  adj_mat(tmp, Gamma);
  mult_mat(tmp2, g0, tmp);
  mult_mat(Gammabar, tmp2, g0);
  mult_mat(tmp, Gammabar, Gamma);
  trace_mat(r, tmp);
  return r / 4.0;
}
char char_t[100];
FILE *fp;
char path[1035];

/* Ugly: read the correlator from the output files ! This is because the function that compute the props, make the contractions and write the output zeros the correlators after writting them.
The purpose of this entire test is to check the function without modifying it */
static double read_output(int t, int i, int re_im_flag)
{
  char char_t[100];
  FILE *fp;
  char path[1035];
  char command[500];
  double out;

  sprintf(char_t, "%d", t + 1);
  strcpy(command, "grep \"TRIPLET ");
  strcat(command, mes_channel_names[i]);
  if (re_im_flag == 0) // get the real part.
  {
    strcat(command, "=\" out_0 | awk -F'=' '{print $3}' | awk '{for (i=1;i<=NF;i++) print $i}' | awk 'NR==");
  }
  if (re_im_flag == 1)
  {
    strcat(command, "_im=\" out_0 | awk -F'=' '{print $3}' | awk  '{for (i=1;i<=NF;i++) print $i}' | awk 'NR==");
  }
  strcat(command, char_t);
  strcat(command, "'");

  fp = popen(command, "r");
  while (fgets(path, sizeof(path) - 1, fp) != NULL)
  {
    sscanf(path, "%lf", &out);
  }

  /* close */
  pclose(fp);

  return out;
}

/* = tr gamma_mu \bar{Gamma} gamma_mu Gamma / 4 \bar{Gamma} = gamma_0 Gamma^\dagger gamma_0 */
double complex get_gmu(double complex Gamma[4][4], int mu)
{
  int sign;
  double complex tmp[4][4], tmp2[4][4], tmp3[4][4];
  double complex Gammabar[4][4];
  double complex r;
  double complex g0[4][4], gmu[4][4];

  g0_debug(g0, &sign);
  if (mu == 1)
    g1_debug(gmu, &sign);
  if (mu == 2)
    g2_debug(gmu, &sign);
  if (mu == 3)
    g3_debug(gmu, &sign);
  if (mu == 0)
    g0_debug(gmu, &sign);
  adj_mat(tmp, Gamma);

  mult_mat(tmp2, g0, tmp);
  mult_mat(Gammabar, tmp2, g0);

  mult_mat(tmp, gmu, Gammabar);
  mult_mat(tmp2, gmu, Gamma);
  mult_mat(tmp3, tmp, tmp2);

  trace_mat(r, tmp3);
  return r / 4.0;
}

double complex get_g0(double complex Gamma[4][4])
{
  double complex tmp[4][4];
  double complex r;
  set_zero_mat(tmp);
  mult_mat(tmp, Gamma, Gamma);
  trace_mat(r, tmp);
  return r / 4.0;
}
int main(int argc, char *argv[])
{
  int i, t, sign;
  double **ex_triplets;
  char pame[256];
  int return_value = 0;
  double tol=1.e-4;
  double complex g[16][4][4];
  g5_debug(g[0], &sign);
  id_debug(g[1], &sign);
  g0_debug(g[2], &sign);
  g1_debug(g[3], &sign);
  g2_debug(g[4], &sign);
  g3_debug(g[5], &sign);
  mult_mat(g[6], g[2], g[3]) //g0g1
  g0g2_debug(g[7], &sign);
  g0g3_debug(g[8], &sign);
  mult_mat(g[9], g[2], g[0]) //g0g5
  g5g1_debug(g[10], &sign);
  g5g2_debug(g[11], &sign);
  g5g3_debug(g[12], &sign);
  mult_mat(g[13], g[2], g[10]); //  g0g5g1
  mult_mat(g[14], g[2], g[11]); //  g0g5g2
  mult_mat(g[15], g[2], g[12]); //  g0g5g3

  char *mes_channel_names[16] = {"g5", "id", "g0", "g1", "g2", "g3", "g0g1", "g0g2", "g0g3", "g0g5", "g5g1", "g5g2", "g5g3", "g0g5g1", "g0g5g2", "g0g5g3"};

  for (i = 0; i < 16; i++)
  {
    gid[i] = get_gid(g[i]);
    g0[i] = get_gmu(g[i], 0);
    g1[i] = get_gmu(g[i], 1);
    g2[i] = get_gmu(g[i], 2);
    g3[i] = get_gmu(g[i], 3);
  }

  logger_map("DEBUG", "debug");
  logger_setlevel(0, 200);

  /* setup process id and communications */
  setup_process(&argc, &argv);

  setup_gauge_fields();

  read_input(mes_ip.read, get_input_filename());

  strcpy(pame, mes_ip.mstring);
  mass = atof(strtok(pame, ";"));
 

  lprintf("MAIN", 0, "mes:masses = %f\n", mass);
  lprintf("MAIN", 0, "mes:nhits_2pt = %d\n", mes_ip.nhits_2pt);

#if defined(WITH_CLOVER) || defined(WITH_EXPCLOVER)
  set_csw(&mes_ip.csw);
#endif
  unit_u(u_gauge);
  represent_gauge_field();
 #ifdef REPR_FUNDAMENTAL 
  apply_BCs_on_represented_gauge_field(); //This is a trick: the BCs are not applied in the case the REPR is fundamental because represent_gauge field assumes that the right BCs are already applied on the fundamental field!
  #endif
  ex_triplets = (double **) malloc(sizeof(double*)*16);
  for (i=0;i<16;i++) ex_triplets[i] = (double *) calloc(GLB_T,sizeof(double*));
  

  lprintf("MAIN", 0, "Measuring Gamma Gamma correlators and PCAC-mass\n");
  init_meson_correlators(0);
  lprintf("MAIN", 0, "Zerocoord{%d,%d,%d,%d}\n", zerocoord[0], zerocoord[1], zerocoord[2], zerocoord[3]);

  error(!(GLB_X == GLB_Y && GLB_X == GLB_Z), 1, "main", "This test works only for GLB_X=GLB_Y=GLB_Z");
  measure_spectrum_semwall(1, &mass, mes_ip.nhits_2pt, 0, 1e-28);

  double complex corr_triplets[16][GLB_T];

  for (i = 0; i < 16; i++)
  {
    for (t = 0; t < GLB_T; t++)
    {
      corr_triplets[i][t] = read_output(t, i, 0) + I * read_output(t, i, 1);
    }
  }


  //  /* CALCOLO ESPLICITO */
  free_correlators(ex_triplets);

  lprintf("TEST", 0, "\nANALITICO\tPOINT-TO-ALL\tERROR (must be less than %e)\n",tol);
  for (i = 0; i < 16; i++)
  {
    lprintf("TEST", 0, "TRIPLET CORRELATOR %s\n", mes_channel_names[i]);
    for (t = 0; t < GLB_T; t++)
    {
      lprintf("TEST", 0, "%e\t%e\t%e\n", ex_triplets[i][t], creal(corr_triplets[i][t]), fabs(ex_triplets[i][t] - creal(corr_triplets[i][t])));
      if (fabs(ex_triplets[i][t] - creal(corr_triplets[i][t])) > tol)
      {
        return_value += 1;
      }
      if (fabs(cimag(corr_triplets[i][t])) > tol)
      {
        return_value += 1;
      }
    }
  }
  
  global_sum_int(&return_value, 1);
  lprintf("MAIN", 0, "return_value= %d\n ", return_value);
  free_meson_observables();
  for (i=0;i<16;i++)     free(ex_triplets[i]);
  free(ex_triplets);

  finalize_process();
  return return_value;
}

double square(double x)
{
  return x * x;
}

double complex ev(double k[4])
{
  double complex z;
  z = mass + 4.0 + cos((2.0 * M_PI * k[0]) / GLB_T) + cos((2.0 * M_PI * k[1]) / GLB_X) + cos((2.0 * M_PI * k[2]) / GLB_Y) + cos((2.0 * M_PI * k[3]) / GLB_Z) + I * sqrt(square(sin((2.0 * M_PI * k[0]) / GLB_T)) + square(sin((2.0 * M_PI * k[1]) / GLB_X)) + square(sin((2.0 * M_PI * k[2]) / GLB_Y)) + square(sin((2.0 * M_PI * k[3]) / GLB_Z)));
  return z;
}

void free_correlators(double **triplets)
{
  double A2[GLB_T], B2[4][GLB_T];
  double complex A[GLB_T], B[4][GLB_T];
  double complex tmp, eit;
  double norm2, z;
  int i, j, t;
  double k[4];
  double sigma[4] = {0., 0., 0., 0.};
#ifdef BC_T_ANTIPERIODIC
  sigma[0] = .5;
#endif
#ifdef BC_X_ANTIPERIODIC
  sigma[1] = .5;
#endif
#ifdef BC_Y_ANTIPERIODIC
  sigma[2] = .5;
#endif
#ifdef BC_Z_ANTIPERIODIC
  sigma[3] = .5;
#endif

  z = -4. * NF / square(GLB_T * GLB_X * GLB_Y * GLB_Z);

  lprintf("FREE", 0, "sigma = (%f,%f,%f,%f)\n", sigma[0], sigma[1], sigma[2], sigma[3]);
  for (t = 0; t < GLB_T; t++)
  {
    A2[t] = 0.;
    B2[0][t] = B2[1][t] = B2[2][t] = B2[3][t] = 0.;
  }

  for (k[1] = sigma[1]; k[1] < GLB_X + sigma[1] - .5; k[1] += 1.)
    for (k[2] = sigma[2]; k[2] < GLB_Y + sigma[2] - .5; k[2] += 1.)
      for (k[3] = sigma[3]; k[3] < GLB_Z + sigma[3] - .5; k[3] += 1.)
      {

        for (t = 0; t < GLB_T; t++)
        {
          A[t] = 0.;
          B[0][t] = 0.;
          B[1][t] = 0.;
          B[2][t] = 0.;
          B[3][t] = 0.;
        }
        for (k[0] = sigma[0]; k[0] < GLB_T + sigma[0] - .5; k[0] += 1.)
        {
          tmp = ev(k);
          norm2 = tmp * conj(tmp);

          for (t = 0; t < GLB_T; t++)
          {
            eit = cos((2.0 * M_PI * t * k[0]) / GLB_T) + I * sin((2.0 * M_PI * t * k[0]) / GLB_T);
            A[t] += creal(tmp) * eit / norm2;
            B[0][t] += sin((2.0 * M_PI * k[0]) / GLB_T) * eit / norm2;
            for (j = 1; j < 4; j++)
              B[j][t] += sin((2.0 * M_PI * k[j]) / GLB_X) * eit / norm2;
          }
        }
        for (t = 0; t < GLB_T; t++)
        {
          A2[t] += A[t] * conj(A[t]) * z;
          for (j = 0; j < 4; j++)
            B2[j][t] += B[j][t] * conj(B[j][t]) * z;
        }
      }

  lprintf("FREE", 10, "A2\tB2[0]\tB2[1]\tB2[2]\tB[3]\n", A2[0]);
  for (t = 0; t < GLB_T; t++)
    lprintf("FREE", 10, "%e\t%e\t%e\t%e\t%e\n", A2[t], B2[0][t], B2[1][t], B2[2][t], B2[3][t]);

  for (i = 0; i < 16; i++)
  {
    for (t = 0; t < GLB_T; t++)
    {
      triplets[i][t] = gb[i] * (gid[i] * A2[t] - g0[i] * B2[0][t] - g1[i] * B2[1][t] - g2[i] * B2[2][t] - g3[i] * B2[3][t]);
    }
  }

  lprintf("FREE", 0, "Exact free correlators computed.\n");
}
