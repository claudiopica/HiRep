/*******************************************************************************
*
* Computation of the SF coupling and observables
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
#include "communications.h"
#include "setup.h"

/* Static Observables variables */
typedef struct _input_static_os_par
{
  /* observable parameters */
  char polymake[256];
  char wilsonmake[256];
  char configlist[256];

  /* for the reading function */
  input_record_t read[4];
} input_static_os_par;

#define init_input_ospar(varname)                                                        \
  {                                                                                      \
    .read = {                                                                            \
      {"Polyakov make", "so:Polyakov:make = %s", STRING_T, &(varname).polymake},         \
      {"Wilson looop make", "so:WilsonLoop:make = %s", STRING_T, &(varname).wilsonmake}, \
      {"Configuration list", "so:configlist = %s", STRING_T, &(varname).configlist},     \
      {NULL, NULL, INT_T, NULL}                                                          \
    }                                                                                    \
  }

input_static_os_par ospar_var = init_input_ospar(ospar_var);

typedef struct _input_HYP
{
  /*  int nsteps;*/
  double weight[3];

  /* for the reading function */
  input_record_t read[4];

} input_HYP;

#define init_input_HYP(varname)                                                          \
  {                                                                                      \
    .read = {                                                                            \
      {"HYP smearing weight[0]", "HYP:weight0 = %lf", DOUBLE_T, &((varname).weight[0])}, \
      {"HYP smearing weight[1]", "HYP:weight1 = %lf", DOUBLE_T, &((varname).weight[1])}, \
      {"HYP smearing weight[2]", "HYP:weight2 = %lf", DOUBLE_T, &((varname).weight[2])}, \
      {NULL, NULL, INT_T, NULL}                                                          \
    }                                                                                    \
  }

typedef struct _input_wilson
{
  int c[10][3];
  int nsteps[10];

  /* for the reading function */
  input_record_t read[41];

} input_WL;

#define init_input_WL(varname)                                                         \
  {                                                                                    \
    .read = {                                                                          \
      {"WL load path[0] delta(x)", "WL[0]:delta.x = %d", INT_T, &((varname).c[0][0])}, \
      {"WL load path[0] delta(y)", "WL[0]:delta.y = %d", INT_T, &((varname).c[0][1])}, \
      {"WL load path[0] delta(z)", "WL[0]:delta.z = %d", INT_T, &((varname).c[0][2])}, \
      {"WL load path[0] nsteps", "WL[0]:nsteps = %d", INT_T, &((varname).nsteps[0])},  \
      {"WL load path[1] delta(x)", "WL[1]:delta.x = %d", INT_T, &((varname).c[1][0])}, \
      {"WL load path[1] delta(y)", "WL[1]:delta.y = %d", INT_T, &((varname).c[1][1])}, \
      {"WL load path[1] delta(z)", "WL[1]:delta.z = %d", INT_T, &((varname).c[1][2])}, \
      {"WL load path[1] nsteps", "WL[1]:nsteps = %d", INT_T, &((varname).nsteps[1])},  \
      {"WL load path[2] delta(x)", "WL[2]:delta.x = %d", INT_T, &((varname).c[2][0])}, \
      {"WL load path[2] delta(y)", "WL[2]:delta.y = %d", INT_T, &((varname).c[2][1])}, \
      {"WL load path[2] delta(z)", "WL[2]:delta.z = %d", INT_T, &((varname).c[2][2])}, \
      {"WL load path[2] nsteps", "WL[2]:nsteps = %d", INT_T, &((varname).nsteps[2])},  \
      {"WL load path[3] delta(x)", "WL[3]:delta.x = %d", INT_T, &((varname).c[3][0])}, \
      {"WL load path[3] delta(y)", "WL[3]:delta.y = %d", INT_T, &((varname).c[3][1])}, \
      {"WL load path[3] delta(z)", "WL[3]:delta.z = %d", INT_T, &((varname).c[3][2])}, \
      {"WL load path[3] nsteps", "WL[3]:nsteps = %d", INT_T, &((varname).nsteps[3])},  \
      {"WL load path[4] delta(x)", "WL[4]:delta.x = %d", INT_T, &((varname).c[4][0])}, \
      {"WL load path[4] delta(y)", "WL[4]:delta.y = %d", INT_T, &((varname).c[4][1])}, \
      {"WL load path[4] delta(z)", "WL[4]:delta.z = %d", INT_T, &((varname).c[4][2])}, \
      {"WL load path[4] nsteps", "WL[4]:nsteps = %d", INT_T, &((varname).nsteps[4])},  \
      {"WL load path[5] delta(x)", "WL[5]:delta.x = %d", INT_T, &((varname).c[5][0])}, \
      {"WL load path[5] delta(y)", "WL[5]:delta.y = %d", INT_T, &((varname).c[5][1])}, \
      {"WL load path[5] delta(z)", "WL[5]:delta.z = %d", INT_T, &((varname).c[5][2])}, \
      {"WL load path[5] nsteps", "WL[5]:nsteps = %d", INT_T, &((varname).nsteps[5])},  \
      {"WL load path[6] delta(x)", "WL[6]:delta.x = %d", INT_T, &((varname).c[6][0])}, \
      {"WL load path[6] delta(y)", "WL[6]:delta.y = %d", INT_T, &((varname).c[6][1])}, \
      {"WL load path[6] delta(z)", "WL[6]:delta.z = %d", INT_T, &((varname).c[6][2])}, \
      {"WL load path[6] nsteps", "WL[6]:nsteps = %d", INT_T, &((varname).nsteps[6])},  \
      {"WL load path[7] delta(x)", "WL[7]:delta.x = %d", INT_T, &((varname).c[7][0])}, \
      {"WL load path[7] delta(y)", "WL[7]:delta.y = %d", INT_T, &((varname).c[7][1])}, \
      {"WL load path[7] delta(z)", "WL[7]:delta.z = %d", INT_T, &((varname).c[7][2])}, \
      {"WL load path[7] nsteps", "WL[7]:nsteps = %d", INT_T, &((varname).nsteps[7])},  \
      {"WL load path[8] delta(x)", "WL[8]:delta.x = %d", INT_T, &((varname).c[8][0])}, \
      {"WL load path[8] delta(y)", "WL[8]:delta.y = %d", INT_T, &((varname).c[8][1])}, \
      {"WL load path[8] delta(z)", "WL[8]:delta.z = %d", INT_T, &((varname).c[8][2])}, \
      {"WL load path[8] nsteps", "WL[8]:nsteps = %d", INT_T, &((varname).nsteps[8])},  \
      {"WL load path[9] delta(x)", "WL[9]:delta.x = %d", INT_T, &((varname).c[9][0])}, \
      {"WL load path[9] delta(y)", "WL[9]:delta.y = %d", INT_T, &((varname).c[9][1])}, \
      {"WL load path[9] delta(z)", "WL[9]:delta.z = %d", INT_T, &((varname).c[9][2])}, \
      {"WL load path[9] nsteps", "WL[9]:nsteps = %d", INT_T, &((varname).nsteps[9])},  \
      {NULL, NULL, INT_T, NULL}                                                        \
    }                                                                                  \
  }

input_HYP HYP_var = init_input_HYP(HYP_var);
input_WL WL_var = init_input_WL(WL_var);




/* BC variables */
typedef struct _input_bcpar
{
  /* rhmc parameters */
  double theta[4];
  double SF_ct;
  double SF_ds;
  char SF_bkg[256];
  /* for the reading function */
  input_record_t read[8];

} input_bcpar;

#define init_input_bcpar(varname)                                  \
{                                                                \
    .read = {                                                      \
      {"theta_T", "theta_T = %lf", DOUBLE_T, &(varname).theta[0]}, \
      {"theta_X", "theta_X = %lf", DOUBLE_T, &(varname).theta[1]}, \
      {"theta_Y", "theta_Y = %lf", DOUBLE_T, &(varname).theta[2]}, \
      {"theta_Z", "theta_Z = %lf", DOUBLE_T, &(varname).theta[3]}, \
      {"SF_ds", "sf:ds = %lf", DOUBLE_T, &(varname).SF_ds},        \
      {"SF_ct", "sf:ct = %lf", DOUBLE_T, &(varname).SF_ct},        \
      {"SF_background", "sf:background = %s", STRING_T, &(varname).SF_bkg},        \
      {NULL, NULL, INT_T, NULL}                                    \
    }                                                              \
  }

input_bcpar bcpar_var = init_input_bcpar(bcpar_var);


int main(int argc, char *argv[])
{

  int i;
  FILE *list = NULL;

  char cnfg_filename[256] = "";

  setup_process(&argc, &argv);
  setup_gauge_fields();

  read_input(ospar_var.read, get_input_filename());
  read_input(bcpar_var.read, get_input_filename());

  HYP_var.weight[0] = HYP_var.weight[1] = HYP_var.weight[2] = 0.;
  for (i = 0; i < 10; i++)
  {
    WL_var.c[i][0] = WL_var.c[i][1] = WL_var.c[i][2] = WL_var.nsteps[i] = 0;
  }
  read_input(HYP_var.read, get_input_filename());
  read_input(WL_var.read, get_input_filename());

  if (strcmp(ospar_var.configlist, "") != 0)
  {
    error((list = fopen(ospar_var.configlist, "r")) == NULL, 1, "main [mk_static_observables.c]",
          "Failed to open list file\n");
  }

  lprintf("MAIN", 0, "Config list filename = %s\n", ospar_var.configlist);

  /* initialize boundary conditions */
  BCs_pars_t BCs_pars = {
      .fermion_twisting_theta = {0., 0., 0., 0.},
      .gauge_boundary_improvement_cs = 1.,
      .gauge_boundary_improvement_ct = 1.,
      .chiSF_boundary_improvement_ds = 1.,
      .SF_BCs = 0};
#ifdef FERMION_THETA
  BCs_pars.fermion_twisting_theta[0] = bcpar_var.theta[0];
  BCs_pars.fermion_twisting_theta[1] = bcpar_var.theta[1];
  BCs_pars.fermion_twisting_theta[2] = bcpar_var.theta[2];
  BCs_pars.fermion_twisting_theta[3] = bcpar_var.theta[3];
#endif
  BCs_pars.gauge_boundary_improvement_ct = bcpar_var.SF_ct;
#ifdef ROTATED_SF
  BCs_pars.chiSF_boundary_improvement_ds = bcpar_var.SF_ds;
#endif
  if (strcmp(bcpar_var.SF_bkg, "true") == 0)
    BCs_pars.SF_BCs = 1;


  init_BCs(&BCs_pars);

  if (strcmp(ospar_var.wilsonmake, "true") == 0)
  {

#ifndef BC_T_PERIODIC
    error(0 == 0, 1, "mk_static_observables.c", "Wilson Loop can bea measured only for T periodic BC");
#endif

    lprintf("MAIN", 0, "HYP smearing weights: %f %f %f\n", HYP_var.weight[0], HYP_var.weight[1], HYP_var.weight[2]);

    WL_initialize();
    for (i = 0; i < 10; i++)
    {
      if (WL_var.c[i][0] * WL_var.c[i][0] + WL_var.c[i][1] * WL_var.c[i][1] + WL_var.c[i][2] * WL_var.c[i][2] != 0 && WL_var.nsteps[i] != 0)
        WL_load_path(WL_var.c[i], WL_var.nsteps[i]);
    }
  }

  i = 0;
  while (++i)
  {
    struct timeval start, end, etime;

    if (list != NULL)
      if (fscanf(list, "%s", cnfg_filename) == 0 || feof(list))
        break;

    lprintf("MAIN", 0, "Configuration from %s\n", cnfg_filename);
    read_gauge_field(cnfg_filename);

    represent_gauge_field();

    gettimeofday(&start, 0);

    if (strcmp(ospar_var.polymake, "true") == 0)
      polyakov();

    if (strcmp(ospar_var.wilsonmake, "true") == 0)
      WL_wilsonloops(HYP_var.weight);

    gettimeofday(&end, 0);
    timeval_subtract(&etime, &end, &start);
    lprintf("MAIN", 0, "Configuration #%d: analysed in [%ld sec %ld usec]\n", i, etime.tv_sec, etime.tv_usec);

    if (list == NULL)
      break;
  }

  if (list != NULL)
    fclose(list);

  finalize_process();

  return 0;
}
