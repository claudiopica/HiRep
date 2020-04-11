/*******************************************************************************
*
* Computation of the observable E(t) evolved with the Wilson flow
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

typedef struct _input_WF_meas
{
  double tmax;
  int nmeas;
  double eps;
  double delta;
  char configlist[256];

  /* for the reading function */
  input_record_t read[6];

} input_WF_meas;

#define init_input_WF_meas(varname)                                                  \
  {                                                                                  \
    .read = {                                                                        \
      {"WF integration time", "WF:tmax = %lf", DOUBLE_T, &((varname).tmax)},         \
      {"Configuration list", "WF:configlist = %s", STRING_T, &(varname).configlist}, \
      {"WF number of measures", "WF:nmeas = %d", DOUBLE_T, &((varname).nmeas)},      \
      {"WF initial epsilon", "WF:eps = %lf", DOUBLE_T, &((varname).eps)},            \
      {"WF delta", "WF:delta = %lf", DOUBLE_T, &((varname).delta)},                  \
      {NULL, NULL, 0, NULL}                                                          \
    }                                                                                \
  }

input_WF_meas WF_var = init_input_WF_meas(WF_var);

#if defined(BASIC_SF)

typedef struct _input_SF
{
  double ct;
  double beta;
  int background;

  /* for the reading function */
  input_record_t read[4];

} input_SF;

#define init_input_SF(varname)                                                 \
  {                                                                            \
    .read = {                                                                  \
      {"SF background", "SF:background = %d", INT_T, &((varname).background)}, \
      {"SF ct", "SF:ct = %lf", DOUBLE_T, &((varname).ct)},                     \
      {"SF beta", "SF:beta = %lf", DOUBLE_T, &((varname).beta)},               \
      {NULL, NULL, 0, NULL}                                                    \
    }                                                                          \
  }

input_SF SF_var = init_input_SF(SF_var);

#endif
#if defined(ROTATED_SF)

typedef struct _input_SF
{
  double ct;
  double beta;
  double zf = 1.;
  double ds = 1.;
  int sign = 1;
  int background;

  /* for the reading function */
  input_record_t read[7];

} input_SF;

#define init_input_SF(varname)                                                 \
  {                                                                            \
    .read = {                                                                  \
      {"SF background", "SF:background = %d", INT_T, &((varname).background)}, \
      {"SF sign", "SF:sign = %d", INT_T, &((varname).sign)},                   \
      {"SF ct", "SF:ct = %lf", DOUBLE_T, &((varname).ct)},                     \
      {"SF zf", "SF:zf = %lf", DOUBLE_T, &((varname).zf)},                     \
      {"SF ds", "SF:ds = %lf", DOUBLE_T, &((varname).ds)},                     \
      {"SF beta", "SF:beta = %lf", DOUBLE_T, &((varname).beta)},               \
      {NULL, NULL, 0, NULL}                                                    \
    }                                                                          \
  }

input_SF SF_var = init_input_SF(SF_var);

#endif

typedef struct _input_bcpar
{
  /* rhmc parameters */
  double theta[4];
  /* for the reading function */
  input_record_t read[5];

} input_bcpar;

#define init_input_bcpar(varname)                                  \
  {                                                                \
    .read = {                                                      \
      {"theta_T", "theta_T = %lf", DOUBLE_T, &(varname).theta[0]}, \
      {"theta_X", "theta_X = %lf", DOUBLE_T, &(varname).theta[1]}, \
      {"theta_Y", "theta_Y = %lf", DOUBLE_T, &(varname).theta[2]}, \
      {"theta_Z", "theta_Z = %lf", DOUBLE_T, &(varname).theta[3]}, \
      {NULL, NULL, INT_T, NULL}                                    \
    }                                                              \
  }

input_bcpar bcpar_var = init_input_bcpar(bcpar_var);

int main(int argc, char *argv[])
{
  int i;
  FILE *list = NULL;

  /* setup process id and communications */
  setup_process(&argc, &argv);

  setup_gauge_fields();

  read_input(WF_var.read, get_input_filename());

#if defined(ROTATED_SF) || defined(BASIC_SF)
  read_input(SF_var.read, get_input_filename());
#endif

  lprintf("MAIN", 0, "list file [%s]\n", WF_var.configlist);

  if (strcmp(WF_var.configlist, "") != 0)
  {
    error((list = fopen(WF_var.configlist, "r")) == NULL, 1, "main [WF_measure.c]",
          "Failed to open list file\n");
  }

  double dt = (double)WF_var.tmax / (double)WF_var.nmeas;

  lprintf("MAIN", 0, "WF tmax: %e\n", WF_var.tmax);
  lprintf("MAIN", 0, "WF number of measures: %d\n", WF_var.nmeas);
  lprintf("MAIN", 0, "WF initial epsilon: %e\n", WF_var.eps);
  lprintf("MAIN", 0, "WF delta: %e\n", WF_var.delta);
  lprintf("MAIN", 0, "WF measurement interval dt : %e\n", dt);

#ifdef ROTATED_SF
  lprintf("MAIN", 0, "beta = %.8f\n rotatedSF ds = %.8f\n rotatedSF ct = %.8f\n", SF_var.beta, SF_var.ds, SF_var.ct);
#elif defined(BASIC_SF)
  lprintf("MAIN", 0, "beta = %.8f ct = %.8f\n", SF_var.beta, SF_var.ct);
#endif

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
#if defined(ROTATED_SF) || defined(BASIC_SF)
  BCs_pars.gauge_boundary_improvement_ct = SF_var.ct;
  error(SF_var.background != 0 && SF_var.background != 1, 0, "init_mc_ghmc" __FILE__, "Wrong value of SF_background\n");
  BCs_pars.SF_BCs = SF_var.background;
#if defined(ROTATED_SF)
  BCs_pars.chiSF_boundary_improvement_ds = SF_var.ds;
#endif
#endif

  init_BCs(&BCs_pars);

  WF_initialize();

  char cnfg_filename[256];
  struct timeval start, end, etime;

  i = 0;
  while (1)
  {

    if (list != NULL)
      if (fscanf(list, "%s", cnfg_filename) == 0 || feof(list))
        break;

    i++;

    lprintf("MAIN", 0, "Configuration from %s\n", cnfg_filename);

    read_gauge_field(cnfg_filename);

    apply_BCs_on_fundamental_gauge_field();

    gettimeofday(&start, 0);

    full_plaquette();

    WF_adaptive_full_measure(u_gauge, &(WF_var.tmax), &(WF_var.eps), &(WF_var.delta), WF_var.nmeas);

    gettimeofday(&end, 0);
    timeval_subtract(&etime, &end, &start);
    lprintf("TIMING", 0, "Wilson Flow evolution and measurements for configuration  [%s] done [%ld sec %ld usec]\n", cnfg_filename, etime.tv_sec, etime.tv_usec);

    if (list == NULL)
      break;
  } // end loop configurations

  if (list != NULL)
    fclose(list);

  WF_free();

  finalize_process();

  return 0;
}
