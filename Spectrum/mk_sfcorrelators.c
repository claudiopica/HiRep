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
#include "clover_tools.h"

typedef struct _input_sfc
{
  double precision;
  double mass;
  char configlist[256];
  char background[256];
  double csw;
  double beta;
  double rho_s;
  double rho_t;

  /* for the reading function */
  input_record_t read[11];

} input_sfc;

#define init_input_sfc(varname)                                                      \
  {                                                                                  \
    .read = {                                                                        \
      {"beta", "sf:beta = %lf", DOUBLE_T, &(varname).beta},                          \
      {"precision", "sf:precision = %lf", DOUBLE_T, &(varname).precision},           \
      {"quark mass", "sf:mass = %lf", DOUBLE_T, &(varname).mass},                    \
      {"csw coefficient", "sf:csw = %lf", DOUBLE_T, &(varname).csw},                 \
      {"smearing space", "sf:rho_s = %lf", DOUBLE_T, &(varname).rho_s},              \
      {"smearing time", "sf:rho_t = %lf", DOUBLE_T, &(varname).rho_t},               \
      {"SF background", "sf:background = %s", STRING_T, &(varname).background},      \
      {"Configuration list", "sf:configlist = %s", STRING_T, &(varname).configlist}, \
      {NULL, NULL, INT_T, NULL}                                                      \
    }                                                                                \
  }

input_sfc SF_var = init_input_sfc(SF_var);

/* BC variables */
typedef struct _input_bcpar
{
  /* rhmc parameters */
  double theta[4];
  double SF_ct;
  double SF_ds;
  /* for the reading function */
  input_record_t read[7];

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
      {NULL, NULL, INT_T, NULL}                                    \
    }                                                              \
  }

input_bcpar bcpar_var = init_input_bcpar(bcpar_var);

int main(int argc, char *argv[])
{
#if !(defined(BASIC_SF)) && !(defined(ROTATED_SF))
  error(1 == 1, 0, "main" __FILE__, "This code is to be used only if some SF BC are defined\n");
#endif

  int i;
  FILE *list = NULL;
  double gsf;
  char cnfg_filename[256] = "";

  setup_process(&argc, &argv);
  setup_gauge_fields();

  read_input(SF_var.read, get_input_filename());
  read_input(bcpar_var.read, get_input_filename());

#if defined(WITH_CLOVER) || defined(WITH_EXPCLOVER)
  set_csw(&SF_var.csw);
#endif

  if (strcmp(SF_var.configlist, "") != 0)
  {
    error((list = fopen(SF_var.configlist, "r")) == NULL, 1, "main [mk_sfcoupling.c]",
          "Failed to open list file\n");
  }

  lprintf("MAIN", 0, "Config list filename = %s\n", SF_var.configlist);
  lprintf("MAIN", 0, "Inverter precision = %e\n", SF_var.precision);
  lprintf("MAIN", 0, "Mass = %f\n", SF_var.mass);
#ifdef ROTATED_SF
  lprintf("MAIN", 0, "beta = %.8f\n rotatedSF ds = %.8f\n rotatedSF ct = %.8f\n", SF_var.beta, bcpar_var.SF_ds, bcpar_var.SF_ct);
#else
  lprintf("MAIN", 0, "beta = %.8f ct = %.8f\n", SF_var.beta, bcpar_var.SF_ct);
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
  BCs_pars.gauge_boundary_improvement_ct = bcpar_var.SF_ct;
#ifdef ROTATED_SF
  BCs_pars.chiSF_boundary_improvement_ds = bcpar_var.SF_ds;
#endif
  if (strcmp(SF_var.background, "true") == 0)
    BCs_pars.SF_BCs = 1;

  init_BCs(&BCs_pars);

  i = 0;
  while (++i)
  {
    struct timeval start, end, etime;

    if (list != NULL)
      if (fscanf(list, "%s", cnfg_filename) == 0 || feof(list))
        break;

    if (strcmp(cnfg_filename, "classical") == 0)
    {
      lprintf("MAIN", 0, "Generating a classical solution interpolating the SF boundaries\n");
      fflush(stdout);
      SF_classical_solution();
    }
    else
    {
      lprintf("MAIN", 0, "Configuration from %s\n", cnfg_filename);
      read_gauge_field(cnfg_filename);
    }

    represent_gauge_field();

    gettimeofday(&start, 0);
    // perform SF measurements.
    gsf = SF_action(SF_var.beta);
    lprintf("SF_action", 10, "gsf = %.10e\n", gsf);

    SF_PCAC_wall_corr(SF_var.mass, SF_var.precision, DONTSTORE);

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
