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

#if !(defined(BASIC_SF)) && !(defined(ROTATED_SF))
#error This main code works only if some SF boundary conditions are enabled
#endif
/* we need the beta for normalization */
typedef struct _input_sfc
{
  double precision;
  char mstring[1024];
  char configlist[256];
  double csw;
  double beta;
  double rho_s;
  double rho_t;

  /* for the reading function */
  input_record_t read[10];

} input_sfc;

#define init_input_sfc(varname)                                                       \
  {                                                                                   \
    .read = {                                                                         \
      {"beta", "sf:beta = %lf", DOUBLE_T, &(varname).beta},                           \
      {"precision", "sf:precision = %lf", DOUBLE_T, &(varname).precision},            \
      {"quark masses", "sf:mass = %s", STRING_T, (varname).mstring},                  \
      {"csw coefficient", "sf:csw = %lg", DOUBLE_T, &(varname).csw},                  \
      {"smearing space", "sf:rho_s = %lg", DOUBLE_T, &(varname).rho_s},               \
      {"smearing time", "sf:rho_t = %lg", DOUBLE_T, &(varname).rho_t},                \
      {"Configuration list:", "sf:configlist = %s", STRING_T, &(varname).configlist}, \
      {NULL, NULL, INT_T, NULL}                                                       \
    }                                                                                 \
  }

input_sfc sfc_var = init_input_sfc(sfc_var);

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
      {"SF_ds", "SF_ds = %lf", DOUBLE_T, &(varname).SF_ds},        \
      {"SF_ct", "SF_ct = %lf", DOUBLE_T, &(varname).SF_ct},        \
      {NULL, NULL, INT_T, NULL}                                    \
    }                                                              \
  }

input_bcpar bcpar_var = init_input_bcpar(bcpar_var);

char cnfg_filename[256] = "";
char list_filename[256] = "";
char input_filename[256] = "input_file";
char output_filename[256] = "sfcoupling.out";
enum
{
  UNKNOWN_CNFG,
  DYNAMICAL_CNFG,
  QUENCHED_CNFG
};

typedef struct
{
  char string[256];
  int t, x, y, z;
  int nc, nf;
  double b, m;
  int n;
  int type;
} filename_t;

int main(int argc, char *argv[])
{
  int i;
  FILE *list;
  int nm;
  double m[256];
  double gsf;

  /* setup process id and communications */

  setup_process(&argc, &argv);

  setup_gauge_fields();

  read_input(glb_var.read, get_input_filename());
  read_input(sfc_var.read, get_input_filename());
  read_input(bcpar_var.read, get_input_filename());

#ifdef WITH_CLOVER
  set_csw(sfc_var.csw);
#endif

  strcpy(list_filename, sfc_var.configlist);

  if (strcmp(list_filename, "") != 0)
  {
    error((list = fopen(list_filename, "r")) == NULL, 1, "main [mk_sfcoupling.c]",
          "Failed to open list file\n");
  }

  nm = 1;
  m[0] = atof(sfc_var.mstring);

  lprintf("MAIN", 0, "list_filename = %s\n", list_filename);
  lprintf("MAIN", 0, "Inverter precision = %e\n", sfc_var.precision);
  lprintf("MAIN", 0, "Mass[%d] = %f\n", 0, m[0]);

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
#ifdef ROTATED_SF
  BCs_pars.gauge_boundary_improvement_ct = bcpar_var.SF_ct;
  BCs_pars.chiSF_boundary_improvement_ds = bcpar_var.SF_ds;
  BCs_pars.SF_BCs = 1;
#endif
#ifdef BASIC_SF
  BCs_pars.SF_BCs = 1;
#endif
  init_BCs(&BCs_pars);

#ifdef ROTATED_SF
  lprintf("MAIN", 0, "beta = %.8f\n rotatedSF ds = %.8f\n rotatedSF ct = %.8f\n", sfc_var.beta, bcpar_var.SF_ds, bcpar_var.SF_ct);
#else
  lprintf("MAIN", 0, "beta = %.8f\n", sfc_var.beta);
#endif

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

    full_plaquette();

    gettimeofday(&start, 0);
    // perform SF measurements.
    gsf = SF_action(sfc_var.beta);
    lprintf("SF_action", 10, "gsf = %.10e\n", gsf);
    SF_PCAC_wall_corr(m[0], sfc_var.precision);

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
