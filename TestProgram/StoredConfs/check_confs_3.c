/*******************************************************************************
*
* NOCOMPILE= !BC_T_PERIODIC
* NOCOMPILE= !NG=3
* NOCOMPILE= !REPR_FUNDAMENTAL
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
#include "wilsonflow.h"
#include "data_storage.h"

double openQCDWFobsl0[3] = {1.472098e+05, 8.273978e+03, -1.37e-01};

double openQCDWFobsl2[3] = {9.293370e+04, 1.122513e+04, 4.21e-01};

typedef struct _input_inverter
{
  double precision;
  double mass;
  double csw;
  double beta;
} input_inverter;

typedef struct _input_WF_meas
{
  double tmax;
  int nmeas;
  double eps;
  double delta;
  WF_integrator_type ittype;
} input_WF_meas;

int main(int argc, char *argv[])
{
  int return_value = 0;
  double test = 0.;
  struct timeval start, end, etime;
  data_storage_array *store;
  int idx[3];

  char cnfg_filename[256] = "cnfg/HiRep_qcd1_pb_L8T8_b6.0_c1.234_k0.13_r0_id4n1";

  setup_process(&argc, &argv);
  setup_gauge_fields();

  input_inverter INV_var = {.precision = 1e-16, .beta = 12.0, .mass = -0.14800000308159977, .csw = 1.13295};
  input_WF_meas WF_var = {.tmax = 0.2, .nmeas = 1, .eps = .8e-5, .delta = 1.0e-5, .ittype = RK3_ADAPTIVE};
#if defined(WITH_CLOVER) || defined(WITH_EXPCLOVER)
  set_csw(&(INV_var.csw));
#endif

  lprintf("MAIN", 0, "Inverter precision = %e\n", INV_var.precision);
  lprintf("MAIN", 0, "Mass = %f\n", INV_var.mass);
  lprintf("MAIN", 0, "beta = %.8f ct = %.8f\n", INV_var.beta, 1.0);

  lprintf("MAIN", 0, "WF tmax: %e\n", WF_var.tmax);
  lprintf("MAIN", 0, "WF number of measures: %d\n", WF_var.nmeas);
  lprintf("MAIN", 0, "WF initial epsilon: %e\n", WF_var.eps);
  lprintf("MAIN", 0, "WF delta: %e\n", WF_var.delta);
  lprintf("MAIN", 0, "WF integrator type: %d (0=Euler 1=3rd order Runge-Kutta 2=Adaptive 3rd order Runge-Kutta)\n", WF_var.ittype);

  /* initialize boundary conditions */
  BCs_pars_t BCs_pars = {.fermion_twisting_theta = {0., 0., 0., 0.}, .gauge_boundary_improvement_cs = 1.0, .gauge_boundary_improvement_ct = 1.0, .chiSF_boundary_improvement_ds = 1.0, .SF_BCs = 1};

  init_BCs(&BCs_pars);
  WF_initialize();

  lprintf("MAIN", 0, "Configuration from %s\n", cnfg_filename);
  read_gauge_field(cnfg_filename);

  represent_gauge_field();

  gettimeofday(&start, 0);
  lprintf("MAIN", 0, "2pt pion function still to be added\n");

  gettimeofday(&end, 0);
  timeval_subtract(&etime, &end, &start);

  gettimeofday(&start, 0);

  store = WF_update_and_measure(WF_var.ittype, u_gauge, &(WF_var.tmax), &(WF_var.eps), &(WF_var.delta), WF_var.nmeas, STORE);
  gettimeofday(&end, 0);
  timeval_subtract(&etime, &end, &start);

  print_data_storage(store);

  lprintf("MAIN", 0, "WF Observables flowed and measured in [%ld sec %ld usec]\n", etime.tv_sec, etime.tv_usec);

  test = 0.0;

  double tavg[3];

  lprintf("MAIN", 0, "Comparing Wl(t=0) Wl(t=.2) Yl(t=0) Yl(t=.2) TC(t=0) TC(t=.2) to the openQCD results found in:\nComparisonLogs/openQCD_qcd1_pb_L8T8_b6.0_c1.234_k0.13_r0_id4.ms3.log\n");
  idx[0] = 0;
  idx[1] = 1;
  tavg[0] = *data_storage_element(store, 0, idx);
  idx[1] = 2;
  tavg[1] = *data_storage_element(store, 0, idx);
  idx[1] = 3;
  tavg[2] = *data_storage_element(store, 0, idx);

  test = fabs(1. - GLB_VOLUME * tavg[0] / openQCDWFobsl0[0]);

  lprintf("TEST", 0, "Wl(t=0) relative difference: %.2e \n(should be around 1*10^(-7) or so)\n\n", test);
  if (test > 1e-6 && PID == 0)
  {
    lprintf("TEST", 0, "Test failed\n");
    return_value += 1;
  }

  test = fabs(1. - GLB_VOLUME * tavg[1] / openQCDWFobsl0[1]);
  lprintf("TEST", 0, "Yl(t=0) relative difference: %.2e \n(should be around 1*10^(-7) or so)\n\n", test);
  if (test > 1e-6 && PID == 0)
  {
    lprintf("TEST", 0, "Test failed\n");
    return_value += 1;
  }

  test = fabs(1. - tavg[2] / openQCDWFobsl0[2]);
  lprintf("TEST", 0, "TC(t=0) relative difference: %.2e \n(should be around 1*10^(-3) or so)\n\n", test);
  if (test > 1e-2 && PID == 0)
  {
    lprintf("TEST", 0, "Test failed\n");
    return_value += 1;
  }

  idx[0] = 1;
  idx[1] = 1;
  tavg[0] = *data_storage_element(store, 0, idx);
  idx[1] = 2;
  tavg[1] = *data_storage_element(store, 0, idx);
  idx[1] = 3;
  tavg[2] = *data_storage_element(store, 0, idx);

  test = fabs(1. - GLB_VOLUME * tavg[0] / openQCDWFobsl2[0]);

  lprintf("TEST", 0, "Wl(t=.2) relative difference: %.2e \n(should be around 1*10^(-7) or so)\n\n", test);
  if (test > 1e-6 && PID == 0)
  {
    lprintf("TEST", 0, "Test failed\n");
    return_value += 1;
  }

  test = fabs(1. - GLB_VOLUME * tavg[1] / openQCDWFobsl2[1]);
  lprintf("TEST", 0, "Yl(t=.2) relative difference: %.2e \n(should be around 1*10^(-7) or so)\n\n", test);
  if (test > 1e-6 && PID == 0)
  {
    lprintf("TEST", 0, "Test failed\n");
    return_value += 1;
  }

  test = fabs(1. - tavg[2] / openQCDWFobsl2[2]);
  lprintf("TEST", 0, "TC(t=.2) relative difference: %.2e \n(should be around 1*10^(-3) or so)\n\n", test);
  if (test > 1e-2 && PID == 0)
  {
    lprintf("TEST", 0, "Test failed\n");
    return_value += 1;
  }

  free_data_storage(store);

  WF_free();

  finalize_process();

  return return_value;
}
