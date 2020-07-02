/*******************************************************************************
*
* NOCOMPILE= !NG=3
* NOCOMPILE= !BC_T_OPEN
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

double openQCDWsl0[8] = {
    2.1909123450308125e+03,
    2.9722361042026669e+03,
    3.0259644616131136e+03,
    3.0367303385187884e+03,
    3.0300970965090023e+03,
    2.9649859783298034e+03,
    2.9359214984384339e+03,
    2.2052896618928303e+03};

double openQCDWsl2[8] = {
    2.6765320194413977e+02,
    4.3236116872750966e+02,
    5.0381778362571441e+02,
    5.4540524750120494e+02,
    5.3406739744445190e+02,
    4.9821098379050790e+02,
    4.4323964021188277e+02,
    2.6762499410912039e+02};

double openQCDYsl0[8] = {
    0.0000000000000000e+00,
    5.9997765048648591e+02,
    6.1282361881016141e+02,
    6.3095249267986230e+02,
    6.1822993203120677e+02,
    6.0704066139211864e+02,
    5.9840893154888931e+02,
    0.0000000000000000e+00};

double openQCDYsl2[8] = {
    0.0000000000000000e+00,
    2.0686068189329623e+02,
    2.4733137411699627e+02,
    2.7787571681954273e+02,
    2.7083362437002722e+02,
    2.4972541162443318e+02,
    2.1416235919599703e+02,
    0.0000000000000000e+00};

double openQCDWFobsl0[3] = {2.2362137484535451e+04, 3.6674332869487243e+03, 3.2066549926534210e-01};

double openQCDWFobsl2[3] = {3.4923804173545318e+03, 1.4667891680202927e+03, 5.5754965636342813e-01};

typedef struct _input_obc
{
  double precision;
  double mass;
  double beta;
} input_obc;

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

  char cnfg_filename[256] = "cnfg/HiRep_obc_openQCD_qcd1_ob_L8T8_b12.0_c1.1329_k0.07n1";

  setup_process(&argc, &argv);
  setup_gauge_fields();

  input_obc OB_var = {.precision = 1e-16, .beta = 12.0, .mass = 0.1};
  input_WF_meas WF_var = {.tmax = 0.2, .nmeas = 1, .eps = 8.0e-05, .delta = .50e-6, .ittype = RK3};

  lprintf("MAIN", 0, "Inverter precision = %e\n", OB_var.precision);
  lprintf("MAIN", 0, "Mass = %f\n", OB_var.mass);
  lprintf("MAIN", 0, "beta = %.8f ct = %.8f\n", OB_var.beta, 1.0);

  lprintf("MAIN", 0, "WF tmax: %e\n", WF_var.tmax);
  lprintf("MAIN", 0, "WF number of measures: %d\n", WF_var.nmeas);
  lprintf("MAIN", 0, "WF initial epsilon: %e\n", WF_var.eps);
  lprintf("MAIN", 0, "WF delta: %e\n", WF_var.delta);
  lprintf("MAIN", 0, "WF integrator type: %d (0=Euler 1=3rd order Runge-Kutta 2=Adaptive 3rd order Runge-Kutta)\n", WF_var.ittype);

  /* initialize boundary conditions */
  BCs_pars_t BCs_pars = {.fermion_twisting_theta = {0., 0., 0., 0.}, .gauge_boundary_improvement_cs = 1.0, .gauge_boundary_improvement_ct = 1.0, .chiSF_boundary_improvement_ds = 1.0, .SF_BCs = 0};

  init_BCs(&BCs_pars);
  WF_initialize();

  lprintf("MAIN", 0, "Configuration from %s\n", cnfg_filename);
  read_gauge_field(cnfg_filename);

  represent_gauge_field();

  gettimeofday(&start, 0);

  store = WF_update_and_measure(WF_var.ittype, u_gauge, &(WF_var.tmax), &(WF_var.eps), &(WF_var.delta), WF_var.nmeas, STORE);
  gettimeofday(&end, 0);
  timeval_subtract(&etime, &end, &start);

  lprintf("MAIN", 0, "WF Observables flowed and measured in [%ld sec %ld usec]\n", etime.tv_sec, etime.tv_usec);

  double comb[2];

  for (int i = 0; i < GLB_T; i++)
  {
    if (i == 0)
    {
      idx[0] = 0;
      idx[1] = i;
      idx[2] = 1;
      comb[0] = *data_storage_element(store, 0, idx);
      idx[2] = 2;
      comb[0] += 2 * (*data_storage_element(store, 0, idx));

      idx[0] = 1;
      idx[1] = i;
      idx[2] = 1;
      comb[1] = *data_storage_element(store, 0, idx);
      idx[2] = 2;
      comb[1] += 2 * (*data_storage_element(store, 0, idx));
    }
    else
    {
      idx[0] = 0;
      idx[1] = i - 1;
      idx[2] = 1;
      comb[0] = *data_storage_element(store, 0, idx);
      idx[1] = i;
      comb[0] += *data_storage_element(store, 0, idx);
      idx[2] = 2;
      comb[0] += 2 * (*data_storage_element(store, 0, idx));

      idx[0] = 1;
      idx[1] = i - 1;
      idx[2] = 1;
      comb[1] = *data_storage_element(store, 0, idx);
      idx[1] = i;
      comb[1] += *data_storage_element(store, 0, idx);
      idx[2] = 2;
      comb[1] += 2 * (*data_storage_element(store, 0, idx));
    }

    comb[0] *= GLB_VOL3 * NG;
    comb[1] *= GLB_VOL3 * NG;
    if (fabs(comb[0]) > 1.e-16)
      test += fabs((comb[0] - openQCDWsl0[i]) / comb[0]);
    if (fabs(comb[1]) > 1.e-16)
      test += fabs((comb[1] - openQCDWsl2[i]) / comb[1]);

    idx[0] = 0;
    idx[1] = i;
    idx[2] = 3;
    comb[0] = *data_storage_element(store, 0, idx);
    idx[2] = 4;
    comb[0] += *data_storage_element(store, 0, idx);

    idx[0] = 1;
    idx[1] = i;
    idx[2] = 3;
    comb[1] = *data_storage_element(store, 0, idx);
    idx[2] = 4;
    comb[1] += *data_storage_element(store, 0, idx);

    comb[0] *= 2 * GLB_VOL3 * NG;
    comb[1] *= 2 * GLB_VOL3 * NG;

    if (fabs(comb[0]) > 1.e-16)
      test += fabs((comb[0] - openQCDYsl0[i]) / comb[0]);
    if (fabs(comb[1]) > 1.e-16)
      test += fabs((comb[1] - openQCDYsl2[i]) / comb[1]);

  }

  test /= 4. * GLB_T;

  lprintf("MAIN", 0, "Comparing Wl(t=0,T) Wl(t=.2,T) Yl(t=0,T) Yl(t=.2,T) to the openQCD results found in:\nComparisonLogs/openQCD_qcd1_sf_L8T16_b12.0_c1.13295_k0.1298027_r0_id5.ms3.log\n");

  lprintf("TEST", 0, "Cumulative relative difference: %.2e \n(should be around 1*10^(-12) or so)\n\n", test);
  if (test > 1e-11 && PID == 0)
  {
    lprintf("TEST", 0, "Test failed\n");
    return_value += 1;
  }

  double tavg[6];

  lprintf("MAIN", 0, "Comparing Wl(t=0) Wl(t=.2) Yl(t=0) Yl(t=.2) TC(t=0) TC(t=.2) to the openQCD results found in:\nComparisonLogs/openQCD_qcd1_sf_L8T16_b12.0_c1.13295_k0.1298027_r0_id5.ms3.log\n");
  idx[0] = 0;
  idx[1] = 0;
  tavg[0] = *data_storage_element(store, 1, idx);
  idx[1] = 1;
  tavg[1] = *data_storage_element(store, 1, idx);
  idx[1] = 2;
  tavg[2] = *data_storage_element(store, 1, idx);
  idx[1] = 3;
  tavg[3] = *data_storage_element(store, 1, idx);
  idx[1] = 4;
  tavg[4] = *data_storage_element(store, 1, idx);
  idx[1] = 5;
  tavg[5] = *data_storage_element(store, 1, idx);

  test = fabs(1. - 2.0 * GLB_VOLUME * NG * (tavg[1] + tavg[2]) / openQCDWFobsl0[0]);
  lprintf("TEST", 0, "Wl(t=0) relative difference: %.2e \n(should be around 1*10^(-13) or so)\n\n", test);

  if (test > 1e-11 && PID == 0)
  {
    lprintf("TEST", 0, "Test failed\n");
    return_value += 1;
  }

  test = fabs(1.0 - 2.0 * GLB_VOLUME * NG * (tavg[3] + tavg[4]) / openQCDWFobsl0[1]);
  lprintf("TEST", 0, "Yl(t=0) relative difference: %.2e \n(should be around 1*10^(-13) or so)\n\n", test);
  if (test > 1e-11 && PID == 0)
  {
    lprintf("TEST", 0, "Test failed\n");
    return_value += 1;
  }

  test = fabs(1.0 - tavg[5] / openQCDWFobsl0[2]);
  lprintf("TEST", 0, "TC(t=0) relative difference: %.2e \n(should be around 1*10^(-11) or so)\n\n", test);
  if (test > 1e-10 && PID == 0)
  {
    lprintf("TEST", 0, "Test failed\n");
    return_value += 1;
  }

  idx[0] = 1;
  idx[1] = 0;
  tavg[0] = *data_storage_element(store, 1, idx);
  idx[1] = 1;
  tavg[1] = *data_storage_element(store, 1, idx);
  idx[1] = 2;
  tavg[2] = *data_storage_element(store, 1, idx);
  idx[1] = 3;
  tavg[3] = *data_storage_element(store, 1, idx);
  idx[1] = 4;
  tavg[4] = *data_storage_element(store, 1, idx);
  idx[1] = 5;
  tavg[5] = *data_storage_element(store, 1, idx);

  test = fabs(1. - 2.0 * GLB_VOLUME * NG * (tavg[1] + tavg[2]) / openQCDWFobsl2[0]);
  lprintf("TEST", 0, "Wl(t=0.2) relative difference: %.2e \n(should be around 1*10^(-13) or so)\n\n", test);
  if (test > 1e-11 && PID == 0)
  {
    lprintf("TEST", 0, "Test failed\n");
    return_value += 1;
  }

  test = fabs(1.0 - 2.0 * GLB_VOLUME * NG * (tavg[3] + tavg[4]) / openQCDWFobsl2[1]);
  lprintf("TEST", 0, "Yl(t=0.2) relative difference: %.2e \n(should be around 1*10^(-13) or so)\n\n", test);
  if (test > 1e-11 && PID == 0)
  {
    lprintf("TEST", 0, "Test failed\n");
    return_value += 1;
  }

  test = fabs(1.0 - tavg[5] / openQCDWFobsl2[2]);
  lprintf("TEST", 0, "TC(t=0.2) relative difference: %.2e \n(should be around 1*10^(-13) or so)\n\n", test);
  if (test > 1e-11 && PID == 0)
  {
    lprintf("TEST", 0, "Test failed\n");
    return_value += 1;
  }

  free_data_storage(store);

  WF_free();

  finalize_process();

  return return_value;
}
