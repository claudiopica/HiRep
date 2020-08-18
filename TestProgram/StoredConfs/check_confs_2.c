/*******************************************************************************
*
* NOCOMPILE= !BASIC_SF
* NOCOMPILE= !NG=3
* NOCOMPILE= !WITH_CLOVER
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

double openQCDfA[16] = {-0.0000000000000000e+00,
                        -3.1886983381171885e-01,
                        -2.8168610838137947e-01,
                        -1.6296972021618447e-01,
                        -2.0464759666502405e-01,
                        -1.5302301558626705e-01,
                        -1.2811706222607241e-01,
                        -1.1176979625577449e-01,
                        -1.1257627355480730e-01,
                        -1.1229029627859480e-01,
                        -1.0876996633725135e-01,
                        -1.0334190170514863e-01,
                        -9.9203586601439922e-02,
                        -9.7345066304528730e-02,
                        -9.4401258139749991e-02,
                        -9.3003914262174150e-02};

double openQCDgA[16] = {0.0000000000000000e+00,
                        -9.6535436097137106e-02,
                        -9.9385252358823548e-02,
                        -1.0138126146229129e-01,
                        -1.0202098874661183e-01,
                        -1.0300082713203582e-01,
                        -1.0281174483401950e-01,
                        -1.0419570448052083e-01,
                        -1.0395126127179871e-01,
                        -1.0565530616260711e-01,
                        -1.0784956649237709e-01,
                        -1.0829355121977642e-01,
                        -1.1650170912005142e-01,
                        -1.2644063781661752e-01,
                        -1.4206418218034111e-01,
                        -1.9553748084642025e-01};

double openQCDfP[16] = {-0.0000000000000000e+00,
                        4.9759122473194157e+00,
                        4.4426544988476966e+00,
                        3.8442054011323230e+00,
                        3.2755695842044039e+00,
                        2.7293295903592725e+00,
                        2.1182614555819881e+00,
                        1.5932087534069062e+00,
                        1.1512953657110654e+00,
                        8.2832271565452364e-01,
                        6.2107095253505562e-01,
                        4.6306581656314150e-01,
                        3.2317110990214815e-01,
                        2.1531570936394953e-01,
                        1.4742159809300009e-01,
                        1.0859343752359941e-01};

double openQCDgP[16] = {-0.0000000000000000e+00,
                        1.0193950640377299e-01,
                        1.0935820823473641e-01,
                        1.1324091176924735e-01,
                        1.1576605463646100e-01,
                        1.1823203593167542e-01,
                        1.2048361693523496e-01,
                        1.2963488331221493e-01,
                        1.4871964015660472e-01,
                        1.8757147386545184e-01,
                        2.5144925452165340e-01,
                        3.6071789110330543e-01,
                        5.5487981824148258e-01,
                        9.2207434147757461e-01,
                        1.6139374768366910e+00,
                        2.9674665793728208e+00};

double openQCDf1 = 8.0350412186568132e-02;

double openQCDWsl0[16] = {6.8269748822088570e+02,
                          3.0783638808513065e+03,
                          3.2431322761917768e+03,
                          3.2481427049613180e+03,
                          3.2395148742307983e+03,
                          3.2156422642098655e+03,
                          3.2934475921063968e+03,
                          3.1966816829788331e+03,
                          3.1619703208773417e+03,
                          3.1608452484588292e+03,
                          3.1915878183812733e+03,
                          3.2088985812543629e+03,
                          3.2729546797836597e+03,
                          3.2208507878528194e+03,
                          3.2227566305839423e+03,
                          3.0346313046911800e+03};

double openQCDWsl2[16] = {9.3412091079710009e+01,
                          3.0854818651585435e+02,
                          3.2619355692323802e+02,
                          3.2558901360119512e+02,
                          3.2932564090843681e+02,
                          3.2531794314030998e+02,
                          3.3105671568595102e+02,
                          3.2270443646197390e+02,
                          3.2061742466075077e+02,
                          3.1723060038309109e+02,
                          3.2007767060738047e+02,
                          3.2460412306725749e+02,
                          3.2968110448055387e+02,
                          3.2012398132057126e+02,
                          3.1669584830412276e+02,
                          3.0766383310510912e+02};

double openQCDYsl0[16] = {0.0000000000000000e+00,
                          4.5515652852019804e+02,
                          5.0374341515676491e+02,
                          4.9834888306821506e+02,
                          5.0732728026591548e+02,
                          4.9801639204602600e+02,
                          5.0183897566338430e+02,
                          4.9250860107779107e+02,
                          4.9764225787923863e+02,
                          4.7930379662469301e+02,
                          4.9378001329798195e+02,
                          5.0301200699232294e+02,
                          5.0542390089340955e+02,
                          4.8935497339842925e+02,
                          4.9395860489544617e+02,
                          4.4862776159253957e+02};

double openQCDYsl2[16] = {0.0000000000000000e+00,
                          1.0208212216268434e+02,
                          1.1091233607192865e+02,
                          1.1261962279721104e+02,
                          1.1691073485407485e+02,
                          1.1626345739620616e+02,
                          1.1674534368987221e+02,
                          1.1252222884621874e+02,
                          1.1157610568334775e+02,
                          1.0678938404798102e+02,
                          1.1167291215356992e+02,
                          1.1472094459891323e+02,
                          1.1501039063720815e+02,
                          1.0925960424790441e+02,
                          1.0649904345718734e+02,
                          1.0079509598393557e+02};

double openQCDWFobsl0[3] = {4.938097e+04, 7.368043e+03, -1.95e-03};

double openQCDWFobsl2[3] = {5.014492e+03, 1.664379e+03, 5.29e-03};

typedef struct _input_sfc
{
  double precision;
  double mass;
  double csw;
  double beta;
} input_sfc;

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
  double gsf, test = 0.;
  struct timeval start, end, etime;
  data_storage_array *store;

  char cnfg_filename[256] = "cnfg/HiRepSF_qcd1_L8T16_b12.0_c1.13295_k0.1298027_r0_id5n1395";

  setup_process(&argc, &argv);
  setup_gauge_fields();

  input_sfc SF_var = {.precision = 1e-16, .beta = 12.0, .mass = -0.14800000308159977, .csw = 1.13295};
  input_WF_meas WF_var = {.tmax = 0.2, .nmeas = 1, .eps = .8e-5, .delta = 1.0e-5, .ittype = RK3_ADAPTIVE};
  set_csw(&SF_var.csw);

  lprintf("MAIN", 0, "Inverter precision = %e\n", SF_var.precision);
  lprintf("MAIN", 0, "Mass = %f\n", SF_var.mass);
  lprintf("MAIN", 0, "beta = %.8f ct = %.8f\n", SF_var.beta, 1.0);

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
  // perform SF measurements.
  gsf = SF_action(SF_var.beta);
  lprintf("SF_action", 10, "gsf = %.10e\n", gsf);

  store = SF_PCAC_wall_corr(SF_var.mass, SF_var.precision, STORE);

  gettimeofday(&end, 0);
  timeval_subtract(&etime, &end, &start);
  lprintf("MAIN", 0, "Correlator measured in [%ld sec %ld usec]\n", etime.tv_sec, etime.tv_usec);

  lprintf("MAIN", 0, "Comparing fA, fP, gA, gP , f1 to the openQCD results found in:\nComparisonLogs/openQCD_qcd1_sf_L8T16_b12.0_c1.13295_k0.1298027_r0_id5.mscf.log\n");

  int idx[3];
  for (int i = 0; i < GLB_T - 2; i++)
  {

    idx[1] = i;
    idx[0] = 1;
    test += fabs(openQCDfA[i] - *data_storage_element(store, 0, idx));

    idx[1] = i;
    idx[0] = 0;
    test += fabs(openQCDfP[i] - *data_storage_element(store, 0, idx));

    idx[1] = GLB_T - 2 - i;
    idx[0] = 1;
    test += fabs(openQCDgA[i] - *data_storage_element(store, 1, idx));

    idx[1] = GLB_T - 2 - i;
    idx[0] = 0;
    test += fabs(openQCDgP[i] - *data_storage_element(store, 1, idx));
  }
  idx[0] = 0;
  test += fabs(openQCDf1 - *data_storage_element(store, 2, idx));

  test /= 4 * (GLB_T - 2) + 1;

  lprintf("TEST", 0, "Cumulative difference: %.2e \n(should be around 1*10^(-10) or so)\n\n", test);
  if (test > 1.e-8 && PID==0)
  {
    lprintf("TEST", 0, "Test failed\n");
    return_value += 1;
  }
  free_data_storage(store);

  gettimeofday(&start, 0);

  store = WF_update_and_measure(WF_var.ittype, u_gauge, &(WF_var.tmax), &(WF_var.eps), &(WF_var.delta), WF_var.nmeas, STORE);
  gettimeofday(&end, 0);
  timeval_subtract(&etime, &end, &start);

  lprintf("MAIN", 0, "WF Observables flowed and measured in [%ld sec %ld usec]\n", etime.tv_sec, etime.tv_usec);

  double comb[2];
  test = 0.0;
  for (int i = 1; i < GLB_T - 1; i++)
  {
    if (i == 1)
    {
      idx[0] = 0;
      idx[1] = i;
      idx[2] = 1;
      comb[0] = *data_storage_element(store, 0, idx);
      idx[2] = 2;
      comb[0] += *data_storage_element(store, 0, idx);

      idx[0] = 1;
      idx[1] = i;
      idx[2] = 1;
      comb[1] = *data_storage_element(store, 0, idx);
      idx[2] = 2;
      comb[1] += *data_storage_element(store, 0, idx);

      comb[0] *= GLB_VOL3 * NG;
      comb[1] *= GLB_VOL3 * NG;
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

      comb[0] *= GLB_VOL3 * NG;
      comb[1] *= GLB_VOL3 * NG;
    }
    test += fabs((comb[0] - openQCDWsl0[i - 1]) / comb[0]);
    test += fabs((comb[1] - openQCDWsl2[i - 1]) / comb[1]);

    if (i > 1)
    {
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

      test += fabs((comb[0] - openQCDYsl0[i - 1]) / comb[0]);
      test += fabs((comb[1] - openQCDYsl2[i - 1]) / comb[1]);
    }
  }

  test /= 4 * (GLB_T - 2) - 1;

  lprintf("MAIN", 0, "Comparing Wl(t=0,T) Wl(t=.2,T) Yl(t=0,T) Yl(t=.2,T) to the openQCD results found in:\nComparisonLogs/openQCD_qcd1_sf_L8T16_b12.0_c1.13295_k0.1298027_r0_id5.ms3.log\n");

  lprintf("TEST", 0, "Cumulative relative difference: %.2e \n(should be around 1*10^(-6) or so)\n\n", test);
  if (test > 1e-6 && PID==0)
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

  test = fabs(1. - 2.0 * GLB_VOL3 * NG * ((GLB_T - 2) * tavg[1] + (GLB_T - 3) * tavg[2]) / openQCDWFobsl0[0]);
  lprintf("TEST", 0, "Wl(t=0) relative difference: %.2e \n(should be around 1*10^(-7) or so)\n\n", test);
  if (test > 1e-6 && PID==0)
  {
    lprintf("TEST", 0, "Test failed\n");
    return_value += 1;
  }

  test = fabs(1.0 - 2 * GLB_VOL3 * NG * (GLB_T - 3) * (tavg[3] + tavg[4]) / openQCDWFobsl0[1]);
  lprintf("TEST", 0, "Yl(t=0) relative difference: %.2e \n(should be around 1*10^(-7) or so)\n\n", test);
  if (test > 1e-6 && PID==0)
  {
    lprintf("TEST", 0, "Test failed\n");
    return_value += 1;
  }

  test = fabs(1.0 - tavg[5] / openQCDWFobsl0[2]);
  lprintf("TEST", 0, "TC(t=0) relative difference: %.2e \n(should be around 1*10^(-3) or so)\n\n", test);
  if (test > 1e-2 && PID==0)
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

  test = fabs(1. - 2.0 * GLB_VOL3 * NG * ((GLB_T - 2) * tavg[1] + (GLB_T - 3) * tavg[2]) / openQCDWFobsl2[0]);
  lprintf("TEST", 0, "Wl(t=0.2) relative difference: %.2e \n(should be around 1*10^(-7) or so)\n\n", test);
  if (test > 1e-5 && PID==0)
  {
    lprintf("TEST", 0, "Test failed\n");
    return_value += 1;
  }

  test = fabs(1.0 - 2 * GLB_VOL3 * NG * (GLB_T - 3) * (tavg[3] + tavg[4]) / openQCDWFobsl2[1]);
  lprintf("TEST", 0, "Yl(t=0.2) relative difference: %.2e \n(should be around 1*10^(-7) or so)\n\n", test);
  if (test > 1e-6 && PID==0)
  {
    lprintf("TEST", 0, "Test failed\n");
    return_value += 1;
  }

  test = fabs(1.0 - tavg[5] / openQCDWFobsl2[2]);
  lprintf("TEST", 0, "TC(t=0.2) relative difference: %.2e \n(should be around 1*10^(-3) or so)\n\n", test);
  if (test > 1e-2 && PID==0)
  {
    lprintf("TEST", 0, "Test failed\n");
    return_value += 1;
  }

  free_data_storage(store);

  WF_free();

  finalize_process();

  return return_value;
}
