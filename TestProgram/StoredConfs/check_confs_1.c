/*******************************************************************************
*
* NOCOMPILE= !BASIC_SF
* NOCOMPILE= !NG=3
* NOCOMPILE= !WITH_EXPCLOVER
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
                        -7.6556680293627077e-01,
                        -6.2388404822446297e-01,
                        -4.2319108709704523e-01,
                        -3.1283235789167463e-01,
                        -2.2847499478249672e-01,
                        -1.6722271330995428e-01,
                        -1.3834298486825891e-01,
                        -1.0660067704385566e-01,
                        -8.2464046163855914e-02,
                        -6.8259655861102858e-02,
                        -5.6403748482926210e-02,
                        -4.9667711537176723e-02,
                        -4.1918419954850746e-02,
                        -3.9428996555828696e-02,
                        -3.6027672398801724e-02};

double openQCDgA[16] = {0.0000000000000000e+00,
                        -3.7973955081079237e-02,
                        -3.9963704219861659e-02,
                        -4.1720758893073515e-02,
                        -4.4276207076463028e-02,
                        -4.6757359008152678e-02,
                        -4.9505714455895702e-02,
                        -5.6220939438667739e-02,
                        -6.0888109286844204e-02,
                        -6.7958127561238646e-02,
                        -7.8956318038499090e-02,
                        -9.0011251543361717e-02,
                        -1.1259591879988308e-01,
                        -1.3948756289204736e-01,
                        -2.0227823087597710e-01,
                        -2.9363475943157114e-01};

double openQCDfP[16] = {-0.0000000000000000e+00,
                        4.3933788150923876e+00,
                        3.7825458812362798e+00,
                        3.0026981593754165e+00,
                        2.2393962436668429e+00,
                        1.6019040817565806e+00,
                        1.0946487326807692e+00,
                        7.7987674019298803e-01,
                        5.8111779626502114e-01,
                        4.2462115913659448e-01,
                        3.0623288608015659e-01,
                        2.1844545662345949e-01,
                        1.5005299814922909e-01,
                        9.9117929568100727e-02,
                        6.5496263231629173e-02,
                        4.3873552687589731e-02};

double openQCDgP[16] = {-0.0000000000000000e+00,
                        4.0468750973430126e-02,
                        4.5025961166660987e-02,
                        5.0162873098094654e-02,
                        5.8687183404137713e-02,
                        7.2497824222916979e-02,
                        9.3695622833957978e-02,
                        1.2404741025170883e-01,
                        1.6249185685331041e-01,
                        2.1849776934094045e-01,
                        2.9688419045501924e-01,
                        4.1930103697620236e-01,
                        6.1054486136239383e-01,
                        9.5718062989575092e-01,
                        1.6073779762057943e+00,
                        2.8549932103844147e+00};

double openQCDf1 = 3.1423987322220417e-02;

double openQCDWsl0[16] = {7.1265572645753355e+02,
                          3.0698932453211769e+03,
                          3.2427634744232314e+03,
                          3.2507359360093219e+03,
                          3.2751138734618721e+03,
                          3.2306453285552625e+03,
                          3.2677284831437205e+03,
                          3.2355890819328647e+03,
                          3.2074931774147403e+03,
                          3.2655558124056656e+03,
                          3.2340493955741285e+03,
                          3.2167146491534227e+03,
                          3.2449495400012497e+03,
                          3.2600416173253202e+03,
                          3.2141719084559754e+03,
                          3.0969574470501811e+03};

double openQCDWsl2[16] = {9.4319583991203288e+01,
                          3.1058547567338132e+02,
                          3.3000203047389283e+02,
                          3.3633189401309716e+02,
                          3.3302055739185244e+02,
                          3.2809663335024464e+02,
                          3.3584590591498943e+02,
                          3.3181597785160051e+02,
                          3.2793135409815528e+02,
                          3.3161310471118350e+02,
                          3.2393143405927827e+02,
                          3.2234730429969630e+02,
                          3.3595629761191327e+02,
                          3.3305918666054788e+02,
                          3.1885061538219423e+02,
                          3.1015218316362416e+02};

double openQCDYsl0[16] = {0.0000000000000000e+00,
                          4.5784744319345668e+02,
                          5.1079139007591067e+02,
                          5.1320649358109449e+02,
                          5.0889096482851846e+02,
                          5.0220267641439256e+02,
                          5.0962828780318904e+02,
                          5.0035172981621406e+02,
                          5.0701588492503748e+02,
                          5.0974549763835472e+02,
                          4.9664782130200120e+02,
                          5.0160302945267517e+02,
                          5.1110414665962162e+02,
                          5.1040941609946248e+02,
                          5.0179010127432372e+02,
                          4.5381536852631200e+02};

double openQCDYsl2[16] = {0.0000000000000000e+00,
                          1.0199465823514454e+02,
                          1.1696027820104261e+02,
                          1.2053185158097313e+02,
                          1.1723167897097471e+02,
                          1.1597044564620620e+02,
                          1.1852207792052910e+02,
                          1.1639467393189501e+02,
                          1.1556100198305930e+02,
                          1.1537261835863023e+02,
                          1.1089367037955670e+02,
                          1.1441508563366952e+02,
                          1.1808764580745937e+02,
                          1.1700000345389182e+02,
                          1.0991134223389042e+02,
                          1.0184094252945448e+02};

double openQCDWFobsl0[3] = {4.972843e+04, 7.495050e+03, 1.53e-01};

double openQCDWFobsl2[3] = {5.098304e+03, 1.710688e+03, 2.72e-02};

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
  int idx[3];

  char cnfg_filename[256] = "cnfg/HiRepSF_qcd1_L8T16_b12.0_c1.13295_k0.1298027_r0_id55n1366";

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

  lprintf("MAIN", 0, "Comparing fA, fP, gA, gP , f1 to the openQCD results found in:\nComparisonLogs/openQCD_qcd1_sf_L8T16_b12.0_c1.13295_k0.1298027_r0_id55.mscf.log\n");

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
  if (test > 1.e-8 && PID == 0)
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

      comb[0] *= VOL3 * NG;
      comb[1] *= VOL3 * NG;
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

      comb[0] *= VOL3 * NG;
      comb[1] *= VOL3 * NG;
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

      comb[0] *= 2 * VOL3 * NG;
      comb[1] *= 2 * VOL3 * NG;

      test += fabs((comb[0] - openQCDYsl0[i - 1]) / comb[0]);
      test += fabs((comb[1] - openQCDYsl2[i - 1]) / comb[1]);
    }
  }

  test /= 4 * (GLB_T - 2) - 1;
  lprintf("MAIN", 0, "Comparing Wl(t=0,T) Wl(t=.2,T) Yl(t=0,T) Yl(t=.2,T) to the openQCD results found in:\nComparisonLogs/openQCD_qcd1_sf_L8T16_b12.0_c1.13295_k0.1298027_r0_id55.ms3.log\n");

  lprintf("TEST", 0, "Cumulative relative difference: %.2e \n(should be around 1*10^(-6) or so)\n\n", test);
  if (test > 1e-6 && PID == 0)
  {
    lprintf("TEST", 0, "Test failed\n");
    return_value += 1;
  }

  double tavg[6];

  lprintf("MAIN", 0, "Comparing Wl(t=0) Wl(t=.2) Yl(t=0) Yl(t=.2) TC(t=0) TC(t=.2) to the openQCD results found in:\nComparisonLogs/openQCD_qcd1_sf_L8T16_b12.0_c1.13295_k0.1298027_r0_id55.ms3.log\n");
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
  if (test > 1e-6 && PID == 0)
  {
    lprintf("TEST", 0, "Test failed\n");
    return_value += 1;
  }

  test = fabs(1.0 - 2 * GLB_VOL3 * NG * (GLB_T - 3) * (tavg[3] + tavg[4]) / openQCDWFobsl0[1]);
  lprintf("TEST", 0, "Yl(t=0) relative difference: %.2e \n(should be around 1*10^(-7) or so)\n\n", test);
  if (test > 1e-6 && PID == 0)
  {
    lprintf("TEST", 0, "Test failed\n");
    return_value += 1;
  }

  test = fabs(1.0 - tavg[5] / openQCDWFobsl0[2]);
  lprintf("TEST", 0, "TC(t=0) relative difference: %.2e \n(should be around 1*10^(-3) or so)\n\n", test);
  if (test > 1e-2 && PID == 0)
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
  if (test > 1e-5 && PID == 0)
  {
    lprintf("TEST", 0, "Test failed\n");
    return_value += 1;
  }

  test = fabs(1.0 - 2 * GLB_VOL3 * NG * (GLB_T - 3) * (tavg[3] + tavg[4]) / openQCDWFobsl2[1]);
  lprintf("TEST", 0, "Yl(t=0.2) relative difference: %.2e \n(should be around 1*10^(-7) or so)\n\n", test);
  if (test > 1e-6)
  {
    lprintf("TEST", 0, "Test failed\n");
    return_value += 1;
  }

  test = fabs(1.0 - tavg[5] / openQCDWFobsl2[2]);
  lprintf("TEST", 0, "TC(t=0.2) relative difference: %.2e \n(should be around 1*10^(-3) or so)\n\n", test);
  if (test > 1e-2)
  {
    lprintf("TEST", 0, "Test failed\n");
    return_value += 1;
  }

  free_data_storage(store);

  WF_free();

  finalize_process();

  return return_value;
}
