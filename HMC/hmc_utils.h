/***************************************************************************\
 * Copyright (c) 2008, Claudio Pica, Agostino Patella                        *   
 * All rights reserved.                                                      * 
 \***************************************************************************/
 
#ifndef HMC_UTILS_H
#define HMC_UTILS_H

#include "update.h"
#include "input_par.h"
#include "gpu.h"

/* HMC variables */
typedef struct _input_hmc {
  /* hmc parameters */
  rhmc_par rhmc_p;

  /* for the reading function */
  input_record_t read[17];
  
} input_hmc;

#define init_input_hmc(varname) \
{ \
  .read={\
    {"beta", "beta = %lf", DOUBLE_T, &(varname).rhmc_p.beta},\
    {"nf", "nf = %d", INT_T, &(varname).rhmc_p.nf},\
    {"n_pf", "n_pf = %d", INT_T, &(varname).rhmc_p.n_pf},\
    {"mass", "mass = %lf", DOUBLE_T, &(varname).rhmc_p.mass},\
    {"SF_zf", "SF_zf = %lf", DOUBLE_T, &(varname).rhmc_p.SF_zf},\
    {"SF_ds", "SF_ds = %lf", DOUBLE_T, &(varname).rhmc_p.SF_ds},\
    {"SF_sign", "SF_sign = %d", INT_T, &(varname).rhmc_p.SF_sign},\
    {"SF_ct", "SF_ct = %lf", DOUBLE_T, &(varname).rhmc_p.SF_ct}, \
    {"MT_prec", "MT_prec = %lf", DOUBLE_T, &(varname).rhmc_p.MT_prec},\
    {"MD_prec", "MD_prec = %lf", DOUBLE_T, &(varname).rhmc_p.MD_prec},\
    {"HB_prec", "HB_prec = %lf", DOUBLE_T, &(varname).rhmc_p.HB_prec},\
    {"force_prec", "force_prec = %lf", DOUBLE_T, &(varname).rhmc_p.force_prec},\
    {"force_prec_flt", "force_prec_flt = %lf", DOUBLE_T, &(varname).rhmc_p.force_prec_flt},\
    {"tlen", "tlen = %lf", DOUBLE_T, &(varname).rhmc_p.tlen},\
    {"nsteps", "nsteps = %u", UNSIGNED_T, &(varname).rhmc_p.nsteps},\
    {"gsteps", "gsteps = %u", UNSIGNED_T, &(varname).rhmc_p.gsteps},\
    {NULL, NULL, 0, NULL}\
  }\
}

/* Flow control variables variables */
typedef struct _hmc_flow {
  char run_name[64]; /* name for this run */
  char g_start[64]; /* for gauge fields => unit, random, file */
  char r_start[64]; /* for ranlux: name of the state file  */ 

  char last_conf[64]; /* last conf: can be a number or of the format "+n" */
  char conf_dir[64]; /* directory to store gconfs */
  
  int save_freq; /* save gauge conf if number%save_freq==0 */
  int meas_freq; /* mk measures if number%meas_freq==0 */

  /* these are not actually read from input
   * but inferred from the above
   */
  int start, end;
  input_hmc *hmc_v;

  /* for the reading function */
  input_record_t read[8];
  
} hmc_flow;

#define init_hmc_flow(varname) \
{ \
  .read={\
    {"run name", "run name = %s", STRING_T, &((varname).run_name[0])},\
    {"gauge start", "gauge start = %s", STRING_T, &((varname).g_start[0])},\
    {"ranlxd start", "ranlxd start = %s", STRING_T, &((varname).r_start[0])},\
    {"gauge last conf", "last conf = %s", STRING_T, &((varname).last_conf[0])},\
    {"conf save frequency", "save freq = %d", INT_T, &((varname).save_freq)},\
    {"measure frequency", "meas freq = %d", INT_T, &((varname).meas_freq)},\
    {"config dir", "conf dir = %s", STRING_T, &((varname).conf_dir[0])},\
    {NULL, NULL, 0, NULL}\
  }\
}

int init_mc(hmc_flow *rf, char *ifile);
int save_conf(hmc_flow *rf, int id);
int end_mc();

#endif /* HMC_UTILS_H */
