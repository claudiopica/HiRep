/***************************************************************************\
 * Copyright (c) 2008, Claudio Pica, Agostino Patella                        *   
 * All rights reserved.                                                      * 
 \***************************************************************************/
 
#ifndef HMC_UTILS_H
#define HMC_UTILS_H

#include "update.h"
#include "input_par.h"

/* HMC variables */
typedef struct _input_hmc {
  /* hmc parameters */
  hmc_par hmc_p;
  /* "Temporary" strings to read list to read variable number of hasenbush parameters  */
  char hasen_steps[2048]; 
  char hasen_dm[2048];
  /* for the reading function */
  input_record_t read[28];
  
} input_hmc;

#define init_input_hmc(varname) \
{ \
  .read={\
    {"beta", "beta = %lf", DOUBLE_T, &(varname).hmc_p.beta},\
    {"nf", "nf = %d", INT_T, &(varname).hmc_p.nf},\
    {"mass", "mass = %lf", DOUBLE_T, &(varname).hmc_p.mass},\
    {"theta_T", "theta_T = %lf", DOUBLE_T, &(varname).hmc_p.theta[0]},\
    {"theta_X", "theta_X = %lf", DOUBLE_T, &(varname).hmc_p.theta[1]},\
    {"theta_Y", "theta_Y = %lf", DOUBLE_T, &(varname).hmc_p.theta[2]},\
    {"theta_Z", "theta_Z = %lf", DOUBLE_T, &(varname).hmc_p.theta[3]},\
    {"hasen_dm", "hasen_dm = %s", STRING_T, &(varname).hasen_dm},\
    {"SF_zf", "SF_zf = %lf", DOUBLE_T, &(varname).hmc_p.SF_zf},\
    {"SF_ds", "SF_ds = %lf", DOUBLE_T, &(varname).hmc_p.SF_ds},\
    {"SF_sign", "SF_sign = %d", INT_T, &(varname).hmc_p.SF_sign},\
    {"SF_ct", "SF_ct = %lf", DOUBLE_T, &(varname).hmc_p.SF_ct}, \
    {"n_MT_prec", "n_MT_prec = %lf", DOUBLE_T, &(varname).hmc_p.n_MT_prec},\
    {"h_MT_prec", "h_MT_prec = %lf", DOUBLE_T, &(varname).hmc_p.h_MT_prec},\
    {"n_MT_prec_flt", "n_MT_prec_flt = %lf", DOUBLE_T, &(varname).hmc_p.n_MT_prec_flt},\
    {"h_MT_prec_flt", "h_MT_prec_flt = %lf", DOUBLE_T, &(varname).hmc_p.h_MT_prec_flt},\
    {"MD_prec", "MD_prec = %lf", DOUBLE_T, &(varname).hmc_p.MD_prec},\
    {"HB_prec", "HB_prec = %lf", DOUBLE_T, &(varname).hmc_p.HB_prec},\
    {"n_force_prec", "n_force_prec = %lf", DOUBLE_T, &(varname).hmc_p.n_force_prec},\
    {"h_force_prec", "h_force_prec = %lf", DOUBLE_T, &(varname).hmc_p.h_force_prec},\
    {"n_force_prec_flt", "n_force_prec_flt = %lf", DOUBLE_T, &(varname).hmc_p.n_force_prec_flt},\
    {"h_force_prec_flt", "h_force_prec_flt = %lf", DOUBLE_T, &(varname).hmc_p.h_force_prec_flt},\
    {"tlen", "tlen = %lf", DOUBLE_T, &(varname).hmc_p.tlen},\
    {"nsteps", "nsteps = %u", UNSIGNED_T, &(varname).hmc_p.nsteps},\
    {"hsteps", "hsteps = %s", STRING_T, &(varname).hasen_steps},	\
    {"gsteps", "gsteps = %u", UNSIGNED_T, &(varname).hmc_p.gsteps},\
    {"hasenbusch", "hasenbusch = %d", INT_T, &(varname).hmc_p.hasenbusch},\
    {NULL, NULL, 0, NULL}\
  }\
}

/* Flow control variables variables */
typedef struct _hmc_flow {
  char run_name[64]; /* name for this run */
  char g_start[64]; /* for gauge fields => unit, random, file */

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
  input_record_t read[7];
  
} hmc_flow;

#define init_hmc_flow(varname) \
{ \
  .read={\
    {"run name", "run name = %s", STRING_T, &((varname).run_name[0])},\
    {"gauge start", "gauge start = %s", STRING_T, &((varname).g_start[0])},\
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
