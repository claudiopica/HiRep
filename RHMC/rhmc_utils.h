/***************************************************************************\
 * Copyright (c) 2008, Claudio Pica, Agostino Patella                        *   
 * All rights reserved.                                                      * 
 \***************************************************************************/
 
#ifndef RHMC_UTILS_H
#define RHMC_UTILS_H

#include "update.h"
#include "input_par.h"

/* RHMC variables */
typedef struct _input_rhmc {
  /* rhmc parameters */
  rhmc_par rhmc_p;

  /* for the reading function */
  input_record_t read[20];
  
} input_rhmc;

#define init_input_rhmc(varname) \
{ \
  .read={\
    {"beta", "beta = %lf", DOUBLE_T, &(varname).rhmc_p.beta},\
    {"nf", "nf = %d", INT_T, &(varname).rhmc_p.nf},\
    {"mass", "mass = %lf", DOUBLE_T, &(varname).rhmc_p.mass},\
    {"theta_T", "theta_T = %lf", DOUBLE_T, &(varname).rhmc_p.theta[0]},\
    {"theta_X", "theta_X = %lf", DOUBLE_T, &(varname).rhmc_p.theta[1]},\
    {"theta_Y", "theta_Y = %lf", DOUBLE_T, &(varname).rhmc_p.theta[2]},\
    {"theta_Z", "theta_Z = %lf", DOUBLE_T, &(varname).rhmc_p.theta[3]},\
    {"SF_zf", "SF_zf = %lf", DOUBLE_T, &(varname).rhmc_p.SF_zf},\
    {"SF_ds", "SF_ds = %lf", DOUBLE_T, &(varname).rhmc_p.SF_ds},\
    {"SF_sign", "SF_sign = %d", INT_T, &(varname).rhmc_p.SF_sign},\
    {"SF_ct", "SF_ct = %lf", DOUBLE_T, &(varname).rhmc_p.SF_ct}, \
    {"MT_prec", "MT_prec = %lf", DOUBLE_T, &(varname).rhmc_p.MT_prec},\
    {"MD_prec", "MD_prec = %lf", DOUBLE_T, &(varname).rhmc_p.MD_prec},\
    {"HB_prec", "HB_prec = %lf", DOUBLE_T, &(varname).rhmc_p.HB_prec},\
    {"force_prec", "force_prec = %lf", DOUBLE_T, &(varname).rhmc_p.force_prec},\
    {"n_pf", "n_pf = %u", UNSIGNED_T, &(varname).rhmc_p.n_pf},\
    {"tlen", "tlen = %lf", DOUBLE_T, &(varname).rhmc_p.tlen},\
    {"nsteps", "nsteps = %u", UNSIGNED_T, &(varname).rhmc_p.nsteps},\
    {"gsteps", "gsteps = %u", UNSIGNED_T, &(varname).rhmc_p.gsteps},\
    {NULL, NULL, 0, NULL}\
  }\
}

/* Flow control variables variables */
typedef struct _rhmc_flow {
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
  input_rhmc *rhmc_v;

  /* for the reading function */
  input_record_t read[8];
  
} rhmc_flow;

#define init_rhmc_flow(varname) \
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

int init_mc(rhmc_flow *rf, char *ifile);
int save_conf(rhmc_flow *rf, int id);
int end_mc();

#endif /* RHMC_UTILS_H */
