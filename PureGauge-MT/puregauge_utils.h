/***************************************************************************\
 * Copyright (c) 2008, Claudio Pica, Agostino Patella                        *   
 * All rights reserved.                                                      * 
 \***************************************************************************/
 
#ifndef RHMC_UTILS_H
#define RHMC_UTILS_H

#include "update.h"
#include "input_par.h"

/* GAUGE variables */
typedef struct _input_puregauge {
  /* puregauge parameters */
  puregauge_par puregauge_p;

  /* for the reading function */
  input_record_t read[6];
  
} input_puregauge;

#define init_input_puregauge(varname) \
{ \
  .read={\
    {"beta", "beta = %lf", DOUBLE_T, &(varname).puregauge_p.beta},\
    {"MT_prec", "MT_prec = %lf", DOUBLE_T, &(varname).puregauge_p.MT_prec},\
    {"MD_prec", "MD_prec = %lf", DOUBLE_T, &(varname).puregauge_p.MD_prec},\
    {"tlen", "tlen = %lf", DOUBLE_T, &(varname).puregauge_p.tlen},\
    {"gsteps", "gsteps = %u", UNSIGNED_T, &(varname).puregauge_p.gsteps},\
    {NULL, NULL, INT_T, NULL}\
  }\
}

/* Flow control variables variables */
typedef struct _puregauge_flow {
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
  input_puregauge *puregauge_v;

  /* for the reading function */
  input_record_t read[8];
  
} puregauge_flow;

#define init_puregauge_flow(varname) \
{ \
  .read={\
    {"run name", "run name = %s", STRING_T, &((varname).run_name[0])},\
    {"gauge start", "gauge start = %s", STRING_T, &((varname).g_start[0])},\
    {"ranlxd start", "ranlxd start = %s", STRING_T, &((varname).r_start[0])},\
    {"gauge last conf", "last conf = %s", STRING_T, &((varname).last_conf[0])},\
    {"conf save frequency", "save freq = %d", INT_T, &((varname).save_freq)},\
    {"measure frequency", "meas freq = %d", INT_T, &((varname).meas_freq)},\
    {"config dir", "conf dir = %s", STRING_T, &((varname).conf_dir[0])},\
    {NULL, NULL, INT_T, NULL}\
  }\
}

int init_mc(puregauge_flow *rf, char *ifile);
int save_conf(puregauge_flow *rf, int id);
int end_mc();

#endif /* RHMC_UTILS_H */
