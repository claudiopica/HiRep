/***************************************************************************\
 * Copyright (c) 2008, Claudio Pica                                          *   
 * All rights reserved.                                                      * 
 \***************************************************************************/
 
#ifndef SUN_UTILS_H
#define SUN_UTILS_H

#include "update.h"
#include "input_par.h"

/* suN variables */
typedef struct _input_pg {

  double beta;
  int nth, nms, nit, nhb, nor;

  /* for the reading function */
  input_record_t read[7];
  
} input_pg;

#define init_input_pg(varname) \
{ \
  .read={\
    {"beta", "beta = %lf", DOUBLE_T, &(varname).beta},\
    {"nit", "nit = %d", INT_T, &(varname).nit},\
    {"nhb", "nhb = %d", INT_T, &(varname).nhb},\
    {"nor", "nor = %d", INT_T, &(varname).nor},\
    {NULL, NULL, INT_T, NULL}\
  }\
}


/* Flow control variables variables */
typedef struct _pg_flow {
  char run_name[64]; /* name for this run */
  char g_start[64]; /* for gauge fields => unit, random, file */

  int therm;

  char last_conf[64]; /* last conf: can be a number or of the format "+n" */
  char conf_dir[64]; /* directory to store gconfs */
  
  int save_freq; /* save gauge conf if number%save_freq==0 */
  int meas_freq; /* mk measures if number%meas_freq==0 */

  /* these are not actually read from input
   * but inferred from the above
   */
  int start, end;
  input_pg *pg_v;

  /* for the reading function */
  input_record_t read[8];
  
} pg_flow;

#define init_pg_flow(varname) \
{ \
  .read={\
    {"run name", "run name = %s", STRING_T, &((varname).run_name[0])},\
    {"gauge start", "gauge start = %s", STRING_T, &((varname).g_start[0])},\
    {"gauge last conf", "last conf = %s", STRING_T, &((varname).last_conf[0])},\
    {"conf save frequency", "save freq = %d", INT_T, &((varname).save_freq)},\
    {"measure frequency", "meas freq = %d", INT_T, &((varname).meas_freq)},\
    {"config dir", "conf dir = %s", STRING_T, &((varname).conf_dir[0])},\
    {"thermalization", "therm = %d", INT_T, &((varname).therm)},\
    {NULL, NULL, INT_T, NULL}\
  }\
}

int init_mc(pg_flow *rf, char *ifile);
int save_conf(pg_flow *rf, int id);
int end_mc();

#endif /* SUN_UTILS_H */
