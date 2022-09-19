/***************************************************************************\
 * Copyright (c) 2019, Antonio Rago                                          *   
 * All rights reserved.                                                      * 
 \***************************************************************************/

#ifndef SUN_UTILS_H
#define SUN_UTILS_H
#include "wilsonflow.h"
#include "input_par.h"
#define GENERIC_MAX(x, y) ((x) > (y) ? (x) : (y))

#define ENSURE_int(i) _Generic((i), int \
                               : (i))
#define ENSURE_float(f) _Generic((f), float \
                                 : (f))

#define MAX(type, x, y) \
  (type) GENERIC_MAX(ENSURE_##type(x), ENSURE_##type(y))

/* suN variables */
typedef struct _input_pg
{

  double beta, anisotropy;
  int nhb, nor;
  /* for the reading function */
  input_record_t read[5];

} input_pg;

#define init_input_pg(varname)                                             \
  {                                                                        \
    .read = {                                                              \
      {"beta", "beta = %lf", DOUBLE_T, &(varname).beta},                   \
      {"anisotropy", "anisotropy = %lf", DOUBLE_T, &(varname).anisotropy}, \
      {"nhb", "nhb = %d", INT_T, &(varname).nhb},                          \
      {"nor", "nor = %d", INT_T, &(varname).nor},                          \
      {NULL, NULL, INT_T, NULL}                                            \
    }                                                                      \
  }

/* Polyakov variables */

typedef struct _input_poly
{
  char make[256];

  /* for the reading function */
  input_record_t read[2];
} input_poly;

#define init_input_poly(varname)                                         \
  {                                                                      \
    .read = {                                                            \
      {"make Polyakov", "polyakov:make = %s", STRING_T, (varname).make}, \
      {NULL, NULL, 0, NULL}                                              \
    }                                                                    \
  }

/* WF variables */

typedef struct _input_WF
{
  char make[256];
  double tmax;
  int nmeas;
  WF_integrator_type ittype;
  double eps;
  double delta;
  double anisotropy;

  /* for the reading function */
  input_record_t read[8];
} input_WF;

#define init_input_WF(varname)                                                     \
  {                                                                                \
    .read = {                                                                      \
      {"make Wilson Flow", "WF:make = %s", STRING_T, (varname).make},              \
      {"WF integrator type", "WF:integrator = %d", INT_T, &((varname).ittype)},    \
      {"WF max integration time", "WF:tmax = %lf", DOUBLE_T, &((varname).tmax)},   \
      {"WF number of measures", "WF:nmeas = %d", INT_T, &((varname).nmeas)},       \
      {"WF initial epsilon", "WF:eps = %lf", DOUBLE_T, &((varname).eps)},          \
      {"WF delta", "WF:delta = %lf", DOUBLE_T, &((varname).delta)},                \
      {"WF anisotropy", "WF:anisotropy = %lf", DOUBLE_T, &((varname).anisotropy)}, \
      {NULL, NULL, 0, NULL}                                                        \
    }                                                                              \
  }

/* Flow control variables variables */
typedef struct _pg_flow
{
  char run_name[64]; /* name for this run */
  char g_start[64];  /* for gauge fields => unit, random, file */

  int therm;

  char last_conf[64]; /* last conf: can be a number or of the format "+n" */
  char conf_dir[64];  /* directory to store gconfs */

  int save_freq; /* save gauge conf if number%save_freq==0 */
  int nit;
  /* these are not actually read from input
   * but inferred from the above
   */
  int start, end;
  input_pg *pg_v;

  input_WF *wf;

  input_poly *poly;
  /* for the reading function */
  input_record_t read[8];

} pg_flow;

#define init_pg_flow(varname)                                                        \
  {                                                                                  \
    .read = {                                                                        \
      {"thermalization", "therm = %d", INT_T, &((varname).therm)},                   \
      {"run name", "run name = %s", STRING_T, &((varname).run_name[0])},             \
      {"number of updates between each measure", "nit = %d", INT_T, &(varname).nit}, \
      {"conf save frequency", "save freq = %d", INT_T, &((varname).save_freq)},      \
      {"gauge start", "gauge start = %s", STRING_T, &((varname).g_start[0])},        \
      {"gauge last conf", "last conf = %s", STRING_T, &((varname).last_conf[0])},    \
      {"config dir", "conf dir = %s", STRING_T, &((varname).conf_dir[0])},           \
      {NULL, NULL, INT_T, NULL}                                                      \
    }                                                                                \
  }

int init_mc(pg_flow *rf, char *ifile);
int save_conf(pg_flow *rf, int id);

#endif /* SUN_UTILS_H */
