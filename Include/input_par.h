/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *   
* All rights reserved.                                                      * 
\***************************************************************************/


#ifndef INPUT_PAR_H
#define INPUT_PAR_H

typedef enum _datatype_t {
  INT_T,
  UNSIGNED_T,
  DOUBLE_T,
  STRING_T
} datatype_t;

typedef struct _input_record_t {
  char *name;
  char *descr;
  datatype_t type;
  void *ptr;
} input_record_t;

/* Global or common variables */
typedef struct _input_glb {
  /* global size of lattice and processing grid */
  /* THIS ARE DEFINED GLOBALLY !!! */
  /* int GLB_T, GLB_X, GLB_Y, GLB_Z; */
  /* int NP_T, NP_X, NP_Y, NP_Z; */
  /* int N_REP; */

  /* for the reading function */
  input_record_t read[10];
  
} input_glb;

#define init_input_glb(varname) \
{ \
  .read={\
    {"GLB_T", "GLB_T = %d", INT_T, &GLB_T},\
    {"GLB_X", "GLB_X = %d", INT_T, &GLB_X},\
    {"GLB_Y", "GLB_Y = %d", INT_T, &GLB_Y},\
    {"GLB_Z", "GLB_Z = %d", INT_T, &GLB_Z},\
    {"NP_T", "NP_T = %d", INT_T, &NP_T},\
    {"NP_X", "NP_X = %d", INT_T, &NP_X},\
    {"NP_Y", "NP_Y = %d", INT_T, &NP_Y},\
    {"NP_Z", "NP_Z = %d", INT_T, &NP_Z},\
    {"N_REP", "N_REP = %d", INT_T, &N_REP},\
    {NULL, NULL, INT_T, NULL}				    \
  }\
}

/* Global or common variables */
typedef struct _input_rlx {
    
    /* random numbers */
    int rlxd_level, rlxd_seed;
    char rlxd_state[256];
    char rlxd_start[256];
    
    /* for the reading function */
    input_record_t read[5];
    
} input_rlx;

#define init_input_rlx(varname) \
{ \
.read={\
{"ranlux level", "rlx_level = %d", INT_T, &(varname).rlxd_level},\
{"ranlux seed", "rlx_seed = %d", INT_T, &(varname).rlxd_seed},\
{"ranlux state", "rlx_state = %s", STRING_T, &(varname).rlxd_state},\
{"ranlux start", "rlx_start = %s", STRING_T, &(varname).rlxd_start},\
{NULL, NULL, INT_T, NULL}				    \
},\
.rlxd_state=""\
}

/* Logger global variables */
typedef struct _input_logger {
/* Logger level defined globally */
/* They are defined at imput level */
/* If you need to separate the log level for a channel insert it here */

  int def_log_lvl;
  int inverter_log_lvl;
  int forcestat_log_lvl;
  /* for the reading function */
  input_record_t read[4];
  
} input_logger;

#define init_input_logger(varname) \
{ \
  .read={\
    {"Default logger level", "log:default = %d", INT_T, &(varname).def_log_lvl},\
    {"Inverter logger level", "log:inverter = %d", INT_T, &(varname).inverter_log_lvl},\
    {"Forcestat logger level", "log:forcestat = %d", INT_T, &(varname).forcestat_log_lvl},\
    {NULL, NULL, INT_T, NULL}\
  }, \
  .def_log_lvl=-1,\
  .inverter_log_lvl=-1,\
  .forcestat_log_lvl=-1\
}


#endif /* INPUT_PAR_H */
