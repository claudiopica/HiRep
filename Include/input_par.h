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

  /* random numbers */
  int rlxd_level, rlxd_seed;

  /* for the reading function */
  input_record_t read[11];
  
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
    {"ranlux level", "level = %d", INT_T, &(varname).rlxd_level},\
    {"ranlux seed", "seed = %d", INT_T, &(varname).rlxd_seed},\
    {NULL, NULL, 0, NULL}\
  }\
}


/* Mesons parameters */
typedef struct _input_mesons {
  char mstring[256];

  /* for the reading function */
  input_record_t read[2];
  
} input_mesons;

#define init_input_mesons(varname) \
{ \
  .read={\
    {"quark quenched masses", "mes:masses = %s", STRING_T, (varname).mstring},\
    {NULL, NULL, 0, NULL}\
  }\
}


#endif /* INPUT_PAR_H */
