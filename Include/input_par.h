/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *   
* All rights reserved.                                                      * 
\***************************************************************************/


#ifndef INPUT_PAR_H
#define INPUT_PAR_H
#include "update.h"

typedef struct _input_par {
  /* global size of lattice and processing grid */
  int GLB_T, GLB_X, GLB_Y, GLB_Z;
  int NP_T, NP_X, NP_Y, NP_Z;

  /* random numbers */
  int rlxd_level, rlxd_seed;
  
  /* simulation parameters */
  rhmc_par rhmc_p;
  int_par int_p;

  /* run parameters */

} input_par;

#endif /* INPUT_PAR_H */
