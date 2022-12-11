/***************************************************************************\
 * Copyright (c) 2008-2014, Claudio Pica                                    *
 * All rights reserved.                                                     *
 \**************************************************************************/

/**
 * @file geometry.h
 * @brief This file contains information on the geometry of the local 
 * 	  lattice, block decomposed geometry, buffers, etc. All geometric
 * 	  properties that are needed to perform operations correctly
 * 	  have to be contained in the struct geometry_descriptor.
 * 	  In order to encapsulate the complicated geometry structure
 * 	  iteration over lattice sites is simplified through macros
 * 	  given in this file. This simplified performing operations
 * 	  on lattice data.
 */

#ifndef GEOMETRY_H
#define GEOMETRY_H
#ifdef __cplusplus
  extern "C" {
#endif 

#include "geometry_descriptor.h"
#include "geometry_omp.h"
#include "geometry_fuse.h"
#include "geometry_init.h"
#include "cpu_geometry.h"
#include "geometry_check.h"
#include "geometry_maskstate.h"
#ifdef WITH_GPU
   #include "gpu_geometry.h"
#endif

/* this define the width of the borders for parallel dimensions
 * For different actions it must be modified
 * FIXME: Put into global variables, or put more geometry variables from global to
 * here.
 */

#define BORDERSIZE 1


#ifdef __cplusplus
  }
#endif
#endif
