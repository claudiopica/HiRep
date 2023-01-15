/***************************************************************************\
* Copyright (c) 2022, Claudio Pica, Sofie Martins                           *
* All rights reserved.                                                      *
\***************************************************************************/

/**
 * @file global_gpu.h
 * @brief Global variables need to be declared once in the main program
 *        and then again as extern variables in the files that use them.
 *        This can be achieved by just including this header 
 *        with any other modification. A macro will automatically
 *        declare the variable in the main program without the 
 *        extern modifier and with it everywhere else.
 */

// These are device global variables
// must be declared with __device__ or __constant__
#ifndef GLOBAL_GPU_H
#define GLOBAL_GPU_H

#ifdef WITH_GPU

#include "core_utils.h"
#include "gpu.h"

GLB_VAR(__device__ __constant__ int,T_EXT_GPU);
GLB_VAR(__device__ __constant__ int,X_EXT_GPU);
GLB_VAR(__device__ __constant__ int,Y_EXT_GPU);
GLB_VAR(__device__ __constant__ int,Z_EXT_GPU);

#if defined(BC_T_THETA) || defined(BC_X_THETA) || defined(BC_Y_THETA) || defined(BC_Z_THETA)
#include "hr_complex.h"
GLB_VAR(__device__ __constant__ hr_complex,eitheta_gpu[4]);
#endif

#endif

#endif