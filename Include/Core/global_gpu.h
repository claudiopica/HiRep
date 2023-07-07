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

#ifndef GLOBAL_GPU_H
#define GLOBAL_GPU_H

#ifdef WITH_GPU

#include "core_utils.h"
#include "gpu.h"
#include "geometry.h"

GLB_VAR(__device__ __constant__ int, T_EXT_GPU);
GLB_VAR(__device__ __constant__ int, X_EXT_GPU);
GLB_VAR(__device__ __constant__ int, Y_EXT_GPU);
GLB_VAR(__device__ __constant__ int, Z_EXT_GPU);
GLB_VAR(__device__ __constant__ int, T_GPU);
GLB_VAR(__device__ __constant__ int, X_GPU);
GLB_VAR(__device__ __constant__ int, Y_GPU);
GLB_VAR(__device__ __constant__ int, Z_GPU);

GLB_VAR(__device__ __constant__ char, UP_MASK, = T_UP_MASK | X_UP_MASK | Y_UP_MASK | Z_UP_MASK);
GLB_VAR(__device__ __constant__ char, DN_MASK, = T_DN_MASK | X_DN_MASK | Y_DN_MASK | Z_DN_MASK);
GLB_VAR(__device__ __constant__ char, T_MASK, = T_UP_MASK | T_DN_MASK);
GLB_VAR(__device__ __constant__ char, X_MASK, = X_UP_MASK | X_DN_MASK);
GLB_VAR(__device__ __constant__ char, Y_MASK, = Y_UP_MASK | Y_DN_MASK);
GLB_VAR(__device__ __constant__ char, Z_MASK, = Z_UP_MASK | Z_DN_MASK);

#if defined(BC_T_THETA) || defined(BC_X_THETA) || defined(BC_Y_THETA) || defined(BC_Z_THETA)
#include "hr_complex.h"
GLB_VAR(__device__ __constant__ hr_complex, eitheta_gpu[4]);
#endif

#endif

#endif