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
#endif

#endif