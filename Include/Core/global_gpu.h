// These are device global variables
// must be declared with __device__ or __constant__

#ifndef GLOBAL_GPU_H
#define GLOBAL_GPU_H

#ifdef WITH_GPU

#include "global.h"

#ifdef GPU_MAIN_PROGRAM
#  define GPU_GLB_VAR(type,name) type name
#else
#  define GPU_GLB_VAR(type,name) extern type name
#endif

GPU_GLB_VAR(__device__ __constant__ int,T_EXT_GPU);
GPU_GLB_VAR(__device__ __constant__ int,X_EXT_GPU);
GPU_GLB_VAR(__device__ __constant__ int,Y_EXT_GPU);
GPU_GLB_VAR(__device__ __constant__ int,Z_EXT_GPU);

#define ipt_ext_gpu(t, x, y, z) ipt_gpu[_lexi(T_EXT_GPU, X_EXT_GPU, Y_EXT_GPU, Z_EXT_GPU, t, x, y, z)]
#endif

#endif