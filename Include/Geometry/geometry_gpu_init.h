/**
 * @file geometry_gpu_init.h
 * @brief Initialization functions, that determine all important parameters of the 
 *        geometry, so that communications and operations can be completed correctly.
 */
/// Headerfile for:
/// - init_gpu.c
/// - geometry_init_gpu.c

#ifndef GEOMETRY_INIT_GPU_H
#define GEOMETRY_INIT_GPU_H

#ifdef WITH_GPU

#include "gpu.h"

#ifdef __cplusplus
   extern "C" {
#endif

//init_gpu.c
/**
 * @brief Call this in an init function to setup available graphics card for
 *        use. This also logs information on available software and hardware.
 * 
 * @param input_gpu             A struct containing information on the current active
 *                              GPU
 */
void init_gpu(input_gpu gpu_var);

// geometry_init_gpu.c
void init_neighbors_gpu(void);

#ifdef __cplusplus
   }
#endif
#endif
#endif


