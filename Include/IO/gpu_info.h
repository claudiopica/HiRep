/***************************************************************************\
* Copyright (c) 2012, 2022, Ari Hietanen, Sofie Martins                     *   
* All rights reserved.                                                      * 
\***************************************************************************/

/// Headerfile for:
/// - gpu_info.c

/**
 * @file gpu_info.h
 * @brief Functions that print collective information on cluster 
 *        drivers, necessary software and available hardware.
 *        To be logged to simulation outfiles.
 */

#ifndef GPU_INFO_H
#define GPU_INFO_H

#ifdef WITH_GPU

#include "gpu.h"

#ifdef __cplusplus
    extern "C" {
#endif


typedef struct cudaDeviceProp cudaDeviceProp;

/**
 * @brief Query number of GPUs and print related information
 *
 * @param input_gpu             A struct containing parameters on the current GPU.
 */
void print_device_count_info(input_gpu);
/**
 * @brief Print CUDA driver version information
 */
void print_driver_info(cudaDeviceProp);
/**
 * @brief Print CUDA runtime version information
 *
 * @param input_gpu             A struct containing parameters on the current GPU.
 */
void print_runtime_info(cudaDeviceProp);
/**
 * @brief Print Global memory information including bandwidth 
 *        paramters and supported features
 *
 * @param cudaDeviceProp        A CUDA class containing information on the device.
 * @param input_gpu             A struct containing parameters on the current GPU.
 */
void print_global_memory_info(cudaDeviceProp, input_gpu);
/**
 * @brief Print information on shared memory
 *
 * @param cudaDeviceProp        A CUDA class containing information on the device.
 */
void print_shared_memory_info(cudaDeviceProp);
/**
 * @brief Print (L2) cache information
 *
 * @param cudaDeviceProp        A CUDA class containing information on the device.
 * @param input_gpu             A struct containing parameters on the current GPU.
 */
void print_cache_info(cudaDeviceProp, input_gpu);
/**
 * @brief Print information on constant memory, in particular
 *        amount available, alignment, texture and layered texture 
 *        memory paramters
 *
 * @param cudaDeviceProp        A CUDA class containing information on the device.
 */
void print_constant_memory_info(cudaDeviceProp);
/**
 * @brief Prints all memory related info, here on global, shared
 *        memory, cache and constant memory, see resp. functions
 *        for more information.
 *
 * @param cudaDeviceProp        A CUDA class containing information on the device.
 * @param input_gpu             A struct containing parameters on the current GPU.
 */
void print_memory_info(cudaDeviceProp, input_gpu);
/**
 * @brief Prints information on compute parameters
 *
 * @param cudaDeviceProp        A CUDA class containing information on the device.
 * @param input_gpu             A struct containing parameters on the current GPU.
 */
void print_compute_info(cudaDeviceProp, input_gpu);
/**
 * @brief Checks for a number of other supported features and
 *        prints information on them.
 *
 * @param cudaDeviceProp        A CUDA class containing information on the device.
 */
void print_supported_features(cudaDeviceProp);
/**
 * @brief Prints all information on hardware, meaning memory, 
 *        compute and features (no driver and runtime info).
 *        This assumes CUDART >= 4000.
 *
 * @param input_gpu             A struct containing parameters on the current GPU.
 */
void print_hardware_info(cudaDeviceProp, input_gpu);
    
#ifdef __cplusplus
    }
#endif
#endif
#endif
