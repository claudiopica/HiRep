/***************************************************************************\
* Copyright (c) 2022, Claudio Pica                                          *   
* All rights reserved.                                                      * 
\***************************************************************************/
/**
 * @file error_gpu.h
 * @brief Error checking for GPU
 */

#ifndef ERROR_GPU_H
#define ERROR_GPU_H

#ifdef WITH_GPU

#include "gpu.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Check last error
 *
 * @param file          File where the last error was raised
 * @param line          Line where the last error was raised
 */
void __cudaCheckError(const char *func, const char *file, int line);

/**
 * @brief Check CUDA call and log error message on failure 
 *        This may be more performant than CHECK_CUDA 
 *           -> TODO: replace (the other function) in the future (SAM)
 *
 * @param err           Function call that should be checked
 */
#define CudaSafeCall(err) __cudaSafeCall(err, __func__, __FILE__, __LINE__)

/**
 * @brief Check last error after CUDA calls
 */
#define CudaCheckError() __cudaCheckError((const char *)__func__, (const char *)__FILE__, (int)__LINE__)

/**
 * @brief Check CUDA call and log error message on failure
 *
 * @param cudaError_t   Error return type from CUDA call
 * @param file          File where the exception was raised
 * @param line          Line where the exception was raised
 */
void __cudaSafeCall(cudaError_t err, const char *func, const char *file, const int line);

/**
 * @brief Check CUDA call and log error message on failure.
 *
 * @param call           Function call that should be checked.
 */
#define CHECK_CUDA(call) CudaSafeCall(call)

#ifdef __cplusplus
}
#endif
#endif
#endif
