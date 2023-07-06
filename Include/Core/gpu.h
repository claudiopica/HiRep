/**
 * @file gpu.h
 * @brief Basic gpu imports and structs. Include this in files that define GPU logic.
 */
#ifndef GPU_H
#define GPU_H

#ifdef WITH_GPU
#define visible __host__ __device__
#define deviceonly __device__
#else
//#define visible inline
#define visible
#define deviceonly
#endif

#ifdef WITH_GPU

#include <stdio.h>
#include <cuda.h>
#include <driver_types.h>
#include <cuda_runtime.h>
#include <cuda_runtime_api.h>

#include "IO/input_par.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    unsigned int gpuID;
    input_record_t read[2]; /* for the reading function */
} input_gpu;

#define init_input_gpu(varname)                                                               \
    {                                                                                         \
        .read = { { "gpuID", "gpuID = %d", INT_T, &(varname) }, { NULL, NULL, INT_T, NULL } } \
    }

#ifdef __cplusplus
}
#endif

#endif
#endif
