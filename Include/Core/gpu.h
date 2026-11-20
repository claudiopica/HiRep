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
#define forceinline __forceinline__
#else
#define forceinline __attribute__((always_inline)) inline
#endif

#ifdef WITH_GPU

#include <stdio.h>
#ifndef HIP
#include <cuda.h>
#include <driver_types.h>
#include <cuda_runtime.h>
#include <cuda_runtime_api.h>
#else
#include <hip/hip_runtime.h>
#include <hip/driver_types.h>
#include <hip/hip_runtime_api.h>
#endif

#include "IO/input_par.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    unsigned int gpuID;
    input_record_t read[2]; /* for the reading function */
} input_gpu;

#ifdef HIP
#define cudaMalloc hipMalloc
#define cudaMemcpy hipMemcpy
#define cudaMemset hipMemset
#define cub hipcub
#define cudaGetLastError hipGetLastError
#define cudaErrorInvalidDevice hipErrorInvalidDevice
#define cudaErrorPeerAccessAlreadyEnabled hipErrorPeerAccessAlreadyEnabled
#define cudaErrorInvalidValue hipErrorInvalidValue
#define cudaStreamDestroy hipStreamDestroy
#define cudaMemcpyHostToDevice hipMemcpyHostToDevice
#define cudaMemcpyDeviceToDevice hipMemcpyDeviceToDevice
#define cudaMemcpyHostToHost hipMemcpyHostToHost
#define cudaMemcpyDeviceToHost hipMemcpyDeviceToHost
#define cudaDeviceSynchronize hipDeviceSynchronize
#define cudaError_t hipError_t
#define cudaDeviceProp hipDeviceProp_t
#define cudaSuccess hipSuccess
#define cudaGetErrorString hipGetErrorString
#define cudaStreamCreate hipStreamCreate
#define cudaStream_t hipStream_t
#define cudaGetDeviceProperties hipGetDeviceProperties
#define cudaSetDevice hipSetDevice
#define cudaGetDevice hipGetDevice
#define cudaDeviceCanAccessPeer hipDeviceCanAccessPeer
#define cudaDeviceEnablePeerAccess hipDeviceEnablePeerAccess
#define cudaFree hipFree
#define cudaGetDeviceCount hipGetDeviceCount
#define cudaDriverGetVersion hipDriverGetVersion
#define cudaRuntimeGetVersion hipRuntimeGetVersion
#define cuDeviceGetAttribute hipDeviceGetAttribute
#define CU_DEVICE_ATTRIBUTE_GLOBAL_MEMORY_BUS_WIDTH hipDeviceAttributeMemoryBusWidth
#define CU_DEVICE_ATTRIBUTE_MEMORY_CLOCK_RATE hipDeviceAttributeMemoryClockRate
#define CU_DEVICE_ATTRIBUTE_L2_CACHE_SIZE hipDeviceAttributeL2CacheSize
#define cudaMemcpyToSymbol hipMemcpyToSymbol
#endif

#ifndef M_PI
// Define M_PI if its not part of the C standard
#define M_PI 3.14159265358979323846264338327
#endif

// Allow also to use PI for M_PI
#define PI M_PI

#define init_input_gpu(varname)                                                               \
    {                                                                                         \
        .read = { { "gpuID", "gpuID = %d", INT_T, &(varname) }, { NULL, NULL, INT_T, NULL } } \
    }

#ifdef __cplusplus
}
#endif

#endif
#endif
