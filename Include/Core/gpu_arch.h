/***************************************************************************\
* Copyright (c) 2024, Sofie Martins                                         *
* All rights reserved.                                                      *
\***************************************************************************/

#ifndef GPU_ARCH_H
#define GPU_ARCH_H

#ifdef WITH_GPU

#include "gpu.h"

#ifdef HIP
#define cudaMalloc hipMalloc
#define cudaMemcpy hipMemcpy
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
#endif

#endif