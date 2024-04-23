/**
 * @file error_gpu.c
 * @brief Error checking on the GPU
 */

#ifdef WITH_GPU
#include "gpu.h"
#include "error.h"
#include "io.h"

void __cudaSafeCall(cudaError_t err, const char *func, const char *file, const int line) {
#ifndef CUDA_NO_CHECK_ERROR
    do {
        if (cudaSuccess != err) {
            lprintf("CUDA", 0, "cudaSafeCall() failed in %s, at %s:%i\n", func, file, line);
            error((cudaSuccess != err), 1, "CudaSafeCall", cudaGetErrorString(err));
        }
    } while (0);
#endif
    return;
}

void __cudaCheckError(const char *func, const char *file, int line) /*TODO: inline void? (SAM) */
{
#ifndef CUDA_NO_CHECK_ERROR
    do {
        cudaError_t err = cudaGetLastError();
        if (cudaSuccess != err) {
            lprintf("CUDA", 0, "cudaCheckError() failed in %s at %s:%i\n", func, file, line);
            error((cudaSuccess != err), 1, "CudaCheckError", cudaGetErrorString(err));
        }

        // More careful checking. However, this will affect performance.
        // Comment if not needed.
        err = cudaDeviceSynchronize();
        if (cudaSuccess != err) {
            lprintf("CUDA", 0, "cudaCheckError() with sync failed in %s at %s:%i\n", func, file, line);
            error((cudaSuccess != err), 1, "CudaCheckError with sync", cudaGetErrorString(err));
        }
    } while (0);
#else
    cudaDeviceSynchronize();
#endif
    return;
}

#endif
