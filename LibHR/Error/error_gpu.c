/**
 * @file error_gpu.c
 * @brief Error checking on the GPU
 */

#ifdef WITH_GPU
#include "gpu.h"
#include "error.h"

void __cudaSafeCall( cudaError_t err, const char *file, const int line )
{
#ifdef CUDA_CHECK_ERROR

#pragma warning( push )
#pragma warning( disable: 4127 ) // Prevent warning on do-while(0);

    do
    {
        if ( cudaSuccess != err )
        {
            lprintf("CUDA",0,"cudaSafeCall() failed at %s:%i\n", file, line);
            error((cudaSuccess != err),1,"CudaSafeCall", cudaGetErrorString( err ));
        }
    } while ( 0 );

#pragma warning( pop )

#endif  // CUDA_CHECK_ERROR

    return;
}

void __cudaCheckError( const char *file, int line ) /*TODO: inline void? (SAM) */ 
{
#ifdef CUDA_CHECK_ERROR

#pragma warning( push )
#pragma warning( disable: 4127 ) // Prevent warning on do-while(0);

    do
    {
        cudaError_t err = cudaGetLastError();
        if ( cudaSuccess != err )
        {
            lprintf("CUDA",0,"cudaCheckError() failed at %s:%i\n", file, line);
            error((cudaSuccess != err),1,"CudaCheckError", cudaGetErrorString( err ));
        }

        // More careful checking. However, this will affect performance.
        // Comment if not needed.
        err = cudaThreadSynchronize();
        if( cudaSuccess != err )
        {
            lprintf("CUDA",0,"cudaCheckError() with sync failed at %s:%i\n", file, line);
            error((cudaSuccess != err),1,"CudaCheckError with sync", cudaGetErrorString( err ));
        }
    } while ( 0 );

#pragma warning( pop )

#endif // CUDA_CHECK_ERROR

    return;
}

#endif
