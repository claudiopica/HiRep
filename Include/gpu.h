/***************************************************************************\
 * Copyright (c) 2008, Claudio Pica                                          *   
 * All rights reserved.                                                      * 
 \***************************************************************************/

#ifndef GPU_H
#define GPU_H
#ifdef WITH_GPU

#include <stdio.h>
#include <cuda.h>
#include <driver_types.h>
#include "error.h"



#define START_SP_ADDRESS_GPU(sf) (sf->gpu_ptr + sf->type->master_start[0])

// Enable this for error checking
#define CUDA_CHECK_ERROR

#define CudaSafeCall( err )     __cudaSafeCall( err, __FILE__, __LINE__ )
#define CudaCheckError()        __cudaCheckError( __FILE__, __LINE__ )

inline void __cudaSafeCall( cudaError err, const char *file, const int line )
{
#ifdef CUDA_CHECK_ERROR

#pragma warning( push )
#pragma warning( disable: 4127 ) // Prevent warning on do-while(0);

    do
    {
        if ( cudaSuccess != err )
        {
            lprintf("CUDA",0,"cudaSafeCall() failed at %s:%i\n",
                     file, line);
            error((cudaSuccess != err),1,"CudaSafeCall", 
                  cudaGetErrorString( err ));
        }
    } while ( 0 );

#pragma warning( pop )

#endif  // CUDA_CHECK_ERROR

    return;
}

inline void __cudaCheckError( const char *file, const int line )
{
#ifdef CUDA_CHECK_ERROR

#pragma warning( push )
#pragma warning( disable: 4127 ) // Prevent warning on do-while(0);

    do
    {
        cudaError_t err = cudaGetLastError();
        if ( cudaSuccess != err )
        {
            lprintf("CUDA",0,"cudaCheckError() failed at %s:%i\n",
                    file, line);

            error((cudaSuccess != err),1,"CudaCheckError", 
                  cudaGetErrorString( err ));
            
        }

        // More careful checking. However, this will affect performance.
        // Comment if not needed.
        err = cudaThreadSynchronize();
        if( cudaSuccess != err )
        {
            lprintf("CUDA",0,"cudaCheckError() with sync failed at %s:%i\n",
                    file, line);
            
            error((cudaSuccess != err),1,"CudaCheckError with sync", 
                  cudaGetErrorString( err ));
        }
    } while ( 0 );

#pragma warning( pop )

#endif // CUDA_CHECK_ERROR

    return;
}

#undef CUDA_CHECK_ERROR

#endif //WITH_GPU
#endif
