/***************************************************************************\
 * Copyright (c) 2008, Claudio Pica                                          *   
 * All rights reserved.                                                      * 
 \***************************************************************************/

#ifndef GPU_H
#define GPU_H
#include <stdio.h>

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
            fprintf( stderr, "cudaSafeCall() failed at %s:%i : %s\n",
                     file, line, cudaGetErrorString( err ) );
            exit( -1 );
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
            fprintf( stderr, "cudaCheckError() failed at %s:%i : %s.\n",
                     file, line, cudaGetErrorString( err ) );
            exit( -1 );
        }

        // More careful checking. However, this will affect performance.
        // Comment if not needed.
        err = cudaThreadSynchronize();
        if( cudaSuccess != err )
        {
            fprintf( stderr, "cudaCheckError() with sync failed at %s:%i : %s.\n",
                     file, line, cudaGetErrorString( err ) );
            exit( -1 );
        }
    } while ( 0 );

#pragma warning( pop )

#endif // CUDA_CHECK_ERROR

    return;
}

#endif
