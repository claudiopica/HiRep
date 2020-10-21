/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *   
* All rights reserved.                                                      * 
\***************************************************************************/

/*******************************************************************************
*
* File error.h
* 
* Error handling functions
*
*******************************************************************************/

#ifndef ERROR_H
#define ERROR_H

#ifdef __cplusplus
extern "C" {
#endif

void error(int test,int no, char *name, char *text);

#ifdef WITH_GPU
#include "gpu.h"

#define CudaSafeCall( err )     __cudaSafeCall( err, __FILE__, __LINE__ )
#define CudaCheckError()        __cudaCheckError( __FILE__, __LINE__ )
void __cudaSafeCall( cudaError_t err, const char *file, const int line );
void __cudaCheckError( const char *file, const int line );

#endif /* WITH_GPU */

#ifdef __cplusplus
}
#endif
#endif
