/***************************************************************************\
 * Copyright (c) 2008, Claudio Pica                                          *   
 * All rights reserved.                                                      * 
 \***************************************************************************/

#ifndef GPU_H
#define GPU_H
#ifdef WITH_GPU

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <cuda.h>
#include <driver_types.h>
#include <cuda_runtime.h>
#include <cuda_runtime_api.h>

#include "logger.h"
#include "error.h"
#include "input_par.h"
#include "hr_complex.h"

/* GPU variables */
typedef struct _input_gpu {
  unsigned int gpuID;
  
  /* for the reading function */
  input_record_t read[2];
  
} input_gpu;

#define init_input_gpu(varname) \
{ \
.read={\
{"gpuID", "gpuID = %d", INT_T, &(varname)},\
{NULL, NULL, INT_T, NULL}\
}\
}
//#define init_input_gpu(varname) \
{ \
.gpuID=0,\
.read={\
{"gpuID", "gpuID = %d", INT_T, &(varname).gpuID},\
{NULL, NULL, INT_T, NULL}\
}\
}

double* alloc_double_sum_field(int n);
hr_complex* alloc_complex_sum_field(int n);

#define START_SP_ADDRESS_GPU(sf) ((sf)->gpu_ptr + (sf)->type->master_start[0])
#define START_GF_ADDRESS_GPU(gf) ((gf)->gpu_ptr + (gf)->type->master_start[0])

#define _GPU_FIELD_BLK(s,i) (((s)->gpu_ptr) + (s)->type->master_start[(i)])
#define _GPU_4FIELD_BLK(s,i) (((s)->gpu_ptr) + 4*(s)->type->master_start[(i)])

#define CudaSafeCall( err )     __cudaSafeCall( err, __FILE__, __LINE__ )
#define CudaCheckError()        __cudaCheckError( __FILE__, __LINE__ )
void __cudaSafeCall( cudaError_t err, const char *file, const int line );
void __cudaCheckError( const char *file, const int line );

#ifdef __cplusplus
}
#endif

#endif //WITH_GPU
#endif
