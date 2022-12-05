/***************************************************************************\
 * Copyright (c) 2008, Claudio Pica                                          *   
 * All rights reserved.                                                      * 
 \***************************************************************************/
/**
 * @file gpu.h
 * @brief Basic gpu imports and structs. Include this in files that define GPU logic.
 */
#ifndef GPU_H
#define GPU_H
#ifdef WITH_GPU

#include <stdio.h>
#include <cuda.h>
#include <driver_types.h>
#include <cuda_runtime.h>
#include <cuda_runtime_api.h>

#include "logger.h"
#include "error.h"
#include "input_par.h"
#include "hr_complex.h"

typedef struct 
{
  unsigned int gpuID;
  input_record_t read[2]; /* for the reading function */
} input_gpu;

#define init_input_gpu(varname) \
  { \
    .read={\
      {"gpuID", "gpuID = %d", INT_T, &(varname)},\
      {NULL, NULL, INT_T, NULL}\
    }\
  }

double* alloc_double_sum_field(int n);
hr_complex* alloc_complex_sum_field(int n);

#endif 
#endif
