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
#include <cuda_runtime.h>
#include <cuda_runtime_api.h>

#include "logger.h"
#include "error.h"



#define START_SP_ADDRESS_GPU(sf) (sf->gpu_ptr + sf->type->master_start[0])

#endif //WITH_GPU
#endif
