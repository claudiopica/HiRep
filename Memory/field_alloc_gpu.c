/***************************************************************************\
 * Copyright (c) 2008, Claudio Pica                                          *   
 * All rights reserved.                                                      * 
 \***************************************************************************/

/*******************************************************************************
 *
 * File gpu_field_alloc.c
 *
 * Functions for fields allocation on GPUs
 *
 *******************************************************************************/

#include <stdlib.h>
#include "suN.h"
#include "error.h"
#include "memory.h"
#include "global.h"
#include "spinor_field.h"
#include "geometry.h"
#ifdef WITH_MPI
#include <mpi.h>
#endif

void gpu_free_gfield(suNg_field *field) {

}


void gpu_alloc_gfield(suNg_field *field) {

    cudaError_t err;
    
	if (field->gpu_ptr!=NULL) {
        gpu_free_gfield(field);
    }
    
    err = cudaMalloc(&field->gpu_prt, 4*type->gsize*sizeof(suNg));
    error(err!=cudaSuccess,1,"gpu_alloc_gfield [gpu_field_alloc.c]",
          "Could not allocate GPU memory space for the gauge field");
    
}