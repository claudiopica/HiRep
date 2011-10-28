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

void free_spinor_field_gpu(spinor_field *field) {
    if (field[0].gpu_ptr!=NULL) {
        cudaFree(field[0].gpu_ptr);
        field[0].gpu_ptr=NULL;
    }
}


void alloc_spinor_field_f_gpu(unsigned int n, spinor_field *field) {

    cudaError_t err;
    unsigned int i;
    suNf_spinor *p
    
// DO SOME TESTS, e.g. deallocate memory if gpu_ptr != NULL
    
    err = cudaMalloc(&p, n*type->gsize*sizeof(suNf_spinor));
    error(err!=cudaSuccess,1,"alloc_spinor_field_f_gpu [spinor_field_alloc_gpu.c]",
          "Could not allocate GPU memory space for the spinor fields");
    
    for(i=0; i<n; ++i) {
        field[i].gpu_ptr=p+i*type->gsize;
    }
    
}

void free_spinor_field_flt_gpu(spinor_field_flt *field) {
    if (field[0].gpu_ptr!=NULL) {
        cudaFree(field[0].gpu_ptr);
        field[0].gpu_ptr=NULL;
    }
}


void alloc_spinor_field_f_flt_gpu(unsigned int n, spinor_field_flt *field) {
    
    cudaError_t err;
    unsigned int i;
    suNf_spinor_flt *p
    
    // DO SOME TESTS, e.g. deallocate memory if gpu_ptr != NULL
    
    err = cudaMalloc(&p, n*type->gsize*sizeof(suNf_spinor_flt));
    error(err!=cudaSuccess,1,"alloc_spinor_field_f_gpu [spinor_field_alloc_gpu.c]",
          "Could not allocate GPU memory space for the spinor fields");
    
    for(i=0; i<n; ++i) {
        field[i].gpu_ptr=p+i*type->gsize;
    }
    
}