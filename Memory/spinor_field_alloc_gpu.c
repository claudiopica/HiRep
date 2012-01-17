/***************************************************************************\
 * Copyright (c) 2008, Claudio Pica                                          *   
 * All rights reserved.                                                      * 
 \***************************************************************************/

/*******************************************************************************
 *
 * File spinor_field_alloc_gpu.c
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
#include "gpu.h"

void free_spinor_field_gpu(spinor_field *field) {
  if (field[0].gpu_ptr!=NULL) {
    cudaFree(field[0].gpu_ptr);
    field[0].gpu_ptr=NULL;
  }
}

void alloc_spinor_field_f_gpu(unsigned int n, spinor_field *field) {

    cudaError_t err;
    unsigned int i;
    suNf_spinor *p;
    
// DO SOME TESTS, e.g. deallocate memory if gpu_ptr != NULL
    
    err = cudaMalloc((void**) &p, n*field->type->gsize*sizeof(suNf_spinor));
    error(err!=cudaSuccess,1,"alloc_spinor_field_f_gpu [spinor_field_alloc_gpu.c]",
          "Could not allocate GPU memory space for the spinor fields");
    
    for(i=0; i<n; ++i) {
      field[i].gpu_ptr=p+i*field->type->gsize;
    }
    
}

void spinor_field_copy_to_gpu_f(spinor_field *field){
    spinor_field *tmp = alloc_spinor_field_f(1, field->type);
    spinor_field_togpuformat(tmp, field);
    cudaMemcpy(field->gpu_ptr,tmp->ptr,field->type->gsize*sizeof(suNf_spinor),cudaMemcpyHostToDevice);
    free_spinor_field(tmp);
}

void spinor_field_copy_from_gpu_f(spinor_field *field){
    spinor_field *tmp = alloc_spinor_field_f(1, field->type);
    cudaMemcpy(tmp->ptr,field->gpu_ptr,field->type->gsize*sizeof(suNf_spinor),cudaMemcpyDeviceToHost);
    spinor_field_tocpuformat(field,tmp);
    free_spinor_field(tmp);
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
    suNf_spinor_flt *p;
    
    // DO SOME TESTS, e.g. deallocate memory if gpu_ptr != NULL
    
    err = cudaMalloc((void **) &p, n*field->type->gsize*sizeof(suNf_spinor_flt));
    error(err!=cudaSuccess,1,"alloc_spinor_field_f_gpu [spinor_field_alloc_gpu.c]",
          "Could not allocate GPU memory space for the spinor fields");
    
    for(i=0; i<n; ++i) {
        field[i].gpu_ptr=p+i*field->type->gsize;
    }
    
}

void spinor_field_copy_to_gpu_f_flt(spinor_field_flt *field){
  spinor_field_flt *tmp = alloc_spinor_field_f_flt(1, field->type);
  spinor_field_togpuformat_flt(tmp, field);
  cudaMemcpy(field->gpu_ptr,tmp->ptr,field->type->gsize*sizeof(suNf_spinor_flt),cudaMemcpyHostToDevice);
  free_spinor_field_flt(tmp);
}

void spinor_field_copy_from_gpu_f_flt(spinor_field_flt *field){
  spinor_field_flt *tmp = alloc_spinor_field_f_flt(1, field->type);
  cudaMemcpy(tmp->ptr,field->gpu_ptr,field->type->gsize*sizeof(suNf_spinor_flt),cudaMemcpyDeviceToHost);
  spinor_field_tocpuformat_flt(field,tmp);
  free_spinor_field_flt(tmp);
}
