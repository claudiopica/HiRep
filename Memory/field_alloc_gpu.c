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
#include "gpu.h"


void free_gfield_gpu(suNg_field *field) {
  if (field->gpu_ptr!=NULL){
    cudaFree(field->gpu_ptr);
  }
}

void alloc_gfield_gpu(suNg_field *field) {
  
    cudaError_t err;
    
	if (field->gpu_ptr!=NULL) {
        free_gfield_gpu(field);
    }
    
    err = cudaMalloc(&field->gpu_ptr, 4*field->type->gsize*sizeof(suNg));
    error(err!=cudaSuccess,1,"gpu_alloc_gfield [gpu_field_alloc.c]",
          "Could not allocate GPU memory space for the gauge field");
    
}

void free_gfield_f_gpu(suNf_field *field) {
  if (field->gpu_ptr!=NULL){
    cudaFree(field->gpu_ptr);
  }
}

void alloc_gfield_f_gpu(suNf_field *field) {
  
    cudaError_t err;
    
	if (field->gpu_ptr!=NULL) {
	  free_gfield_f_gpu(field);
    }
    
    err = cudaMalloc(&field->gpu_ptr, 4*field->type->gsize*sizeof(suNf));
    error(err!=cudaSuccess,1,"gpu_alloc_gfield [gpu_field_alloc.c]",
          "Could not allocate GPU memory space for the gauge field");
    
}

void free_gfield_flt_gpu(suNg_field_flt *field) {
  if (field->gpu_ptr!=NULL){
    cudaFree(field->gpu_ptr);
  }
}

void alloc_gfield_flt_gpu(suNg_field_flt *field) {
  
    cudaError_t err;
    
	if (field->gpu_ptr!=NULL) {
        free_gfield_flt_gpu(field);
    }
    
    err = cudaMalloc(&field->gpu_ptr, 4*field->type->gsize*sizeof(suNg_flt));
    error(err!=cudaSuccess,1,"gpu_alloc_gfield [gpu_field_alloc.c]",
          "Could not allocate GPU memory space for the gauge field");
    
}

void free_gfield_f_flt_gpu(suNf_field_flt *field) {
  if (field->gpu_ptr!=NULL){
    cudaFree(field->gpu_ptr);
  }
}

void alloc_gfield_f_flt_gpu(suNf_field_flt *field) {
  
    cudaError_t err;
    
	if (field->gpu_ptr!=NULL) {
        free_gfield_f_flt_gpu(field);
    }
    
    err = cudaMalloc(&field->gpu_ptr, 4*field->type->gsize*sizeof(suNf_flt));
    error(err!=cudaSuccess,1,"gpu_alloc_gfield [gpu_field_alloc.c]",
          "Could not allocate GPU memory space for the gauge field");
    
}

void gfield_copy_to_gpu(suNg_field *field){
  suNg_field *tmp=alloc_gfield(field->type);
  gfield_togpuformat(tmp,field);
  cudaMemcpy(field->gpu_ptr,tmp->ptr,field->type->gsize*4*sizeof(suNg),cudaMemcpyHostToDevice);
  free_gfield(tmp);
}

void gfield_copy_from_gpu(suNg_field *field){
  suNg_field *tmp=alloc_gfield(field->type);
  cudaMemcpy(tmp->ptr,field->gpu_ptr,field->type->gsize*4*sizeof(suNg),cudaMemcpyDeviceToHost);
  gfield_tocpuformat(field,tmp);
  free_gfield(tmp);
}

void gfield_copy_to_gpu_f(suNf_field *field){
  suNf_field *tmp=alloc_gfield_f(field->type);
  gfield_togpuformat_f(tmp,field);
  cudaMemcpy(field->gpu_ptr,tmp->ptr,field->type->gsize*4*sizeof(suNf),cudaMemcpyHostToDevice);
  free_gfield_f(tmp);
}

void gfield_copy_from_gpu_f(suNf_field *field){
  suNf_field *tmp=alloc_gfield_f(field->type);
  cudaMemcpy(tmp->ptr,field->gpu_ptr,field->type->gsize*4*sizeof(suNf),cudaMemcpyDeviceToHost);
  gfield_tocpuformat_f(field,tmp);
  free_gfield_f(tmp);
}

void gfield_copy_to_gpu_flt(suNg_field_flt *field){
  cudaMemcpy(field->gpu_ptr,field->ptr,field->type->gsize*4*sizeof(suNg_flt),cudaMemcpyHostToDevice);
}

void gfield_copy_from_gpu_flt(suNg_field_flt *field){
  cudaMemcpy(field->ptr,field->gpu_ptr,field->type->gsize*4*sizeof(suNg_flt),cudaMemcpyDeviceToHost);
}

void gfield_copy_to_gpu_f_flt(suNf_field_flt *field){
  cudaMemcpy(field->gpu_ptr,field->ptr,field->type->gsize*4*sizeof(suNf_flt),cudaMemcpyHostToDevice);
}

void gfield_copy_from_gpu_f_flt(suNf_field_flt *field){
  cudaMemcpy(field->ptr,field->gpu_ptr,field->type->gsize*4*sizeof(suNf_flt),cudaMemcpyDeviceToHost);
}

