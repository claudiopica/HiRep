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
  suNg_field_flt *tmp=alloc_gfield_flt(field->type);
  gfield_togpuformat_flt(tmp,field);
  cudaMemcpy(field->gpu_ptr,tmp->ptr,field->type->gsize*4*sizeof(suNg_flt),cudaMemcpyHostToDevice);
	free_gfield_flt(tmp);
}

void gfield_copy_from_gpu_flt(suNg_field_flt *field){
  suNg_field_flt *tmp=alloc_gfield_flt(field->type);
  cudaMemcpy(tmp->ptr,field->gpu_ptr,field->type->gsize*4*sizeof(suNg_flt),cudaMemcpyDeviceToHost);
  gfield_tocpuformat_flt(field,tmp);
  free_gfield_flt(tmp);
}

void gfield_copy_to_gpu_f_flt(suNf_field_flt *field){
  suNf_field_flt *tmp=alloc_gfield_flt(field->type);
  gfield_togpuformat_f_flt(tmp,field);
  cudaMemcpy(field->gpu_ptr,tmp->ptr,field->type->gsize*4*sizeof(suNf_flt),cudaMemcpyHostToDevice);
	free_gfield_f_flt(tmp);
}

void gfield_copy_from_gpu_f_flt(suNf_field_flt *field){
  suNf_field_flt *tmp=alloc_gfield_flt(field->type);
  cudaMemcpy(tmp->ptr,field->gpu_ptr,field->type->gsize*4*sizeof(suNf_flt),cudaMemcpyDeviceToHost);
  gfield_tocpuformat_f_flt(field,tmp);
  free_gfield_f_flt(tmp);
}

