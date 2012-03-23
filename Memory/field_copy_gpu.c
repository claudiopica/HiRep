/***************************************************************************\
 * Copyright (c) 2008, Claudio Pica                                          *   
 * All rights reserved.                                                      * 
 \***************************************************************************/

/*******************************************************************************
 *
 * File field_copy_gpu.c
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
#ifdef WITH_GPU
/* GAUGE FIELDS */
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
  suNf_field_flt *tmp=alloc_gfield_f_flt(field->type);
  gfield_togpuformat_f_flt(tmp,field);
  cudaMemcpy(field->gpu_ptr,tmp->ptr,field->type->gsize*4*sizeof(suNf_flt),cudaMemcpyHostToDevice);
	free_gfield_f_flt(tmp);
}

void gfield_copy_from_gpu_f_flt(suNf_field_flt *field){
  suNf_field_flt *tmp=alloc_gfield_f_flt(field->type);
  cudaMemcpy(tmp->ptr,field->gpu_ptr,field->type->gsize*4*sizeof(suNf_flt),cudaMemcpyDeviceToHost);
  gfield_tocpuformat_f_flt(field,tmp);
  free_gfield_f_flt(tmp);
}

/* SPINORS */
void spinor_field_copy_to_gpu_f(spinor_field *field){
  spinor_field *tmp = alloc_spinor_field_f(1, field->type);
  spinor_field_togpuformat(tmp, field);
  cudaMemcpy(field->gpu_ptr,tmp->ptr,field->type->gsize*sizeof(suNf_spinor),cudaMemcpyHostToDevice);
  free_spinor_field_f(tmp);
}

void spinor_field_copy_from_gpu_f(spinor_field *field){
  spinor_field *tmp = alloc_spinor_field_f(1, field->type);
  cudaMemcpy(tmp->ptr,field->gpu_ptr,field->type->gsize*sizeof(suNf_spinor),cudaMemcpyDeviceToHost);
  spinor_field_tocpuformat(field,tmp);
  free_spinor_field_f(tmp);
}

void spinor_field_copy_to_gpu_f_flt(spinor_field_flt *field){
  spinor_field_flt *tmp = alloc_spinor_field_f_flt(1, field->type);
  spinor_field_togpuformat_flt(tmp, field);
  cudaMemcpy(field->gpu_ptr,tmp->ptr,field->type->gsize*sizeof(suNf_spinor_flt),cudaMemcpyHostToDevice);
  free_spinor_field_f_flt(tmp);
}

void spinor_field_copy_from_gpu_f_flt(spinor_field_flt *field){
  spinor_field_flt *tmp = alloc_spinor_field_f_flt(1, field->type);
  cudaMemcpy(tmp->ptr,field->gpu_ptr,field->type->gsize*sizeof(suNf_spinor_flt),cudaMemcpyDeviceToHost);
  spinor_field_tocpuformat_flt(field,tmp);
  free_spinor_field_f_flt(tmp);
}

/* SCALAR_FIELDS */
void sfield_copy_to_gpu(scalar_field *field){
  cudaMemcpy(field->gpu_ptr,field->ptr,field->type->gsize*sizeof(double),cudaMemcpyHostToDevice);
}
void sfield_copy_from_gpu(scalar_field *field){
  cudaMemcpy(field->ptr,field->gpu_ptr,field->type->gsize*sizeof(double),cudaMemcpyDeviceToHost);
}

//avfield
void suNg_av_field_copy_to_gpu(suNg_av_field *field){
  suNg_av_field* tmp = alloc_avfield(field->type);
  avfield_togpuformat(tmp,field);
  cudaMemcpy(field->gpu_ptr,tmp->ptr,field->type->gsize*4*sizeof(suNg_algebra_vector),cudaMemcpyHostToDevice);
  free_avfield(tmp);
}

void suNg_av_field_copy_from_gpu(suNg_av_field *field){
  suNg_av_field* tmp = alloc_avfield(field->type);
  cudaMemcpy(tmp->ptr,field->gpu_ptr,field->type->gsize*4*sizeof(suNg_algebra_vector),cudaMemcpyDeviceToHost);
  avfield_tocpuformat(field,tmp);
  free_avfield(tmp);
}

#endif
