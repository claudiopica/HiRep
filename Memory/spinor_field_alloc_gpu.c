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
