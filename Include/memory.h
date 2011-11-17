/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *   
* All rights reserved.                                                      * 
\***************************************************************************/

/*******************************************************************************
*
* File memory.h
* 
* Memory handling functions
*
*******************************************************************************/


#ifndef MEMORY_H
#define MEMORY_H

#include <stdlib.h>
#include "suN.h"
#include "spinor_field.h"

#ifdef P4
#  define ALIGN 7
#else
#  define ALIGN 5 
#endif  

#ifdef __cplusplus
extern "C" {
#endif
void *amalloc(size_t size,int p);
void afree(void *addr);

void free_gfield(suNg_field *u);
suNg_field* alloc_gfield(geometry_descriptor* type);
void free_gfield_f(suNf_field *u);
suNf_field* alloc_gfield_f(geometry_descriptor* type);
void free_gfield_flt(suNg_field_flt *u);
suNg_field_flt* alloc_gfield_flt(geometry_descriptor* type);
void free_gfield_f_flt(suNf_field_flt *u);
suNf_field_flt* alloc_gfield_f_flt(geometry_descriptor* type);
void free_gtransf(suNg_field *u);
suNg_field* alloc_gtransf(geometry_descriptor* type);

spinor_field* alloc_spinor_field_f(unsigned int n, geometry_descriptor *type);
spinor_field_flt* alloc_spinor_field_f_flt(unsigned int n, geometry_descriptor *type);
void free_spinor_field(spinor_field *s);
void free_spinor_field_flt(spinor_field_flt *s);

#ifdef WITH_GPU

void free_gfield_gpu(suNg_field *field);
void alloc_gfield_gpu(suNg_field *field);
void free_gfield_f_gpu(suNf_field *field);
void alloc_gfield_f_gpu(suNf_field *field);
void free_gfield_flt_gpu(suNg_field_flt *field);
void alloc_gfield_flt_gpu(suNg_field_flt *field);
void free_gfield_f_flt_gpu(suNf_field_flt *field);
void alloc_gfield_f_flt_gpu(suNf_field_flt *field);

void free_spinor_field_gpu(spinor_field *field);
void free_spinor_field_flt_gpu(spinor_field_flt *field);
void alloc_spinor_field_f_gpu(unsigned int n, spinor_field *field);
void alloc_spinor_field_f_flt_gpu(unsigned int n, spinor_field_flt *field);
void spinor_field_copy_to_gpu_f(spinor_field *field);
void spinor_field_copy_from_gpu_f(spinor_field *field);
void spinor_field_copy_to_gpu_f_flt(spinor_field_flt *field);
void spinor_field_copy_from_gpu_f_flt(spinor_field_flt *field);

void gfield_copy_to_gpu(suNg_field *field);
void gfield_copy_from_gpu(suNg_field *field);
void gfield_copy_to_gpu_f(suNf_field *field);
void gfield_copy_from_gpu_f(suNf_field *field);
void gfield_copy_to_gpu_flt(suNg_field_flt *field);
void gfield_copy_from_gpu_flt(suNg_field_flt *field);
void gfield_copy_to_gpu_f_flt(suNf_field_flt *field);
void gfield_copy_from_gpu_f_flt(suNf_field_flt *field);
    
void spinor_field_togpuformat(spinor_field *out, spinor_field* in);
void spinor_field_tocpuformat(spinor_field *out, spinor_field* in);

#endif

void free_avfield(suNg_av_field *u);
suNg_av_field *alloc_avfield(geometry_descriptor* type);
void free_sfield(scalar_field *u);
scalar_field *alloc_sfield(geometry_descriptor* type);

spinor_field* create_spinor_mask(spinor_field* s, geometry_descriptor* masktype);
void free_spinor_mask(spinor_field* s);
#ifdef __cplusplus
}
#endif


#endif
