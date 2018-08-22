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

void *amalloc(size_t size,int p);
void afree(void *addr);

void free_gfield(suNg_field *u);
suNg_field* alloc_gfield(geometry_descriptor* type);
void free_scalar_field(suNg_scalar_field *u);
suNg_scalar_field* alloc_scalar_field(geometry_descriptor* type);
void free_gfield_f(suNf_field *u);
suNf_field* alloc_gfield_f(geometry_descriptor* type);
void free_gfield_flt(suNg_field_flt *u);
suNg_field_flt* alloc_gfield_flt(geometry_descriptor* type);
void free_gfield_f_flt(suNf_field_flt *u);
suNf_field_flt* alloc_gfield_f_flt(geometry_descriptor* type);
void free_gtransf(suNg_field *u);
suNg_field* alloc_gtransf(geometry_descriptor* type);
void free_avfield(suNg_av_field *u);
suNg_av_field *alloc_avfield(geometry_descriptor* type);

void free_clover_ldl(ldl_field *u);
ldl_field *alloc_clover_ldl(geometry_descriptor* type);
void free_clover_term(suNfc_field *u);
suNfc_field* alloc_clover_term(geometry_descriptor* type);
void free_clover_force(suNf_field *u);
suNf_field* alloc_clover_force(geometry_descriptor* type);

void free_spinor_field_f(spinor_field *s);
spinor_field* alloc_spinor_field_f(unsigned int n, geometry_descriptor *type);
void free_spinor_field_f_flt(spinor_field_flt *s);
spinor_field_flt* alloc_spinor_field_f_flt(unsigned int n, geometry_descriptor *type);
void free_sfield(scalar_field *s);
scalar_field *alloc_sfield(unsigned int n, geometry_descriptor* type);


spinor_field* create_spinor_mask(spinor_field* s, geometry_descriptor* masktype);
void free_spinor_mask(spinor_field* s);


#endif
