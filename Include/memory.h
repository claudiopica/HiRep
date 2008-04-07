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

#ifdef P4
#  define ALIGN 7
#else
#  define ALIGN 5 
#endif  

void *amalloc(size_t size,int p);
void afree(void *addr);

void free_field(void *u);
suNg* alloc_gfield();
suNf* alloc_gfield_f();
suNg_flt* alloc_gfield_flt();
suNf_flt* alloc_gfield_f_flt();
spinor_field* alloc_spinor_field_f(unsigned int n, spinor_descriptor *type);
spinor_field_flt* alloc_spinor_field_f_flt(unsigned int n, spinor_descriptor *type);
void free_spinor_field(spinor_field *s);
void free_spinor_field_flt(spinor_field_flt *s);
suNg_algebra_vector* alloc_momenta();

spinor_field* create_spinor_mask(spinor_field* s, spinor_descriptor* masktype);
void free_spinor_mask(spinor_field* s);


#endif
