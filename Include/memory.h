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
#define ALIGN 8
#endif

#ifdef __cplusplus
extern "C" {
#endif

void *amalloc(size_t size,int p);
void afree(void *addr);


#define _DECLARE_MEMORY_FUNC_GAUGE(_name, _field_type, _site_type, _size)                                               \
        void copy_to_gpu_##_name(_field_type*);                                                                   \
        void copy_from_gpu_##_name(_field_type*);                                                                 \
        void to_gpu_format_##_name(_field_type*, _field_type*);                                                   \
        void to_cpu_format_##_name(_field_type*, _field_type*);                                                   \
    void free_##_name(_field_type*);                                                                              \
    _field_type *alloc_##_name(geometry_descriptor*);

#define _DECLARE_MEMORY_FUNC_MATTER(_name, _field_type, _site_type, _size)                                               \
        void copy_to_gpu_##_name(_field_type*);                                                                   \
        void copy_from_gpu_##_name(_field_type*);                                                                 \
        void to_gpu_format_##_name(_field_type*, _field_type*);                                                   \
        void to_cpu_format_##_name(_field_type*, _field_type*);                                                   \
    void free_##_name(_field_type*);                                                                              \
    _field_type *alloc_##_name(unsigned int, geometry_descriptor*);

// Gauge Fields
_DECLARE_MEMORY_FUNC_GAUGE(gfield, suNg_field, suNg, 4);
_DECLARE_MEMORY_FUNC_GAUGE(gfield_flt, suNg_field_flt, suNg_flt, 4);
_DECLARE_MEMORY_FUNC_GAUGE(gfield_f, suNf_field, suNf, 4);
_DECLARE_MEMORY_FUNC_GAUGE(gfield_f_flt, suNf_field_flt, suNf_flt, 4);
_DECLARE_MEMORY_FUNC_GAUGE(scalar_field, suNg_scalar_field, suNg_vector, 1);
_DECLARE_MEMORY_FUNC_GAUGE(avfield, suNg_av_field, suNg_algebra_vector, 4);
_DECLARE_MEMORY_FUNC_GAUGE(gtransf, suNg_field, suNg, 1);
_DECLARE_MEMORY_FUNC_GAUGE(clover_ldl, ldl_field, ldl_t, 1);
_DECLARE_MEMORY_FUNC_GAUGE(clover_term, suNfc_field, suNfc, 4);
_DECLARE_MEMORY_FUNC_GAUGE(clover_force, suNf_field, suNf, 6);

// Matter Fields
_DECLARE_MEMORY_FUNC_MATTER(spinor_field_f, spinor_field, suNf_spinor, 1);
_DECLARE_MEMORY_FUNC_MATTER(spinor_field_f_flt, spinor_field_flt, suNf_spinor_flt, 1);
_DECLARE_MEMORY_FUNC_MATTER(sfield, scalar_field, double, 1);

void spinor_field_togpuformat(spinor_field*, spinor_field*);

#undef _DECLARE_MEMORY_FUNC

#ifdef __cplusplus
}
#endif

#endif
