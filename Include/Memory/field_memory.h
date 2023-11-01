/***************************************************************************\
* Copyright (c) 2008, 2022, Claudio Pica, Sofie Martins                     *
* All rights reserved.                                                      *
\***************************************************************************/

/// Headerfile for:
/// -  amalloc.c
/// -  memory.c
/// -  convert.cu

/**
 * @file field_memory.h
 * @brief Field allocation, GPU geometry conversion and 
 *        host-device/device-host copy functions
 */

#ifndef FIELD_MEMORY_H
#define FIELD_MEMORY_H

#include <stdlib.h>
#include "libhr_core.h"

#ifdef __cplusplus
extern "C" {
#endif

#ifdef P4
#define ALIGN 7
#else
#define ALIGN 8
#endif

//amalloc.c

/**
 * @brief Allocated memory aligned, because this improves bandwidth
 *
 * @param size			size to be allocated
 * @param p			    alignment
 */
void *amalloc(size_t size, int p);

/**
 * @brief Free memory that was allocated aligned using amalloc
 *
 * @param addr			Free this pointer
 */
void afree(void *addr);

//memory.c
//convert.cu

#define _FIELD_NAME_READABLE "Spinor field"
#define _FIELD_TYPE spinor_field
#define _IS_SPINOR_LIKE 1
#include "TMPL/field_memory.h.tmpl"

#define _FIELD_NAME_READABLE "Single precision spinor field"
#define _FIELD_TYPE spinor_field_flt
#define _IS_SPINOR_LIKE 1
#include "TMPL/field_memory.h.tmpl"

#define _FIELD_NAME_READABLE "Scalar field"
#define _FIELD_TYPE scalar_field
#define _IS_SPINOR_LIKE 1
#include "TMPL/field_memory.h.tmpl"

#define _FIELD_NAME_READABLE "Gauge field"
#define _FIELD_TYPE suNg_field
#define _IS_SPINOR_LIKE 0
#include "TMPL/field_memory.h.tmpl"

#define _FIELD_NAME_READABLE "Single precision gauge field"
#define _FIELD_TYPE suNg_field_flt
#define _IS_SPINOR_LIKE 0
#include "TMPL/field_memory.h.tmpl"

#define _FIELD_NAME_READABLE "Represented gauge field"
#define _FIELD_TYPE suNf_field
#define _IS_SPINOR_LIKE 0
#include "TMPL/field_memory.h.tmpl"

#define _FIELD_NAME_READABLE "Represented single precision gauge field"
#define _FIELD_TYPE suNf_field_flt
#define _IS_SPINOR_LIKE 0
#include "TMPL/field_memory.h.tmpl"

#define _FIELD_NAME_READABLE "SU(NG) scalar field"
#define _FIELD_TYPE suNg_scalar_field
#define _IS_SPINOR_LIKE 0
#include "TMPL/field_memory.h.tmpl"

#define _FIELD_NAME_READABLE "SU(NG) algebra vector field"
#define _FIELD_TYPE suNg_av_field
#define _IS_SPINOR_LIKE 0
#include "TMPL/field_memory.h.tmpl"

#define _FIELD_NAME_READABLE "Gauge transformation"
#define _FIELD_TYPE gtransf
#define _IS_SPINOR_LIKE 0
#include "TMPL/field_memory.h.tmpl"

#define _FIELD_NAME_READABLE "Clover ldl field"
#define _FIELD_TYPE ldl_field
#define _IS_SPINOR_LIKE 0
#include "TMPL/field_memory.h.tmpl"

#define _FIELD_NAME_READABLE "Clover term"
#define _FIELD_TYPE clover_term
#define _IS_SPINOR_LIKE 0
#include "TMPL/field_memory.h.tmpl"

#define _FIELD_NAME_READABLE "Clover force"
#define _FIELD_TYPE clover_force
#define _IS_SPINOR_LIKE 0
#include "TMPL/field_memory.h.tmpl"

#define _FIELD_NAME_READABLE "Staple field"
#define _FIELD_TYPE staple_field
#define _IS_SPINOR_LIKE 0
#include "TMPL/field_memory.h.tmpl"

#ifdef __cplusplus
}
#endif
#endif
