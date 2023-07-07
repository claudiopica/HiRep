/***************************************************************************\
* Copyright (c) 2022, Claudio Pica, Sofie Martins                           *   
* All rights reserved.                                                      * 
\***************************************************************************/

/// Headerfile for:
/// - communications.c
/// - sync_to_buffer.cu

/**
 * @file communications.h
 * @brief Communications to send and receive sites from other nodes using MPI with GPU
 */

#ifndef COMMUNICATIONS_REDUCED_H
#define COMMUNICATIONS_REDUCED_H

#include "spinor_field.h"
#include "geometry.h"
#include "Utils/generics.h"

#ifdef __cplusplus
extern "C" {
#endif

#define _FIELD_NAME_READABLE "Spinor field"
#define _FIELD_NAME spinor_field_f
#define _FIELD_TYPE spinor_field
#include "TMPL/communications_reduced.h.tmpl"

#define _FIELD_NAME_READABLE "Single precision spinor field"
#define _FIELD_NAME spinor_field_f_flt
#define _FIELD_TYPE spinor_field_flt
#include "TMPL/communications_reduced.h.tmpl"

#ifdef __cplusplus
}
#endif
#endif