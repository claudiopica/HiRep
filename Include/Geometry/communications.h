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

#ifndef COMMUNICATIONS_H
#define COMMUNICATIONS_H

#include "spinor_field.h"
#include "geometry.h"
#include "Utils/generics.h"

#ifdef __cplusplus
extern "C" {
#endif

void probe_mpi(void);

#define _FIELD_NAME_READABLE "Spinor field"
#define _FIELD_TYPE spinor_field
#include "TMPL/communications.h.tmpl"

#define _FIELD_NAME_READABLE "Single precision spinor field"
#define _FIELD_TYPE spinor_field_flt
#include "TMPL/communications.h.tmpl"

#define _FIELD_NAME_READABLE "Scalar field"
#define _FIELD_TYPE scalar_field
#include "TMPL/communications.h.tmpl"

#define _FIELD_NAME_READABLE "Gauge field"
#define _FIELD_TYPE suNg_field
#include "TMPL/communications.h.tmpl"

#define _FIELD_NAME_READABLE "Single precision gauge field"
#define _FIELD_TYPE suNg_field_flt
#include "TMPL/communications.h.tmpl"

#define _FIELD_NAME_READABLE "Represented gauge field"
#define _FIELD_TYPE suNf_field
#include "TMPL/communications.h.tmpl"

#define _FIELD_NAME_READABLE "Represented single precision gauge field"
#define _FIELD_TYPE suNf_field_flt
#include "TMPL/communications.h.tmpl"

#define _FIELD_NAME_READABLE "SU(NG) scalar field"
#define _FIELD_TYPE suNg_scalar_field
#include "TMPL/communications.h.tmpl"

#define _FIELD_NAME_READABLE "SU(NG) algebra vector field"
#define _FIELD_TYPE suNg_av_field
#include "TMPL/communications.h.tmpl"

#define _FIELD_NAME_READABLE "Gauge transformation"
#define _FIELD_TYPE gtransf
#include "TMPL/communications.h.tmpl"

#define _FIELD_NAME_READABLE "Clover ldl field"
#define _FIELD_TYPE ldl_field
#include "TMPL/communications.h.tmpl"

#define _FIELD_NAME_READABLE "Clover term"
#define _FIELD_TYPE clover_term
#include "TMPL/communications.h.tmpl"

#define _FIELD_NAME_READABLE "Clover force"
#define _FIELD_TYPE clover_force
#include "TMPL/communications.h.tmpl"

#define _FIELD_NAME_READABLE "Staple field"
#define _FIELD_TYPE staple_field
#include "TMPL/communications.h.tmpl"

#ifdef __cplusplus
}
#endif
#endif