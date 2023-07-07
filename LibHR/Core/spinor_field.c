/***************************************************************************\
* Copyright (c) 2023, Sofie Martins                                         *   
* All rights reserved.                                                      * 
\***************************************************************************/

#include "libhr_core.h"

#define _FIELD_NAME spinor_field_f
#define _FIELD_TYPE spinor_field
#define _FIELD_DIM 1
#include "TMPL/spinor_field.c.tmpl"

#define _FIELD_NAME spinor_field_f_flt
#define _FIELD_TYPE spinor_field_flt
#define _FIELD_DIM 1
#include "TMPL/spinor_field.c.tmpl"

#define _FIELD_NAME sfield
#define _FIELD_TYPE scalar_field
#define _FIELD_DIM 1
#include "TMPL/spinor_field.c.tmpl"

#define _FIELD_NAME gfield
#define _FIELD_TYPE suNg_field
#define _FIELD_DIM 4
#include "TMPL/spinor_field.c.tmpl"

#define _FIELD_NAME gfield_flt
#define _FIELD_TYPE suNg_field_flt
#define _FIELD_DIM 4
#include "TMPL/spinor_field.c.tmpl"

#define _FIELD_NAME gfield_f
#define _FIELD_TYPE suNf_field
#define _FIELD_DIM 4
#include "TMPL/spinor_field.c.tmpl"

#define _FIELD_NAME gfield_f_flt
#define _FIELD_TYPE suNf_field_flt
#define _FIELD_DIM 4
#include "TMPL/spinor_field.c.tmpl"

#define _FIELD_NAME suNg_scalar_field
#define _FIELD_TYPE suNg_scalar_field
#define _FIELD_DIM 1
#include "TMPL/spinor_field.c.tmpl"

#define _FIELD_NAME avfield
#define _FIELD_TYPE suNg_av_field
#define _FIELD_DIM 4
#include "TMPL/spinor_field.c.tmpl"

#define _FIELD_NAME gtransf
#define _FIELD_TYPE suNg_field
#define _FIELD_DIM 1
#include "TMPL/spinor_field.c.tmpl"

#define _FIELD_NAME clover_ldl
#define _FIELD_TYPE ldl_field
#define _FIELD_DIM 1
#include "TMPL/spinor_field.c.tmpl"

#define _FIELD_NAME clover_term
#define _FIELD_TYPE suNfc_field
#define _FIELD_DIM 4
#include "TMPL/spinor_field.c.tmpl"

#define _FIELD_NAME clover_force
#define _FIELD_TYPE suNf_field
#define _FIELD_DIM 6
#include "TMPL/spinor_field.c.tmpl"

#define _FIELD_NAME staple_field
#define _FIELD_TYPE suNg_field
#define _FIELD_DIM 3
#include "TMPL/spinor_field.c.tmpl"