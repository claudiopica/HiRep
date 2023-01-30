/***************************************************************************\
* Copyright (c) 2022, Sofie Martins                                         *   
* All rights reserved.                                                      * 
\***************************************************************************/

#include "geometry.h"
#include "libhr_core.h"

#define _FIELD_NAME spinor_field_f
#define _FIELD_TYPE spinor_field
#define _SITE_TYPE suNf_spinor
#define _FIELD_DIM 1
#define _GEOM_TYPE spinor
#define _COMPLEX hr_complex
#define _REAL double
#include "TMPL/sync_to_buffer.cu.tmpl"

#define _FIELD_NAME spinor_field_f_flt
#define _FIELD_TYPE spinor_field_flt
#define _SITE_TYPE suNf_spinor_flt
#define _FIELD_DIM 1
#define _GEOM_TYPE spinor
#define _COMPLEX hr_complex_flt
#define _REAL float
#include "TMPL/sync_to_buffer.cu.tmpl"

#define _FIELD_NAME sfield
#define _FIELD_TYPE scalar_field
#define _SITE_TYPE double
#define _FIELD_DIM 1
#define _GEOM_TYPE spinor
#define _COMPLEX hr_complex
#define _REAL double
#include "TMPL/sync_to_buffer.cu.tmpl"

#define _FIELD_NAME gfield
#define _FIELD_TYPE suNg_field
#define _SITE_TYPE suNg
#define _FIELD_DIM 4
#define _GEOM_TYPE gauge
#define _COMPLEX hr_complex
#define _REAL double
#include "TMPL/sync_to_buffer.cu.tmpl"

#define _FIELD_NAME gfield_flt
#define _FIELD_TYPE suNg_field_flt
#define _SITE_TYPE suNg_flt
#define _FIELD_DIM 4
#define _GEOM_TYPE gauge
#define _COMPLEX hr_complex_flt
#define _REAL float
#include "TMPL/sync_to_buffer.cu.tmpl"

#define _FIELD_NAME gfield_f
#define _FIELD_TYPE suNf_field
#define _SITE_TYPE suNf
#define _FIELD_DIM 4
#define _GEOM_TYPE gauge
#define _COMPLEX hr_complex
#define _REAL double
#include "TMPL/sync_to_buffer.cu.tmpl"

#define _FIELD_NAME gfield_f_flt
#define _FIELD_TYPE suNf_field_flt
#define _SITE_TYPE suNf_flt
#define _FIELD_DIM 4
#define _GEOM_TYPE gauge
#define _COMPLEX hr_complex_flt
#define _REAL float
#include "TMPL/sync_to_buffer.cu.tmpl"

#define _FIELD_NAME suNg_scalar_field
#define _FIELD_TYPE suNg_scalar_field
#define _SITE_TYPE suNg_vector
#define _FIELD_DIM 1
#define _GEOM_TYPE gauge
#define _COMPLEX hr_complex
#define _REAL double
#include "TMPL/sync_to_buffer.cu.tmpl"

#define _FIELD_NAME avfield
#define _FIELD_TYPE suNg_av_field
#define _SITE_TYPE suNg_algebra_vector
#define _FIELD_DIM 4
#define _GEOM_TYPE gauge
#define _COMPLEX hr_complex
#define _REAL double
#include "TMPL/sync_to_buffer.cu.tmpl"

#define _FIELD_NAME gtransf
#define _FIELD_TYPE suNg_field
#define _SITE_TYPE suNg
#define _FIELD_DIM 1
#define _GEOM_TYPE gauge
#define _COMPLEX hr_complex
#define _REAL double
#include "TMPL/sync_to_buffer.cu.tmpl"

#define _FIELD_NAME clover_ldl
#define _FIELD_TYPE ldl_field
#define _SITE_TYPE ldl_t
#define _FIELD_DIM 1
#define _GEOM_TYPE gauge
#define _COMPLEX hr_complex
#define _REAL double
#include "TMPL/sync_to_buffer.cu.tmpl"

#define _FIELD_NAME clover_term
#define _FIELD_TYPE suNfc_field
#define _SITE_TYPE suNfc
#define _FIELD_DIM 4
#define _GEOM_TYPE gauge
#define _COMPLEX hr_complex
#define _REAL double
#include "TMPL/sync_to_buffer.cu.tmpl"

#define _FIELD_NAME clover_force
#define _FIELD_TYPE suNf_field
#define _SITE_TYPE suNf
#define _FIELD_DIM 6
#define _GEOM_TYPE gauge
#define _COMPLEX hr_complex
#define _REAL double
#include "TMPL/sync_to_buffer.cu.tmpl"

#define _FIELD_NAME staple_field
#define _FIELD_TYPE suNg_field
#define _SITE_TYPE suNg
#define _FIELD_DIM 3
#define _GEOM_TYPE gauge
#define _COMPLEX hr_complex
#define _REAL double
#include "TMPL/sync_to_buffer.cu.tmpl"