/***************************************************************************\
* Copyright (c) 2023, Sofie Martins                                         *   
* All rights reserved.                                                      * 
\***************************************************************************/

#ifdef WITH_GPU

#include "memory.h"
#include "geometry.h"
#include "libhr_core.h"

/* Spinor-like fields */
#define _GEOM_TYPE spinor
#define _IS_SPINOR_LIKE 1
#define _n n

#define _FIELD_NAME spinor_field_f
#define _FIELD_TYPE spinor_field
#define _SITE_TYPE suNf_spinor
#define _FIELD_DIM 1
#include "TMPL/field_convert.c.tmpl"
#include "TMPL/field_device_transfer.c.tmpl"
#include "TMPL/field_alloc.c.tmpl"

#define _FIELD_NAME spinor_field_f_flt
#define _FIELD_TYPE spinor_field_flt
#define _SITE_TYPE suNf_spinor_flt
#define _FIELD_DIM 1
#include "TMPL/field_convert.c.tmpl"
#include "TMPL/field_device_transfer.c.tmpl"
#include "TMPL/field_alloc.c.tmpl"

#define _FIELD_NAME sfield
#define _FIELD_TYPE scalar_field
#define _SITE_TYPE double
#define _FIELD_DIM 1
#include "TMPL/field_convert.c.tmpl"
#include "TMPL/field_device_transfer.c.tmpl"
#include "TMPL/field_alloc.c.tmpl"

#undef _GEOM_TYPE
#undef _IS_SPINOR_LIKE
#undef _n

/* Gauge fields */
#define _GEOM_TYPE gauge
#define _IS_SPINOR_LIKE 0
#define _n 1

#define _FIELD_NAME gfield
#define _FIELD_TYPE suNg_field
#define _SITE_TYPE suNg
#define _FIELD_DIM 4
#include "TMPL/field_convert.c.tmpl"
#include "TMPL/field_device_transfer.c.tmpl"
#include "TMPL/field_alloc.c.tmpl"

#define _FIELD_NAME gfield_flt
#define _FIELD_TYPE suNg_field_flt
#define _SITE_TYPE suNg_flt
#define _FIELD_DIM 4
#include "TMPL/field_convert.c.tmpl"
#include "TMPL/field_device_transfer.c.tmpl"
#include "TMPL/field_alloc.c.tmpl"

#define _FIELD_NAME gfield_f
#define _FIELD_TYPE suNf_field
#define _SITE_TYPE suNf
#define _FIELD_DIM 4
#include "TMPL/field_convert.c.tmpl"
#include "TMPL/field_device_transfer.c.tmpl"
#include "TMPL/field_alloc.c.tmpl"

#define _FIELD_NAME gfield_f_flt
#define _FIELD_TYPE suNf_field_flt
#define _SITE_TYPE suNf_flt
#define _FIELD_DIM 4
#include "TMPL/field_convert.c.tmpl"
#include "TMPL/field_device_transfer.c.tmpl"
#include "TMPL/field_alloc.c.tmpl"

#define _FIELD_NAME suNg_scalar_field
#define _FIELD_TYPE suNg_scalar_field
#define _SITE_TYPE suNg_vector
#define _FIELD_DIM 1
#include "TMPL/field_convert.c.tmpl"
#include "TMPL/field_device_transfer.c.tmpl"
#include "TMPL/field_alloc.c.tmpl"

#define _FIELD_NAME avfield
#define _FIELD_TYPE suNg_av_field
#define _SITE_TYPE suNg_algebra_vector
#define _FIELD_DIM 4
#include "TMPL/field_convert.c.tmpl"
#include "TMPL/field_device_transfer.c.tmpl"
#include "TMPL/field_alloc.c.tmpl"

#define _FIELD_NAME gtransf
#define _FIELD_TYPE suNg_field
#define _SITE_TYPE suNg
#define _FIELD_DIM 1
#include "TMPL/field_convert.c.tmpl"
#include "TMPL/field_device_transfer.c.tmpl"
#include "TMPL/field_alloc.c.tmpl"

#define _FIELD_NAME clover_ldl
#define _FIELD_TYPE ldl_field
#define _SITE_TYPE ldl_t
#define _FIELD_DIM 1
#include "TMPL/field_convert.c.tmpl"
#include "TMPL/field_device_transfer.c.tmpl"
#include "TMPL/field_alloc.c.tmpl"

#define _FIELD_NAME clover_term
#define _FIELD_TYPE suNfc_field
#define _SITE_TYPE suNfc
#define _FIELD_DIM 4
#include "TMPL/field_convert.c.tmpl"
#include "TMPL/field_device_transfer.c.tmpl"
#include "TMPL/field_alloc.c.tmpl"

#define _FIELD_NAME clover_force
#define _FIELD_TYPE suNf_field
#define _SITE_TYPE suNf
#define _FIELD_DIM 6
#include "TMPL/field_convert.c.tmpl"
#include "TMPL/field_device_transfer.c.tmpl"
#include "TMPL/field_alloc.c.tmpl"

#define _FIELD_NAME staple_field
#define _FIELD_TYPE suNg_field
#define _SITE_TYPE suNg
#define _FIELD_DIM 3
#include "TMPL/field_convert.c.tmpl"
#include "TMPL/field_device_transfer.c.tmpl"
#include "TMPL/field_alloc.c.tmpl"

#undef _GEOM_TYPE
#undef _IS_SPINOR_LIKE
#undef _FIELD_NAME
#undef _FIELD_TYPE
#undef _SITE_TYPE
#undef _FIELD_DIM
#undef _ALLOCATE
#undef _n

#endif