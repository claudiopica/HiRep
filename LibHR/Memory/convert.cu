/***************************************************************************\
* Copyright (c) 2023, Sofie Martins                                         *   
* All rights reserved.                                                      * 
\***************************************************************************/

#include "memory.h"
#include "geometry.h"
#include "libhr_core.h"

#ifdef WITH_GPU

#define _FIELD_TYPE spinor_field
#define _SITE_TYPE suNf_spinor
#define _FIELD_DIM 1
#define _IS_SPINOR_LIKE 1
#define _REAL double
#include "TMPL/convert.cu.tmpl"

#define _FIELD_TYPE spinor_field_flt
#define _SITE_TYPE suNf_spinor_flt
#define _FIELD_DIM 1
#define _IS_SPINOR_LIKE 1
#define _REAL float
#include "TMPL/convert.cu.tmpl"

#define _FIELD_TYPE scalar_field
#define _SITE_TYPE double
#define _FIELD_DIM 1
#define _IS_SPINOR_LIKE 1
#define _REAL double
#include "TMPL/convert.cu.tmpl"

#define _FIELD_TYPE suNg_field
#define _SITE_TYPE suNg
#define _FIELD_DIM 4
#define _IS_SPINOR_LIKE 0
#define _REAL double
#include "TMPL/convert.cu.tmpl"

#define _FIELD_TYPE suNg_field_flt
#define _SITE_TYPE suNg_flt
#define _FIELD_DIM 4
#define _IS_SPINOR_LIKE 0
#define _REAL float
#include "TMPL/convert.cu.tmpl"

#define _FIELD_TYPE suNf_field
#define _SITE_TYPE suNf
#define _FIELD_DIM 4
#define _IS_SPINOR_LIKE 0
#define _REAL double
#include "TMPL/convert.cu.tmpl"

#define _FIELD_TYPE suNf_field_flt
#define _SITE_TYPE suNf_flt
#define _FIELD_DIM 4
#define _IS_SPINOR_LIKE 0
#define _REAL float
#include "TMPL/convert.cu.tmpl"

#define _FIELD_TYPE suNg_scalar_field
#define _SITE_TYPE suNg_vector
#define _FIELD_DIM 1
#define _IS_SPINOR_LIKE 0
#define _REAL double
#include "TMPL/convert.cu.tmpl"

#define _FIELD_TYPE suNg_av_field
#define _SITE_TYPE suNg_algebra_vector
#define _FIELD_DIM 4
#define _IS_SPINOR_LIKE 0
#define _REAL double
#include "TMPL/convert.cu.tmpl"

#define _FIELD_TYPE gtransf
#define _SITE_TYPE suNg
#define _FIELD_DIM 1
#define _IS_SPINOR_LIKE 0
#define _REAL double
#include "TMPL/convert.cu.tmpl"

#define _FIELD_TYPE ldl_field
#define _SITE_TYPE ldl_t
#define _FIELD_DIM 1
#define _IS_SPINOR_LIKE 0
#define _REAL double
#include "TMPL/convert.cu.tmpl"

#define _FIELD_TYPE clover_term
#define _SITE_TYPE suNfc
#define _FIELD_DIM 4
#define _IS_SPINOR_LIKE 0
#define _REAL double
#include "TMPL/convert.cu.tmpl"

#define _FIELD_TYPE clover_force
#define _SITE_TYPE suNf
#define _FIELD_DIM 6
#define _IS_SPINOR_LIKE 0
#define _REAL double
#include "TMPL/convert.cu.tmpl"

#define _FIELD_TYPE staple_field
#define _SITE_TYPE suNg
#define _FIELD_DIM 3
#define _IS_SPINOR_LIKE 0
#define _REAL double
#include "TMPL/convert.cu.tmpl"

#endif
