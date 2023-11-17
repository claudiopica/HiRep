/***************************************************************************\
* Copyright (c) 2022, Sofie Martins                                         *   
* All rights reserved.                                                      * 
\***************************************************************************/

#include "geometry.h"
#include "libhr_core.h"
#include "update.h"

#define _FIELD_TYPE spinor_field
#define _SITE_TYPE suNf_spinor
#define _HSPINOR_TYPE suNf_hspinor
#define _GAUGE_TYPE suNf
#define _FIELD_DIM 1
#define _GEOM_TYPE spinor
#define _COMPLEX hr_complex
#define _REAL double
#define _REP_SUFFIX _f
#include "TMPL/sync_to_buffer_reduced.cu.tmpl"

#define _FIELD_TYPE spinor_field_flt
#define _SITE_TYPE suNf_spinor_flt
#define _HSPINOR_TYPE suNf_hspinor_flt
#define _GAUGE_TYPE suNf_flt
#define _FIELD_DIM 1
#define _GEOM_TYPE spinor
#define _COMPLEX hr_complex_flt
#define _REAL float
#define _REP_SUFFIX _f_flt
#include "TMPL/sync_to_buffer_reduced.cu.tmpl"