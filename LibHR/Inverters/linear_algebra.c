/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *
* All rights reserved.                                                      *
\***************************************************************************/

#include "inverters.h"
#include "libhr_core.h"
#include "Utils/generics.h"

// Linear Algebra functions are generic
// They are parametrized over the input types for double/single precision

// double precision
#define _FIELD_TYPE spinor_field
#define _SITE_TYPE suNf_spinor
#define _REAL double
#define _COMPLEX hr_complex
#define _FIELD_DIM 1
#include "TMPL/linear_algebra_reduction.c.tmpl"
#include "TMPL/linear_algebra_base_operations.c.tmpl"
#include "TMPL/linear_algebra_lc.c.tmpl"
#include "TMPL/linear_algebra_gamma.c.tmpl"
#include "TMPL/linear_algebra_base.c.tmpl" // This one has to go last

// single precision
#define _FIELD_TYPE spinor_field_flt
#define _SITE_TYPE suNf_spinor_flt
#define _REAL float
#define _COMPLEX hr_complex_flt
#define _FIELD_DIM 1
#include "TMPL/linear_algebra_reduction.c.tmpl"
#include "TMPL/linear_algebra_base_operations.c.tmpl"
#include "TMPL/linear_algebra_lc.c.tmpl"
#include "TMPL/linear_algebra_gamma.c.tmpl"
#include "TMPL/linear_algebra_base.c.tmpl"

// double precision
#define _FIELD_TYPE scalar_field
#define _SITE_TYPE double
#define _REAL double
#define _COMPLEX hr_complex
#define _FIELD_DIM 1
#include "TMPL/linear_algebra_reduction.c.tmpl"
#include "TMPL/linear_algebra_base_operations.c.tmpl"
#include "TMPL/linear_algebra_base.c.tmpl"

#define _FIELD_TYPE suNg_field
#define _SITE_TYPE suNg
#define _REAL double
#define _COMPLEX hr_complex
#define _FIELD_DIM 4
#include "TMPL/linear_algebra_reduction.c.tmpl"
#include "TMPL/linear_algebra_base_operations.c.tmpl"
#include "TMPL/linear_algebra_base.c.tmpl"

#define _FIELD_TYPE suNf_field
#define _SITE_TYPE suNf
#define _REAL double
#define _COMPLEX hr_complex
#define _FIELD_DIM 4
#include "TMPL/linear_algebra_reduction.c.tmpl"
#include "TMPL/linear_algebra_base_operations.c.tmpl"
#include "TMPL/linear_algebra_base.c.tmpl"

#define _FIELD_TYPE suNfc_field
#define _SITE_TYPE suNfc
#define _REAL double
#define _COMPLEX hr_complex
#define _FIELD_DIM 4
#include "TMPL/linear_algebra_reduction.c.tmpl"
#include "TMPL/linear_algebra_base_operations.c.tmpl"
#include "TMPL/linear_algebra_base.c.tmpl"

#define _FIELD_TYPE suNg_field_flt
#define _SITE_TYPE suNg_flt
#define _REAL float
#define _COMPLEX hr_complex_flt
#define _FIELD_DIM 4
#include "TMPL/linear_algebra_reduction.c.tmpl"
#include "TMPL/linear_algebra_base_operations.c.tmpl"
#include "TMPL/linear_algebra_base.c.tmpl"

#define _FIELD_TYPE suNf_field_flt
#define _SITE_TYPE suNf_flt
#define _REAL float
#define _COMPLEX hr_complex_flt
#define _FIELD_DIM 4
#include "TMPL/linear_algebra_reduction.c.tmpl"
#include "TMPL/linear_algebra_base_operations.c.tmpl"
#include "TMPL/linear_algebra_base.c.tmpl"

#define _FIELD_TYPE suNg_scalar_field
#define _SITE_TYPE suNg_vector
#define _REAL double
#define _COMPLEX hr_complex
#define _FIELD_DIM 1
#include "TMPL/linear_algebra_reduction.c.tmpl"
#include "TMPL/linear_algebra_base_operations.c.tmpl"
#include "TMPL/linear_algebra_base.c.tmpl"

#define _FIELD_TYPE suNg_av_field
#define _SITE_TYPE suNg_algebra_vector
#define _REAL double
#define _COMPLEX hr_complex
#define _FIELD_DIM 4
#include "TMPL/linear_algebra_reduction.c.tmpl"
#include "TMPL/linear_algebra_base_operations.c.tmpl"
#include "TMPL/linear_algebra_base.c.tmpl"

#define _FIELD_TYPE gtransf
#define _SITE_TYPE suNg
#define _REAL double
#define _COMPLEX hr_complex
#define _FIELD_DIM 1
#include "TMPL/linear_algebra_reduction.c.tmpl"
#include "TMPL/linear_algebra_base_operations.c.tmpl"
#include "TMPL/linear_algebra_base.c.tmpl"

#define _FIELD_TYPE ldl_field
#define _SITE_TYPE ldl_t
#define _REAL double
#define _COMPLEX hr_complex
#define _FIELD_DIM 1
#include "TMPL/linear_algebra_base.c.tmpl"

#define _FIELD_TYPE clover_term
#define _SITE_TYPE suNfc
#define _REAL double
#define _COMPLEX hr_complex
#define _FIELD_DIM 4
#include "TMPL/linear_algebra_reduction.c.tmpl"
#include "TMPL/linear_algebra_base_operations.c.tmpl"
#include "TMPL/linear_algebra_base.c.tmpl"

#define _FIELD_TYPE clover_force
#define _SITE_TYPE suNf
#define _REAL double
#define _COMPLEX hr_complex
#define _FIELD_DIM 6
#include "TMPL/linear_algebra_reduction.c.tmpl"
#include "TMPL/linear_algebra_base_operations.c.tmpl"
#include "TMPL/linear_algebra_base.c.tmpl"

#define _FIELD_TYPE staple_field
#define _SITE_TYPE suNg
#define _REAL double
#define _COMPLEX hr_complex
#define _FIELD_DIM 3
#include "TMPL/linear_algebra_reduction.c.tmpl"
#include "TMPL/linear_algebra_base_operations.c.tmpl"
#include "TMPL/linear_algebra_base.c.tmpl"
