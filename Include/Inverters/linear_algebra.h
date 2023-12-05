/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *
* All rights reserved.                                                      *
\***************************************************************************/

/// Headerfile for:
/// - linear_algebra.c
/// - linear_algebra_gpu.cu

/**
 * @file linear_algebra.h
 * @brief Linear algebra operations on spinors both for CPU and with GPU
 */

#ifndef LINEAR_ALGEBRA_H
#define LINEAR_ALGEBRA_H

#include "libhr_core.h"
#include "hr_complex.h"
#include "Utils/generics.h"

#ifdef __cplusplus
extern "C" {
#endif

// Linear Algebra functions are generic
// They are parametrized over the input types for double/single precision
// The template is in TMPL/linear_algebra.h.tmpl

#ifdef WITH_GPU
double *alloc_double_sum_field(int n);
hr_complex *alloc_complex_sum_field(int n);
#endif

// double precision
#define _FIELD_TYPE spinor_field
#define _REAL double
#define _COMPLEX hr_complex
#include "TMPL/linear_algebra_def.h.tmpl" // This has to go first
#include "TMPL/linear_algebra_reduction.h.tmpl"
#include "TMPL/linear_algebra_base_operations.h.tmpl"
#include "TMPL/linear_algebra_lc.h.tmpl"
#include "TMPL/linear_algebra_gamma.h.tmpl"
#include "TMPL/linear_algebra_base.h.tmpl" // This has to come last

// single precision
#define _FIELD_TYPE spinor_field_flt
#define _REAL float
#define _COMPLEX hr_complex_flt
#include "TMPL/linear_algebra_def.h.tmpl"
#include "TMPL/linear_algebra_reduction.h.tmpl"
#include "TMPL/linear_algebra_base_operations.h.tmpl"
#include "TMPL/linear_algebra_lc.h.tmpl"
#include "TMPL/linear_algebra_gamma.h.tmpl"
#include "TMPL/linear_algebra_base.h.tmpl"

// double precision
#define _FIELD_TYPE scalar_field
#define _REAL double
#define _COMPLEX hr_complex
#include "TMPL/linear_algebra_def.h.tmpl"
#include "TMPL/linear_algebra_reduction.h.tmpl"
#include "TMPL/linear_algebra_base_operations.h.tmpl"
#include "TMPL/linear_algebra_base.h.tmpl"

#define _FIELD_TYPE suNg_field
#define _REAL double
#define _COMPLEX hr_complex
#include "TMPL/linear_algebra_def.h.tmpl"
#include "TMPL/linear_algebra_reduction.h.tmpl"
#include "TMPL/linear_algebra_base_operations.h.tmpl"
#include "TMPL/linear_algebra_base.h.tmpl"

#define _FIELD_TYPE suNf_field
#define _REAL double
#define _COMPLEX hr_complex
#include "TMPL/linear_algebra_def.h.tmpl"
#include "TMPL/linear_algebra_reduction.h.tmpl"
#include "TMPL/linear_algebra_base_operations.h.tmpl"
#include "TMPL/linear_algebra_base.h.tmpl"

#define _FIELD_TYPE suNfc_field
#define _REAL double
#define _COMPLEX hr_complex
#include "TMPL/linear_algebra_def.h.tmpl"
#include "TMPL/linear_algebra_reduction.h.tmpl"
#include "TMPL/linear_algebra_base_operations.h.tmpl"
#include "TMPL/linear_algebra_base.h.tmpl"

#define _FIELD_TYPE suNg_field_flt
#define _REAL float
#define _COMPLEX hr_complex_flt
#include "TMPL/linear_algebra_def.h.tmpl"
#include "TMPL/linear_algebra_reduction.h.tmpl"
#include "TMPL/linear_algebra_base_operations.h.tmpl"
#include "TMPL/linear_algebra_base.h.tmpl"

#define _FIELD_TYPE suNf_field_flt
#define _REAL float
#define _COMPLEX hr_complex_flt
#include "TMPL/linear_algebra_def.h.tmpl"
#include "TMPL/linear_algebra_reduction.h.tmpl"
#include "TMPL/linear_algebra_base_operations.h.tmpl"
#include "TMPL/linear_algebra_base.h.tmpl"

#define _FIELD_TYPE suNg_scalar_field
#define _REAL double
#define _COMPLEX hr_complex
#include "TMPL/linear_algebra_def.h.tmpl"
#include "TMPL/linear_algebra_reduction.h.tmpl"
#include "TMPL/linear_algebra_base_operations.h.tmpl"
#include "TMPL/linear_algebra_base.h.tmpl"

#define _FIELD_TYPE suNg_av_field
#define _REAL double
#define _COMPLEX hr_complex
#include "TMPL/linear_algebra_def.h.tmpl"
#include "TMPL/linear_algebra_reduction.h.tmpl"
#include "TMPL/linear_algebra_base_operations.h.tmpl"
#include "TMPL/linear_algebra_base.h.tmpl"

#define _FIELD_TYPE gtransf
#define _REAL double
#define _COMPLEX hr_complex
#include "TMPL/linear_algebra_def.h.tmpl"
#include "TMPL/linear_algebra_reduction.h.tmpl"
#include "TMPL/linear_algebra_base_operations.h.tmpl"
#include "TMPL/linear_algebra_base.h.tmpl"

#define _FIELD_TYPE ldl_field
#define _REAL double
#define _COMPLEX hr_complex
#include "TMPL/linear_algebra_def.h.tmpl"
#include "TMPL/linear_algebra_base.h.tmpl"

#define _FIELD_TYPE clover_term
#define _REAL double
#define _COMPLEX hr_complex
#include "TMPL/linear_algebra_def.h.tmpl"
#include "TMPL/linear_algebra_reduction.h.tmpl"
#include "TMPL/linear_algebra_base_operations.h.tmpl"
#include "TMPL/linear_algebra_base.h.tmpl"

#define _FIELD_TYPE clover_force
#define _REAL double
#define _COMPLEX hr_complex
#include "TMPL/linear_algebra_def.h.tmpl"
#include "TMPL/linear_algebra_reduction.h.tmpl"
#include "TMPL/linear_algebra_base_operations.h.tmpl"
#include "TMPL/linear_algebra_base.h.tmpl"

#define _FIELD_TYPE staple_field
#define _REAL double
#define _COMPLEX hr_complex
#include "TMPL/linear_algebra_def.h.tmpl"
#include "TMPL/linear_algebra_reduction.h.tmpl"
#include "TMPL/linear_algebra_base_operations.h.tmpl"
#include "TMPL/linear_algebra_base.h.tmpl"

#undef _DECLARE_LINA_HEADER

#ifdef __cplusplus
}
#endif
#endif
