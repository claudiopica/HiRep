/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *
* All rights reserved.                                                      *
\***************************************************************************/

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

// double precision 
#define _SPINOR_FIELD_TYPE spinor_field
#define _SPINOR_TYPE suNf_spinor
#define _REAL double
#define _COMPLEX hr_complex
#define _SUFFIX _f
#include "TMPL/linear_algebra.h.tmpl"

// single precision 
#define _SPINOR_FIELD_TYPE spinor_field_flt
#define _SPINOR_TYPE suNf_spinor_flt
#define _REAL float
#define _COMPLEX hr_complex_flt
#define _SUFFIX _f_flt
#include "TMPL/linear_algebra.h.tmpl"


#ifdef __cplusplus
  }
#endif
#endif
