/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *
* All rights reserved.                                                      *
\***************************************************************************/
#ifdef WITH_GPU
//This file should not be compiled if !WITH_GPU

#include "global.h"
#include <string.h>
#include "spinor_field.h"
#include "gamma_spinor.h"
#include "hr_complex.h"
// #include "communications.h"
#include "suN.h"
#include "suN_types.h"
#include "error.h"
#include "linear_algebra.h"

#include "./linear_algebra_gpu_kernels.hpp"

/*
 * LINEAR ALGEBRA FUNCTIONS ARE DEFINED IN THE TEMPLATE
 *
 * TMPL/linear_algebra.c.sdtmpl
 *
 */

/* double precision */
  // Declare types for double precision template parsing
#define _SPINOR_FIELD_TYPE spinor_field
#define _SPINOR_TYPE suNf_spinor
#define _REAL double
#define _COMPLEX hr_complex

// Use GPU template to declare double precision functions (C++)
#define _FUNC(a,b,c) a b##_f_gpu c
#define _BODY(a) a
#include "TMPL/linear_algebra_gpu.cu.sdtmpl"
#undef _FUNC
#undef _BODY

#undef _FUNC
#undef _BODY
#undef _SPINOR_FIELD_TYPE
#undef _SPINOR_TYPE
#undef _REAL
#undef _COMPLEX

/* single precision */
//Declare types for single precision
#define _SPINOR_FIELD_TYPE spinor_field_flt
#define _SPINOR_TYPE suNf_spinor_flt
#define _REAL float
#define _COMPLEX hr_complex_flt

// C++ GPU function template
#define _FUNC(a,b,c) a b##_f_flt_gpu c
#define _BODY(a) a
#include "TMPL/linear_algebra_gpu.cu.sdtmpl"
#undef _FUNC
#undef _BODY

  //Undefine single precision definitions.
#undef _FUNC
#undef _BODY
#undef _SPINOR_FIELD_TYPE
#undef _SPINOR_TYPE
#undef _REAL
#undef _COMPLEX

#endif
