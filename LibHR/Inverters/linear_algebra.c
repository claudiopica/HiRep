/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *
* All rights reserved.                                                      *
\***************************************************************************/

#include "inverters.h"
#include "libhr_core.h"

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

//set double precision function aliases
#ifdef WITH_GPU
  #define _ALIAS_FUNC(a,b,c) a (*b##_f) c = b##_f_gpu
#else
  #define _ALIAS_FUNC(a,b,c) a (*b##_f) c = b##_f_cpu
#endif

// Use CPU templates for double precision functions & aliasing (C)
#define _FUNC(a,b,c) _ALIAS_FUNC(a,b,c); a b##_f_cpu c
#include "TMPL/linear_algebra.c.sdtmpl"
#undef _FUNC
#undef _ALIAS_FUNC

 //Undefine all definitions to move on to single precision#undef _SPINOR_FIELD_TYPE
#undef _SPINOR_FIELD_TYPE
#undef _SPINOR_TYPE
#undef _REAL
#undef _COMPLEX

#define _SPINOR_FIELD_TYPE spinor_field_flt
#define _SPINOR_TYPE suNf_spinor_flt
#define _REAL float
#define _COMPLEX hr_complex_flt

//set single precision function aliases
#ifdef WITH_GPU
  #define _ALIAS_FUNC(a,b,c) a (*b##_f_flt) c = b##_f_flt_gpu
#else
  #define _ALIAS_FUNC(a,b,c) a (*b##_f_flt) c = b##_f_flt_cpu
#endif

#define _FUNC(a,b,c) _ALIAS_FUNC(a,b,c); a b##_f_flt_cpu c
#include "TMPL/linear_algebra.c.sdtmpl"
#undef _FUNC
#undef _ALIAS_FUNC

#undef _SPINOR_FIELD_TYPE
#undef _SPINOR_TYPE
#undef _REAL
#undef _COMPLEX
