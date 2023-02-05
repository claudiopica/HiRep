/***************************************************************************\
* Copyright (c) 2023, Sofie Martins                                         *
* All rights reserved.                                                      *
\***************************************************************************/
#ifdef WITH_GPU

#include "libhr_core.h"
#include "inverters.h"
#include "geometry.h"
#include "update.h"
#include "memory.h"

#ifdef ROTATED_SF
    extern rhmc_par _update_par;
#endif

// double precision
#define _FIELD_TYPE spinor_field
#define _FIELD_NAME spinor_field_f
#define _SPINOR_TYPE suNf_spinor
#define _HSPINOR_TYPE suNf_hspinor
#define _VECTOR_TYPE suNf_vector
#define _COMPLEX hr_complex
#define _REAL double
#define _GAUGE_TYPE suNf
#define _SUFFIX 
#define _REP_SUFFIX _f
#include "TMPL/Dphi_gpu.cu.tmpl"

// single precision
#define _FIELD_TYPE spinor_field_flt
#define _FIELD_NAME spinor_field_f_flt
#define _SPINOR_TYPE suNf_spinor_flt
#define _HSPINOR_TYPE suNf_hspinor_flt
#define _VECTOR_TYPE suNf_vector_flt
#define _COMPLEX hr_complex_flt
#define _REAL float
#define _GAUGE_TYPE suNf_flt
#define _SUFFIX _flt
#define _REP_SUFFIX _f_flt
#include "TMPL/Dphi_gpu.cu.tmpl"

#endif 