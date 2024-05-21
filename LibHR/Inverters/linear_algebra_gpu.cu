/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *
* All rights reserved.                                                      *
\***************************************************************************/
#ifdef WITH_GPU
//This file should not be compiled if !WITH_GPU

#include "inverters.h"
#include "libhr_core.h"
#include "geometry.h"
#include "Utils/generics.h"

#define _CUDA_FOR(s, ixp, body)                                                            \
    do {                                                                                   \
        _PIECE_FOR((s)->type, (ixp)) {                                                     \
            if (!(s)->type->SAP || (s)->type->parity == PARITY) {                          \
                int N = (s)->type->master_end[(ixp)] - (s)->type->master_start[(ixp)] + 1; \
                unsigned int grid_size = (N - 1) / BLOCK_SIZE_LINEAR_ALGEBRA + 1;          \
                body;                                                                      \
                CudaCheckError();                                                          \
            }                                                                              \
        }                                                                                  \
    } while (0)

#define _GENERIC_DECLARATION(_type, _name, _args, _in_args) \
    _type _name _args {                                     \
        _F_NAME(_name##_, _FIELD_TYPE) _in_args;            \
    }
#define _GENERIC_DECLARATION_RED(_type, _name, _args, _in_args) \
    _type _name _args {                                         \
        return _F_NAME(_name##_, _FIELD_TYPE) _in_args;         \
    }

// Internal macro for defining generic functions and alias function pointers
#define _DECLARE_LINEAR_ALGEBRA_GPU_OP(_type, _name, _args, _in_args) \
    _GENERIC_DECLARATION(_type, _name, _args, _in_args)               \
    _GPU_FUNC(_type, _name##_, _FIELD_TYPE, _args)

#define _DECLARE_LINEAR_ALGEBRA_GPU_RED(_type, _name, _args, _in_args) \
    _GENERIC_DECLARATION_RED(_type, _name, _args, _in_args)            \
    _GPU_FUNC(_type, _name##_, _FIELD_TYPE, _args)

#define _TOCOMPLEX(sp) ((_COMPLEX *)(sp))
#define _NCOM (sizeof(_SITE_TYPE) / sizeof(_COMPLEX))

// Linear Algebra functions are generic
// They are parametrized over the input types for double/single precision
// The template for GPU is in TMPL/linear_algebra_gpu.cu.tmpl

// double precision
#define _FIELD_TYPE spinor_field
#define _SITE_TYPE suNf_spinor
#define _REAL double
#define _COMPLEX hr_complex
#define _FIELD_DIM 1
#define _ISREAL 0
#include "TMPL/linear_algebra_base_operations_gpu.cu.tmpl"
#include "TMPL/linear_algebra_reduction_gpu.cu.tmpl"
#include "TMPL/linear_algebra_lc_gpu.cu.tmpl"
#include "TMPL/linear_algebra_gamma_gpu.cu.tmpl"
#include "TMPL/linear_algebra_base_gpu.cu.tmpl" // This one always needs to go last

// single precision
#define _FIELD_TYPE spinor_field_flt
#define _SITE_TYPE suNf_spinor_flt
#define _REAL float
#define _COMPLEX hr_complex_flt
#define _FIELD_DIM 1
#define _ISREAL 0
#include "TMPL/linear_algebra_base_operations_gpu.cu.tmpl"
#include "TMPL/linear_algebra_reduction_gpu.cu.tmpl"
#include "TMPL/linear_algebra_lc_gpu.cu.tmpl"
#include "TMPL/linear_algebra_gamma_gpu.cu.tmpl"
#include "TMPL/linear_algebra_base_gpu.cu.tmpl"

// double precision
#define _FIELD_TYPE scalar_field
#define _SITE_TYPE double
#define _REAL double
#define _COMPLEX hr_complex
#define _FIELD_DIM 1
#define _ISREAL 1
#include "TMPL/linear_algebra_base_operations_gpu.cu.tmpl"
#include "TMPL/linear_algebra_reduction_gpu.cu.tmpl"
#include "TMPL/linear_algebra_base_gpu.cu.tmpl"

#define _FIELD_TYPE suNg_field
#define _SITE_TYPE suNg
#define _REAL double
#define _COMPLEX hr_complex
#define _FIELD_DIM 4
#define _ISREAL 0
#include "TMPL/linear_algebra_base_operations_gpu.cu.tmpl"
#include "TMPL/linear_algebra_reduction_gpu.cu.tmpl"
#include "TMPL/linear_algebra_base_gpu.cu.tmpl"

#define _FIELD_TYPE suNf_field
#define _SITE_TYPE suNf
#define _REAL double
#define _COMPLEX hr_complex
#define _FIELD_DIM 4
#ifdef REPR_IS_REAL
#define _ISREAL 1
#else
#define _ISREAL 0
#endif
#include "TMPL/linear_algebra_base_operations_gpu.cu.tmpl"
#include "TMPL/linear_algebra_reduction_gpu.cu.tmpl"
#include "TMPL/linear_algebra_base_gpu.cu.tmpl"

#define _FIELD_TYPE suNfc_field
#define _SITE_TYPE suNfc
#define _REAL double
#define _COMPLEX hr_complex
#define _FIELD_DIM 4
#define _ISREAL 0
#include "TMPL/linear_algebra_base_operations_gpu.cu.tmpl"
#include "TMPL/linear_algebra_reduction_gpu.cu.tmpl"
#include "TMPL/linear_algebra_base_gpu.cu.tmpl"

#define _FIELD_TYPE suNg_field_flt
#define _SITE_TYPE suNg_flt
#define _REAL float
#define _COMPLEX hr_complex_flt
#define _FIELD_DIM 4
#define _ISREAL 0
#include "TMPL/linear_algebra_base_operations_gpu.cu.tmpl"
#include "TMPL/linear_algebra_reduction_gpu.cu.tmpl"
#include "TMPL/linear_algebra_base_gpu.cu.tmpl"

#define _FIELD_TYPE suNf_field_flt
#define _SITE_TYPE suNf_flt
#define _REAL float
#define _COMPLEX hr_complex_flt
#define _FIELD_DIM 4
#ifdef REPR_IS_REAL
#define _ISREAL 1
#else
#define _ISREAL 0
#endif
#include "TMPL/linear_algebra_base_operations_gpu.cu.tmpl"
#include "TMPL/linear_algebra_reduction_gpu.cu.tmpl"
#include "TMPL/linear_algebra_base_gpu.cu.tmpl"

#define _FIELD_TYPE suNg_scalar_field
#define _SITE_TYPE suNg_vector
#define _REAL double
#define _COMPLEX hr_complex
#define _FIELD_DIM 1
#ifdef REPR_IS_REAL
#define _ISREAL 1
#else
#define _ISREAL 0
#endif
#include "TMPL/linear_algebra_base_operations_gpu.cu.tmpl"
#include "TMPL/linear_algebra_reduction_gpu.cu.tmpl"
#include "TMPL/linear_algebra_base_gpu.cu.tmpl"

#define _FIELD_TYPE suNg_av_field
#define _SITE_TYPE suNg_algebra_vector
#define _REAL double
#define _COMPLEX hr_complex
#define _FIELD_DIM 4
#define _ISREAL 1
#include "TMPL/linear_algebra_base_operations_gpu.cu.tmpl"
#include "TMPL/linear_algebra_reduction_gpu.cu.tmpl"
#include "TMPL/linear_algebra_base_gpu.cu.tmpl"

#define _FIELD_TYPE gtransf
#define _SITE_TYPE suNg
#define _REAL double
#define _COMPLEX hr_complex
#define _FIELD_DIM 1
#define _ISREAL 0
#include "TMPL/linear_algebra_base_operations_gpu.cu.tmpl"
#include "TMPL/linear_algebra_reduction_gpu.cu.tmpl"
#include "TMPL/linear_algebra_base_gpu.cu.tmpl"

#define _FIELD_TYPE ldl_field
#define _SITE_TYPE ldl_t
#define _REAL double
#define _COMPLEX hr_complex
#define _FIELD_DIM 1
#define _ISREAL 0
#include "TMPL/linear_algebra_base_gpu.cu.tmpl"

#define _FIELD_TYPE clover_term
#define _SITE_TYPE suNfc
#define _REAL double
#define _COMPLEX hr_complex
#define _FIELD_DIM 4
#define _ISREAL 0
#include "TMPL/linear_algebra_base_operations_gpu.cu.tmpl"
#include "TMPL/linear_algebra_reduction_gpu.cu.tmpl"
#include "TMPL/linear_algebra_base_gpu.cu.tmpl"

#define _FIELD_TYPE clover_force
#define _SITE_TYPE suNf
#define _REAL double
#define _COMPLEX hr_complex
#define _FIELD_DIM 6
#ifdef _REPR_IS_REAL
#define _ISREAL 1
#else
#define _ISREAL 0
#endif
#include "TMPL/linear_algebra_base_operations_gpu.cu.tmpl"
#include "TMPL/linear_algebra_reduction_gpu.cu.tmpl"
#include "TMPL/linear_algebra_base_gpu.cu.tmpl"

#define _FIELD_TYPE staple_field
#define _SITE_TYPE suNg
#define _REAL double
#define _COMPLEX hr_complex
#define _FIELD_DIM 3
#define _ISREAL 0
#include "TMPL/linear_algebra_base_operations_gpu.cu.tmpl"
#include "TMPL/linear_algebra_reduction_gpu.cu.tmpl"
#include "TMPL/linear_algebra_base_gpu.cu.tmpl"

#endif