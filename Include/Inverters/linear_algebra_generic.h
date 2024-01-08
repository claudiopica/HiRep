/***************************************************************************\
* Copyright (c) 2023, Sofie Martins                                         *
* All rights reserved.                                                      *
\***************************************************************************/

/**
 * @file linear_algebra_generic.h
 * @brief Linear algebra operations on spinors both for CPU and with GPU
 */

#ifndef LINEAR_ALGEBRA_GENERIC_H
#define LINEAR_ALGEBRA_GENERIC_H

#include "libhr_core.h"
#include <stdio.h>
#define template_error error(1, 1, __func__, "Complex multiplication of real-valued field \n");

#ifdef REPR_IS_REAL
#define __mul_add_assign_suNfc(s1, rho, s2) suNfc * : ({ _suNfc_mul_add(*(suNfc *)s1, 1.0, *(suNfc *)s1, rho, *(suNfc *)s2); }),
#else
#define __mul_add_assign_suNfc(s1, rho, s2)
#endif

#define mul_add_assign(s1, rho, s2)                                                                                        \
    _Generic((s2),                                                                                                         \
        spinor_field *: mul_add_assign_spinor_field((spinor_field *)s1, rho, (spinor_field *)s2),                          \
        spinor_field_flt *: mul_add_assign_spinor_field_flt((spinor_field_flt *)s1, rho, (spinor_field_flt *)s2),          \
        scalar_field *: mul_add_assign_scalar_field((scalar_field *)s1, rho, (scalar_field *)s2),                          \
        suNg_field *: mul_add_assign_suNg_field((suNg_field *)s1, rho, (suNg_field *)s2),                                  \
        suNf_field *: mul_add_assign_suNf_field((suNf_field *)s1, rho, (suNf_field *)s2),                                  \
        suNfc_field *: mul_add_assign_suNfc_field((suNfc_field *)s1, rho, (suNfc_field *)s2),                              \
        suNg_field_flt *: mul_add_assign_suNg_field_flt((suNg_field_flt *)s1, rho, (suNg_field_flt *)s2),                  \
        suNf_field_flt *: mul_add_assign_suNf_field_flt((suNf_field_flt *)s1, rho, (suNf_field_flt *)s2),                  \
        suNg_scalar_field *: mul_add_assign_suNg_scalar_field((suNg_scalar_field *)s1, rho, (suNg_scalar_field *)s2),      \
        suNg_av_field *: mul_add_assign_suNg_av_field((suNg_av_field *)s1, rho, (suNg_av_field *)s2),                      \
        gtransf *: mul_add_assign_gtransf((gtransf *)s1, rho, (gtransf *)s2),                                              \
        clover_term *: mul_add_assign_clover_term((clover_term *)s1, rho, (clover_term *)s2),                              \
        clover_force *: mul_add_assign_clover_force((clover_force *)s1, rho, (clover_force *)s2),                          \
        staple_field *: mul_add_assign_staple_field((staple_field *)s1, rho, (staple_field *)s2),                          \
        suNf_spinor *: ({ _spinor_mul_add_assign_f(*(suNf_spinor *)s1, rho, *(suNf_spinor *)s2); }),                       \
        suNf_spinor_flt *: ({ _spinor_mul_add_assign_f(*(suNf_spinor_flt *)s1, rho, *(suNf_spinor_flt *)s2); }),           \
        __mul_add_assign_suNfc(s1, rho, s2) suNf *: ({ _suNf_mul_add(*(suNf *)s1, 1.0, *(suNf *)s1, rho, *(suNf *)s2); }), \
        suNg *: ({ _suNg_mul_add(*(suNg *)s1, 1.0, *(suNg *)s1, rho, *(suNg *)s2); }),                                     \
        suNf_flt *: ({ _suNf_mul_add(*(suNf_flt *)s1, 1.0, *(suNf_flt *)s1, rho, *(suNf_flt *)s2); }),                     \
        suNg_flt *: ({ _suNg_mul_add(*(suNg_flt *)s1, 1.0, *(suNg_flt *)s1, rho, *(suNg_flt *)s2); }),                     \
        suNf_vector *: ({ _vector_mul_add_assign_f(*(suNf_vector *)s1, rho, *(suNf_vector *)s2); }),                       \
        suNg_vector *: ({ _vector_mul_add_assign_g(*(suNg_vector *)s1, rho, *(suNg_vector *)s2); }),                       \
        suNg_algebra_vector *: ({                                                                                          \
            _algebra_vector_mul_add_assign_g(*(suNg_algebra_vector *)s1, rho, *(suNg_algebra_vector *)s2);                 \
        }),                                                                                                                \
        double *: ({ (*(double *)s1) += (rho) * (*(double *)s2); }),                                                       \
        float *: ({ (*(double *)s1) += (rho) * (*(double *)s2); }))

#ifdef REPR_IS_REAL
#define __mulc_add_assign_suNfc(s1, rho, s2)    \
    suNfc * : ({                                \
        suNfc __tmp;                            \
        _suNfc_mulc(__tmp, rho, *(suNfc *)s2);  \
        _suNfc_add_assign(*(suNfc *)s1, __tmp); \
    }),

#define __mulc_add_assign_suNf(s1, rho, s2) suNf * : ({ template_error; }),
#define __mulc_add_assign_suNf_flt(s1, rho, s2) suNf_flt * : ({ template_error; }),
#else
#define __mulc_add_assign_suNfc(s1, rho, s2)
#define __mulc_add_assign_suNf(s1, rho, s2)   \
    suNf * : ({                               \
        suNf __tmp;                           \
        _suNf_mulc(__tmp, rho, *(suNf *)s2);  \
        _suNf_add_assign(*(suNf *)s1, __tmp); \
    }),

#define __mulc_add_assign_suNf_flt(s1, rho, s2)   \
    suNf_flt * : ({                               \
        suNf_flt __tmp;                           \
        _suNf_mulc(__tmp, rho, *(suNf_flt *)s2);  \
        _suNf_add_assign(*(suNf_flt *)s1, __tmp); \
    }),
#endif

#define mulc_add_assign(s1, rho, s2)                                                                                     \
    _Generic((s2),                                                                                                       \
        spinor_field *: mulc_add_assign_spinor_field((spinor_field *)s1, rho, (spinor_field *)s2),                       \
        spinor_field_flt *: mulc_add_assign_spinor_field_flt((spinor_field_flt *)s1, rho, (spinor_field_flt *)s2),       \
        scalar_field *: mulc_add_assign_scalar_field((scalar_field *)s1, rho, (scalar_field *)s2),                       \
        suNg_field *: mulc_add_assign_suNg_field((suNg_field *)s1, rho, (suNg_field *)s2),                               \
        suNf_field *: mulc_add_assign_suNf_field((suNf_field *)s1, rho, (suNf_field *)s2),                               \
        suNfc_field *: mulc_add_assign_suNfc_field((suNfc_field *)s1, rho, (suNfc_field *)s2),                           \
        suNg_field_flt *: mulc_add_assign_suNg_field_flt((suNg_field_flt *)s1, rho, (suNg_field_flt *)s2),               \
        suNf_field_flt *: mulc_add_assign_suNf_field_flt((suNf_field_flt *)s1, rho, (suNf_field_flt *)s2),               \
        suNg_scalar_field *: mulc_add_assign_suNg_scalar_field((suNg_scalar_field *)s1, rho, (suNg_scalar_field *)s2),   \
        suNg_av_field *: mulc_add_assign_suNg_av_field((suNg_av_field *)s1, rho, (suNg_av_field *)s2),                   \
        gtransf *: mulc_add_assign_gtransf((gtransf *)s1, rho, (gtransf *)s2),                                           \
        clover_term *: mulc_add_assign_clover_term((clover_term *)s1, rho, (clover_term *)s2),                           \
        clover_force *: mulc_add_assign_clover_force((clover_force *)s1, rho, (clover_force *)s2),                       \
        staple_field *: mulc_add_assign_staple_field((staple_field *)s1, rho, (staple_field *)s2),                       \
        suNf_spinor *: ({ _spinor_mulc_add_assign_f(*(suNf_spinor *)s1, (hr_complex)rho, *(suNf_spinor *)s2); }),        \
        suNf_spinor_flt *: ({                                                                                            \
            _spinor_mulc_add_assign_f(*(suNf_spinor_flt *)s1, (hr_complex_flt)rho, *(suNf_spinor_flt *)s2);              \
        }),                                                                                                              \
        __mulc_add_assign_suNf(s1, rho, s2) __mulc_add_assign_suNf_flt(s1, rho, s2) __mulc_add_assign_suNfc(s1, rho, s2) \
                suNg *: ({                                                                                               \
                    suNg __tmp;                                                                                          \
                    _suNg_mulc(__tmp, rho, *(suNg *)s2);                                                                 \
                    _suNg_add_assign(*(suNg *)s1, __tmp);                                                                \
                }),                                                                                                      \
                                                                                                                         \
        suNg_flt *: ({                                                                                                   \
            suNg_flt __tmp;                                                                                              \
            _suNg_mulc(__tmp, rho, *(suNg_flt *)s2);                                                                     \
            _suNg_add_assign(*(suNg_flt *)s1, __tmp);                                                                    \
        }),                                                                                                              \
        suNf_vector *: ({ _vector_mulc_add_assign_f(*(suNf_vector *)s1, rho, *(suNf_vector *)s2); }),                    \
        suNg_vector *: ({ _vector_mulc_add_assign_g(*(suNg_vector *)s1, rho, *(suNg_vector *)s2); }),                    \
        suNg_algebra_vector *: ({ template_error; }),                                                                    \
        double *: ({ template_error; }),                                                                                 \
        float *: ({ template_error; }))

#ifdef REPR_IS_REAL
#define __mul_suNfc(s1, rho, s2) suNfc * : ({ _suNfc_mul((*((suNfc *)s1)), rho, (*((suNfc *)s2))); }),
#else
#define __mul_suNfc(s1, rho, s2)
#endif

#define mul(s1, rho, s2)                                                                                                  \
    _Generic((s2),                                                                                                        \
        spinor_field *: mul_spinor_field((spinor_field *)s1, rho, (spinor_field *)s2),                                    \
        spinor_field_flt *: mul_spinor_field_flt((spinor_field_flt *)s1, rho, (spinor_field_flt *)s2),                    \
        scalar_field *: mul_scalar_field((scalar_field *)s1, rho, (scalar_field *)s2),                                    \
        suNg_field *: mul_suNg_field((suNg_field *)s1, rho, (suNg_field *)s2),                                            \
        suNf_field *: mul_suNf_field((suNf_field *)s1, rho, (suNf_field *)s2),                                            \
        suNfc_field *: mul_suNfc_field((suNfc_field *)s1, rho, (suNfc_field *)s2),                                        \
        suNg_field_flt *: mul_suNg_field_flt((suNg_field_flt *)s1, rho, (suNg_field_flt *)s2),                            \
        suNf_field_flt *: mul_suNf_field_flt((suNf_field_flt *)s1, rho, (suNf_field_flt *)s2),                            \
        suNg_scalar_field *: mul_suNg_scalar_field((suNg_scalar_field *)s1, rho, (suNg_scalar_field *)s2),                \
        suNg_av_field *: mul_suNg_av_field((suNg_av_field *)s1, rho, (suNg_av_field *)s2),                                \
        gtransf *: mul_gtransf((gtransf *)s1, rho, (gtransf *)s2),                                                        \
        clover_term *: mul_clover_term((clover_term *)s1, rho, (clover_term *)s2),                                        \
        clover_force *: mul_clover_force((clover_force *)s1, rho, (clover_force *)s2),                                    \
        staple_field *: mul_staple_field((staple_field *)s1, rho, (staple_field *)s2),                                    \
        suNf_spinor *: ({ _spinor_mul_f(*(suNf_spinor *)s1, (hr_complex)rho, *(suNf_spinor *)s2); }),                     \
        suNf_spinor_flt *: ({ _spinor_mul_f(*(suNf_spinor_flt *)s1, (hr_complex_flt)rho, *(suNf_spinor_flt *)s2); }),     \
        __mul_suNfc(s1, rho, s2) suNf *: ({ _suNf_mul((*((suNf *)s1)), rho, (*((suNf *)s2))); }),                         \
        suNg *: ({ _suNg_mul(*(suNg *)s1, rho, *(suNg *)s2); }),                                                          \
        suNf_flt *: ({ _suNf_mul(*(suNf_flt *)s1, rho, *(suNf_flt *)s2); }),                                              \
        suNg_flt *: ({ _suNg_mul(*(suNg_flt *)s1, rho, *(suNg_flt *)s2); }),                                              \
        suNf_vector *: ({ _vector_mul_f(*(suNf_vector *)s1, rho, *(suNf_vector *)s2); }),                                 \
        suNg_vector *: ({ _vector_mul_g(*(suNg_vector *)s1, rho, *(suNg_vector *)s2); }),                                 \
        suNg_algebra_vector *: ({ _algebra_vector_mul_g(*(suNg_algebra_vector *)s1, rho, *(suNg_algebra_vector *)s2); }), \
        double *: ({ (*(double *)s1) = (rho) * (*(double *)s2); }),                                                       \
        float *: ({ (*(double *)s1) = (rho) * (*(double *)s2); }))

#ifdef REPR_IS_REAL
#define __mulc_suNfc(s1, rho, s2) suNfc * : ({ _suNfc_mulc(*(suNfc *)s1, rho, *(suNfc *)s2); }),
#define __mulc_suNf(s1, rho, s2) suNf * : ({ template_error; }),
#define __mulc_suNf_flt(s1, rho, s2) suNf_flt * : ({ template_error; }),
#else
#define __mulc_suNfc(s1, rho, s2)
#define __mulc_suNf(s1, rho, s2) suNf * : ({ _suNf_mulc((*((suNf *)s1)), rho, (*((suNf *)s2))); }),
#define __mulc_suNf_flt(s1, rho, s2) suNf_flt * : ({ _suNf_mulc(*(suNf_flt *)s1, rho, *(suNf_flt *)s2); }),
#endif

#define mulc(s1, rho, s2)                                                                                              \
    _Generic((s2),                                                                                                     \
        spinor_field *: mulc_spinor_field((spinor_field *)s1, rho, (spinor_field *)s2),                                \
        spinor_field_flt *: mulc_spinor_field_flt((spinor_field_flt *)s1, rho, (spinor_field_flt *)s2),                \
        scalar_field *: mulc_scalar_field((scalar_field *)s1, rho, (scalar_field *)s2),                                \
        suNg_field *: mulc_suNg_field((suNg_field *)s1, rho, (suNg_field *)s2),                                        \
        suNf_field *: mulc_suNf_field((suNf_field *)s1, rho, (suNf_field *)s2),                                        \
        suNfc_field *: mulc_suNfc_field((suNfc_field *)s1, rho, (suNfc_field *)s2),                                    \
        suNg_field_flt *: mulc_suNg_field_flt((suNg_field_flt *)s1, rho, (suNg_field_flt *)s2),                        \
        suNf_field_flt *: mulc_suNf_field_flt((suNf_field_flt *)s1, rho, (suNf_field_flt *)s2),                        \
        suNg_scalar_field *: mulc_suNg_scalar_field((suNg_scalar_field *)s1, rho, (suNg_scalar_field *)s2),            \
        suNg_av_field *: mulc_suNg_av_field((suNg_av_field *)s1, rho, (suNg_av_field *)s2),                            \
        gtransf *: mulc_gtransf((gtransf *)s1, rho, (gtransf *)s2),                                                    \
        clover_term *: mulc_clover_term((clover_term *)s1, rho, (clover_term *)s2),                                    \
        clover_force *: mulc_clover_force((clover_force *)s1, rho, (clover_force *)s2),                                \
        staple_field *: mulc_staple_field((staple_field *)s1, rho, (staple_field *)s2),                                \
        suNf_spinor *: ({ _spinor_mulc_f(*(suNf_spinor *)s1, (hr_complex)rho, *(suNf_spinor *)s2); }),                 \
        suNf_spinor_flt *: ({ _spinor_mulc_f(*(suNf_spinor_flt *)s1, (hr_complex_flt)rho, *(suNf_spinor_flt *)s2); }), \
        __mulc_suNfc(s1, rho, s2) __mulc_suNf(s1, rho, s2) __mulc_suNf_flt(s1, rho, s2)                                \
                suNg *: ({ _suNg_mulc(*(suNg *)s1, rho, *(suNg *)s2); }),                                              \
        suNg_flt *: ({ _suNg_mulc(*(suNg_flt *)s1, rho, *(suNg_flt *)s2); }),                                          \
        suNf_vector *: ({ _vector_mulc_f(*(suNf_vector *)s1, rho, *(suNf_vector *)s2); }),                             \
        suNg_vector *: ({ _vector_mulc_g(*(suNg_vector *)s1, rho, *(suNg_vector *)s2); }),                             \
        suNg_algebra_vector *: ({                                                                                      \
            _algebra_vector_mul_g(*(suNg_algebra_vector *)s1, creal(rho), *(suNg_algebra_vector *)s2);                 \
        }),                                                                                                            \
        double *: ({ template_error; }),                                                                               \
        float *: ({ template_error; }))

#ifdef REPR_IS_REAL
#define __add_suNfc(r, s1, s2) suNfc * : ({ _suNfc_mul_add(*(suNfc *)r, 1.0, *(suNfc *)s1, 1.0, *(suNfc *)s2); }),
#else
#define __add_suNfc(r, s1, s2)
#endif

#define add(r, s1, s2)                                                                                                        \
    _Generic((s2),                                                                                                            \
        spinor_field *: add_spinor_field((spinor_field *)r, (spinor_field *)s1, (spinor_field *)s2),                          \
        spinor_field_flt *: add_spinor_field_flt((spinor_field_flt *)r, (spinor_field_flt *)s1, (spinor_field_flt *)s2),      \
        scalar_field *: add_scalar_field((scalar_field *)r, (scalar_field *)s1, (scalar_field *)s2),                          \
        suNg_field *: add_suNg_field((suNg_field *)r, (suNg_field *)s1, (suNg_field *)s2),                                    \
        suNf_field *: add_suNf_field((suNf_field *)r, (suNf_field *)s1, (suNf_field *)s2),                                    \
        suNfc_field *: add_suNfc_field((suNfc_field *)r, (suNfc_field *)s1, (suNfc_field *)s2),                               \
        suNg_field_flt *: add_suNg_field_flt((suNg_field_flt *)r, (suNg_field_flt *)s1, (suNg_field_flt *)s2),                \
        suNf_field_flt *: add_suNf_field_flt((suNf_field_flt *)r, (suNf_field_flt *)s1, (suNf_field_flt *)s2),                \
        suNg_scalar_field *: add_suNg_scalar_field((suNg_scalar_field *)r, (suNg_scalar_field *)s1, (suNg_scalar_field *)s2), \
        suNg_av_field *: add_suNg_av_field((suNg_av_field *)r, (suNg_av_field *)s1, (suNg_av_field *)s2),                     \
        gtransf *: add_gtransf((gtransf *)r, (gtransf *)s1, (gtransf *)s2),                                                   \
        clover_term *: add_clover_term((clover_term *)r, (clover_term *)s1, (clover_term *)s2),                               \
        clover_force *: add_clover_force((clover_force *)r, (clover_force *)s1, (clover_force *)s2),                          \
        staple_field *: add_staple_field((staple_field *)r, (staple_field *)s1, (staple_field *)s2),                          \
        suNf_spinor *: ({ _spinor_add_f(*(suNf_spinor *)r, *(suNf_spinor *)s1, *(suNf_spinor *)s2); }),                       \
        suNf_spinor_flt *: ({ _spinor_add_f(*(suNf_spinor_flt *)r, *(suNf_spinor_flt *)s1, *(suNf_spinor_flt *)s2); }),       \
        __add_suNfc(r, s1, s2) suNf *: ({ _suNf_mul_add(*(suNf *)r, 1.0, *(suNf *)s1, 1.0, *(suNf *)s2); }),                  \
        suNg *: ({ _suNg_mul_add(*(suNg *)r, 1.0, *(suNg *)s1, 1.0, *(suNg *)s2); }),                                         \
        suNf_flt *: ({ _suNf_mul_add(*(suNf_flt *)r, 1.0, *(suNf_flt *)s1, 1.0, *(suNf_flt *)s2); }),                         \
        suNg_flt *: ({ _suNg_mul_add(*(suNg_flt *)r, 1.0, *(suNg_flt *)s1, 1.0, *(suNg_flt *)s2); }),                         \
        suNf_vector *: ({ _vector_add_f(*(suNf_vector *)r, *(suNf_vector *)s1, *(suNf_vector *)s2); }),                       \
        suNg_vector *: ({ _vector_add_g(*(suNg_vector *)r, *(suNg_vector *)s1, *(suNg_vector *)s2); }),                       \
        suNg_algebra_vector *: ({                                                                                             \
            _algebra_vector_zero_g(*(suNg_algebra_vector *)r);                                                                \
            _algebra_vector_add_assign_g(*(suNg_algebra_vector *)r, *(suNg_algebra_vector *)s1);                              \
            _algebra_vector_add_assign_g(*(suNg_algebra_vector *)r, *(suNg_algebra_vector *)s2);                              \
        }),                                                                                                                   \
        double *: ({ (*(double *)r) = (*(double *)s1) + (*(double *)s2); }),                                                  \
        float *: ({ (*(double *)r) = (*(double *)s1) + (*(double *)s2); }))

#ifdef REPR_IS_REAL
#define __sub_suNfc(r, s1, s2) suNfc * : ({ _suNfc_mul_add(*(suNfc *)r, 1.0, *(suNfc *)s1, -1.0, *(suNfc *)s2); }),
#else
#define __sub_suNfc(r, s1, s2)
#endif

#define sub(r, s1, s2)                                                                                                        \
    _Generic((s2),                                                                                                            \
        spinor_field *: sub_spinor_field((spinor_field *)r, (spinor_field *)s1, (spinor_field *)s2),                          \
        spinor_field_flt *: sub_spinor_field_flt((spinor_field_flt *)r, (spinor_field_flt *)s1, (spinor_field_flt *)s2),      \
        scalar_field *: sub_scalar_field((scalar_field *)r, (scalar_field *)s1, (scalar_field *)s2),                          \
        suNg_field *: sub_suNg_field((suNg_field *)r, (suNg_field *)s1, (suNg_field *)s2),                                    \
        suNf_field *: sub_suNf_field((suNf_field *)r, (suNf_field *)s1, (suNf_field *)s2),                                    \
        suNfc_field *: sub_suNfc_field((suNfc_field *)r, (suNfc_field *)s1, (suNfc_field *)s2),                               \
        suNg_field_flt *: sub_suNg_field_flt((suNg_field_flt *)r, (suNg_field_flt *)s1, (suNg_field_flt *)s2),                \
        suNf_field_flt *: sub_suNf_field_flt((suNf_field_flt *)r, (suNf_field_flt *)s1, (suNf_field_flt *)s2),                \
        suNg_scalar_field *: sub_suNg_scalar_field((suNg_scalar_field *)r, (suNg_scalar_field *)s1, (suNg_scalar_field *)s2), \
        suNg_av_field *: sub_suNg_av_field((suNg_av_field *)r, (suNg_av_field *)s1, (suNg_av_field *)s2),                     \
        gtransf *: sub_gtransf((gtransf *)r, (gtransf *)s1, (gtransf *)s2),                                                   \
        clover_term *: sub_clover_term((clover_term *)r, (clover_term *)s1, (clover_term *)s2),                               \
        clover_force *: sub_clover_force((clover_force *)r, (clover_force *)s1, (clover_force *)s2),                          \
        staple_field *: sub_staple_field((staple_field *)r, (staple_field *)s1, (staple_field *)s2),                          \
        suNf_spinor *: ({ _spinor_sub_f(*(suNf_spinor *)r, *(suNf_spinor *)s1, *(suNf_spinor *)s2); }),                       \
        suNf_spinor_flt *: ({ _spinor_sub_f(*(suNf_spinor_flt *)r, *(suNf_spinor_flt *)s1, *(suNf_spinor_flt *)s2); }),       \
        __sub_suNfc(r, s1, s2) suNf *: ({ _suNf_mul_add(*(suNf *)r, 1.0, *(suNf *)s1, -1.0, *(suNf *)s2); }),                 \
        suNg *: ({ _suNg_mul_add(*(suNg *)r, 1.0, *(suNg *)s1, -1.0, *(suNg *)s2); }),                                        \
        suNf_flt *: ({ _suNf_mul_add(*(suNf_flt *)r, (float)(1.0), *(suNf_flt *)s1, -(float)(1.0), *(suNf_flt *)s2); }),      \
        suNg_flt *: ({ _suNg_mul_add(*(suNg_flt *)r, (float)(1.0), *(suNg_flt *)s1, -(float)(1.0), *(suNg_flt *)s2); }),      \
        suNf_vector *: ({ _vector_sub_f(*(suNf_vector *)r, *(suNf_vector *)s1, *(suNf_vector *)s2); }),                       \
        suNg_vector *: ({ _vector_sub_g(*(suNg_vector *)r, *(suNg_vector *)s1, *(suNg_vector *)s2); }),                       \
        suNg_algebra_vector *: ({                                                                                             \
            _algebra_vector_zero_g(*(suNg_algebra_vector *)r);                                                                \
            _algebra_vector_add_assign_g(*(suNg_algebra_vector *)r, *(suNg_algebra_vector *)s1);                              \
            _algebra_vector_sub_assign_g(*(suNg_algebra_vector *)r, *(suNg_algebra_vector *)s2);                              \
        }),                                                                                                                   \
        double *: ({ (*(double *)r) = (*(double *)s1) - (*(double *)s2); }),                                                  \
        float *: ({ (*(double *)r) = (*(double *)s1) - (*(double *)s2); }))

#ifdef REPR_IS_REAL
#define __sub_assign_suNfc(s1, s2) suNfc * : ({ _suNfc_sub_assign(*(suNfc *)s1, *(suNfc *)s2); }),
#else
#define __sub_assign_suNfc(s1, s2)
#endif

#define sub_assign(s1, s2)                                                                                                  \
    _Generic((s2),                                                                                                          \
        spinor_field *: sub_assign_spinor_field((spinor_field *)s1, (spinor_field *)s2),                                    \
        spinor_field_flt *: sub_assign_spinor_field_flt((spinor_field_flt *)s1, (spinor_field_flt *)s2),                    \
        scalar_field *: sub_assign_scalar_field((scalar_field *)s1, (scalar_field *)s2),                                    \
        suNg_field *: sub_assign_suNg_field((suNg_field *)s1, (suNg_field *)s2),                                            \
        suNf_field *: sub_assign_suNf_field((suNf_field *)s1, (suNf_field *)s2),                                            \
        suNfc_field *: sub_assign_suNfc_field((suNfc_field *)s1, (suNfc_field *)s2),                                        \
        suNg_field_flt *: sub_assign_suNg_field_flt((suNg_field_flt *)s1, (suNg_field_flt *)s2),                            \
        suNf_field_flt *: sub_assign_suNf_field_flt((suNf_field_flt *)s1, (suNf_field_flt *)s2),                            \
        suNg_scalar_field *: sub_assign_suNg_scalar_field((suNg_scalar_field *)s1, (suNg_scalar_field *)s2),                \
        suNg_av_field *: sub_assign_suNg_av_field((suNg_av_field *)s1, (suNg_av_field *)s2),                                \
        gtransf *: sub_assign_gtransf((gtransf *)s1, (gtransf *)s2),                                                        \
        clover_term *: sub_assign_clover_term((clover_term *)s1, (clover_term *)s2),                                        \
        clover_force *: sub_assign_clover_force((clover_force *)s1, (clover_force *)s2),                                    \
        staple_field *: sub_assign_staple_field((staple_field *)s1, (staple_field *)s2),                                    \
        suNf_spinor *: ({ _spinor_sub_assign_f(*(suNf_spinor *)s1, *(suNf_spinor *)s2); }),                                 \
        suNf_spinor_flt *: ({ _spinor_sub_assign_f(*(suNf_spinor_flt *)s1, *(suNf_spinor_flt *)s2); }),                     \
        __sub_assign_suNfc(s1, s2) suNf *: ({ _suNf_sub_assign(*(suNf *)s1, *(suNf *)s2); }),                               \
        suNg *: ({ _suNg_sub_assign(*(suNg *)s1, *(suNg *)s2); }),                                                          \
        suNf_flt *: ({ _suNf_sub_assign(*(suNf_flt *)s1, *(suNf_flt *)s2); }),                                              \
        suNg_flt *: ({ _suNg_sub_assign(*(suNg_flt *)s1, *(suNg_flt *)s2); }),                                              \
        suNf_vector *: ({ _vector_sub_assign_f(*(suNf_vector *)s1, *(suNf_vector *)s2); }),                                 \
        suNg_vector *: ({ _vector_sub_assign_g(*(suNg_vector *)s1, *(suNg_vector *)s2); }),                                 \
        suNg_algebra_vector *: ({ _algebra_vector_sub_assign_g(*(suNg_algebra_vector *)s1, *(suNg_algebra_vector *)s2); }), \
        double *: ({ (*(double *)s1) -= (*(double *)s2); }),                                                                \
        float *: ({ (*(double *)s1) -= (*(double *)s2); }))

#ifdef REPR_IS_REAL
#define __minus_suNfc(s1, s2) suNfc * : ({ _suNfc_minus(*(suNfc *)s1, *(suNfc *)s2); }),
#else
#define __minus_suNfc(s1, s2)
#endif

#define minus(s1, s2)                                                                                                      \
    _Generic((s2),                                                                                                         \
        spinor_field *: minus_spinor_field((spinor_field *)s1, (spinor_field *)s2),                                        \
        spinor_field_flt *: minus_spinor_field_flt((spinor_field_flt *)s1, (spinor_field_flt *)s2),                        \
        scalar_field *: minus_scalar_field((scalar_field *)s1, (scalar_field *)s2),                                        \
        suNg_field *: minus_suNg_field((suNg_field *)s1, (suNg_field *)s2),                                                \
        suNf_field *: minus_suNf_field((suNf_field *)s1, (suNf_field *)s2),                                                \
        suNfc_field *: minus_suNfc_field((suNfc_field *)s1, (suNfc_field *)s2),                                            \
        suNg_field_flt *: minus_suNg_field_flt((suNg_field_flt *)s1, (suNg_field_flt *)s2),                                \
        suNf_field_flt *: minus_suNf_field_flt((suNf_field_flt *)s1, (suNf_field_flt *)s2),                                \
        suNg_scalar_field *: minus_suNg_scalar_field((suNg_scalar_field *)s1, (suNg_scalar_field *)s2),                    \
        suNg_av_field *: minus_suNg_av_field((suNg_av_field *)s1, (suNg_av_field *)s2),                                    \
        gtransf *: minus_gtransf((gtransf *)s1, (gtransf *)s2),                                                            \
        clover_term *: minus_clover_term((clover_term *)s1, (clover_term *)s2),                                            \
        clover_force *: minus_clover_force((clover_force *)s1, (clover_force *)s2),                                        \
        staple_field *: minus_staple_field((staple_field *)s1, (staple_field *)s2),                                        \
        suNf_spinor *: ({ _spinor_minus_f(*(suNf_spinor *)s1, *(suNf_spinor *)s2); }),                                     \
        suNf_spinor_flt *: ({ _spinor_minus_f(*(suNf_spinor_flt *)s1, *(suNf_spinor_flt *)s2); }),                         \
        __minus_suNfc(s1, s2) suNf *: ({ _suNf_minus(*(suNf *)s1, *(suNf *)s2); }),                                        \
        suNg *: ({ _suNg_minus(*(suNg *)s1, *(suNg *)s2); }),                                                              \
        suNf_flt *: ({ _suNf_minus(*(suNf_flt *)s1, *(suNf_flt *)s2); }),                                                  \
        suNg_flt *: ({ _suNg_minus(*(suNg_flt *)s1, *(suNg_flt *)s2); }),                                                  \
        suNf_vector *: ({ _vector_minus_f(*(suNf_vector *)s1, *(suNf_vector *)s2); }),                                     \
        suNg_vector *: ({ _vector_minus_g(*(suNg_vector *)s1, *(suNg_vector *)s2); }),                                     \
        suNg_algebra_vector *: ({ _algebra_vector_mul_g(*(suNg_algebra_vector *)s1, -1.0, *(suNg_algebra_vector *)s2); }), \
        double *: ({ (*(double *)s1) = -(*(double *)s2); }),                                                               \
        float *: ({ (*(double *)s1) = -(*(double *)s2); }))

#ifdef REPR_IS_REAL
#define __add_assign_suNfc(s1, s2) suNfc * : ({ _suNfc_add_assign(*(suNfc *)s1, *(suNfc *)s2); }),
#else
#define __add_assign_suNfc(s1, s2)
#endif

#define add_assign(s1, s2)                                                                                                  \
    _Generic((s2),                                                                                                          \
        spinor_field *: add_assign_spinor_field((spinor_field *)s1, (spinor_field *)s2),                                    \
        spinor_field_flt *: add_assign_spinor_field_flt((spinor_field_flt *)s1, (spinor_field_flt *)s2),                    \
        scalar_field *: add_assign_scalar_field((scalar_field *)s1, (scalar_field *)s2),                                    \
        suNg_field *: add_assign_suNg_field((suNg_field *)s1, (suNg_field *)s2),                                            \
        suNf_field *: add_assign_suNf_field((suNf_field *)s1, (suNf_field *)s2),                                            \
        suNfc_field *: add_assign_suNfc_field((suNfc_field *)s1, (suNfc_field *)s2),                                        \
        suNg_field_flt *: add_assign_suNg_field_flt((suNg_field_flt *)s1, (suNg_field_flt *)s2),                            \
        suNf_field_flt *: add_assign_suNf_field_flt((suNf_field_flt *)s1, (suNf_field_flt *)s2),                            \
        suNg_scalar_field *: add_assign_suNg_scalar_field((suNg_scalar_field *)s1, (suNg_scalar_field *)s2),                \
        suNg_av_field *: add_assign_suNg_av_field((suNg_av_field *)s1, (suNg_av_field *)s2),                                \
        gtransf *: add_assign_gtransf((gtransf *)s1, (gtransf *)s2),                                                        \
        clover_term *: add_assign_clover_term((clover_term *)s1, (clover_term *)s2),                                        \
        clover_force *: add_assign_clover_force((clover_force *)s1, (clover_force *)s2),                                    \
        staple_field *: add_assign_staple_field((staple_field *)s1, (staple_field *)s2),                                    \
        suNf_spinor *: ({ _spinor_add_assign_f(*(suNf_spinor *)s1, *(suNf_spinor *)s2); }),                                 \
        suNf_spinor_flt *: ({ _spinor_add_assign_f(*(suNf_spinor_flt *)s1, *(suNf_spinor_flt *)s2); }),                     \
        __add_assign_suNfc(s1, s2) suNf *: ({ _suNf_add_assign(*(suNf *)s1, *(suNf *)s2); }),                               \
        suNg *: ({ _suNg_add_assign(*(suNg *)s1, *(suNg *)s2); }),                                                          \
        suNf_flt *: ({ _suNf_add_assign(*(suNf_flt *)s1, *(suNf_flt *)s2); }),                                              \
        suNg_flt *: ({ _suNg_add_assign(*(suNg_flt *)s1, *(suNg_flt *)s2); }),                                              \
        suNf_vector *: ({ _vector_add_assign_f(*(suNf_vector *)s1, *(suNf_vector *)s2); }),                                 \
        suNg_vector *: ({ _vector_add_assign_g(*(suNg_vector *)s1, *(suNg_vector *)s2); }),                                 \
        suNg_algebra_vector *: ({ _algebra_vector_add_assign_g(*(suNg_algebra_vector *)s1, *(suNg_algebra_vector *)s2); }), \
        double *: ({ (*(double *)s1) += (*(double *)s2); }),                                                                \
        float *: ({ (*(double *)s1) += (*(double *)s2); }))

double spinor_prod_re_f(suNf_spinor *r, suNf_spinor *s);
double spinor_prod_re_f_flt(suNf_spinor_flt *r, suNf_spinor_flt *s);
double suNg_prod_re(suNg *u, suNg *v);
double suNf_prod_re(suNf *u, suNf *v);
double suNg_flt_prod_re(suNg_flt *u, suNg_flt *v);
double suNf_flt_prod_re(suNf_flt *u, suNf_flt *v);
double suNgc_prod_re(suNgc *u, suNgc *v);
double suNfc_prod_re(suNfc *u, suNfc *v);
double suNf_vector_prod_re(suNf_vector *r, suNf_vector *s);
double suNg_vector_prod_re(suNg_vector *r, suNg_vector *s);
double suNg_algebra_vector_prod_re(suNg_algebra_vector *r, suNg_algebra_vector *s);

#ifdef REPR_IS_REAL
#define __prod_re_suNfc(s1, s2) suNfc * : suNfc_prod_re((suNfc *)s1, (suNfc *)s2),
#else
#define __prod_re_suNfc(s1, s2)
#endif

#define prod_re(s1, s2)                                                                                           \
    _Generic((s2),                                                                                                \
        spinor_field *: prod_re_spinor_field((spinor_field *)s1, (spinor_field *)s2),                             \
        spinor_field_flt *: prod_re_spinor_field_flt((spinor_field_flt *)s1, (spinor_field_flt *)s2),             \
        scalar_field *: prod_re_scalar_field((scalar_field *)s1, (scalar_field *)s2),                             \
        suNg_field *: prod_re_suNg_field((suNg_field *)s1, (suNg_field *)s2),                                     \
        suNf_field *: prod_re_suNf_field((suNf_field *)s1, (suNf_field *)s2),                                     \
        suNfc_field *: prod_re_suNfc_field((suNfc_field *)s1, (suNfc_field *)s2),                                 \
        suNg_field_flt *: prod_re_suNg_field_flt((suNg_field_flt *)s1, (suNg_field_flt *)s2),                     \
        suNf_field_flt *: prod_re_suNf_field_flt((suNf_field_flt *)s1, (suNf_field_flt *)s2),                     \
        suNg_scalar_field *: prod_re_suNg_scalar_field((suNg_scalar_field *)s1, (suNg_scalar_field *)s2),         \
        suNg_av_field *: prod_re_suNg_av_field((suNg_av_field *)s1, (suNg_av_field *)s2),                         \
        gtransf *: prod_re_gtransf((gtransf *)s1, (gtransf *)s2),                                                 \
        clover_term *: prod_re_clover_term((clover_term *)s1, (clover_term *)s2),                                 \
        clover_force *: prod_re_clover_force((clover_force *)s1, (clover_force *)s2),                             \
        staple_field *: prod_re_staple_field((staple_field *)s1, (staple_field *)s2),                             \
        suNf_spinor *: spinor_prod_re_f((suNf_spinor *)s1, (suNf_spinor *)s2),                                    \
        suNf_spinor_flt *: spinor_prod_re_f_flt((suNf_spinor_flt *)s1, (suNf_spinor_flt *)s2),                    \
        __prod_re_suNfc(s1, s2) suNf *: suNf_prod_re((suNf *)s1, (suNf *)s2),                                     \
        suNg *: suNg_prod_re((suNg *)s1, (suNg *)s2),                                                             \
        suNf_flt *: suNf_flt_prod_re((suNf_flt *)s1, (suNf_flt *)s2),                                             \
        suNg_flt *: suNg_flt_prod_re((suNg_flt *)s1, (suNg_flt *)s2),                                             \
        suNf_vector *: suNf_vector_prod_re((suNf_vector *)s1, (suNf_vector *)s2),                                 \
        suNg_vector *: suNg_vector_prod_re((suNg_vector *)s1, (suNg_vector *)s2),                                 \
        suNg_algebra_vector *: suNg_algebra_vector_prod_re((suNg_algebra_vector *)s1, (suNg_algebra_vector *)s2), \
        double *: ({ (*(double *)s1) * (*(double *)s2); }),                                                       \
        float *: ({ (*(float *)s1) * (*(float *)s2); }))

double spinor_prod_im_f(suNf_spinor *r, suNf_spinor *s);
double spinor_prod_im_f_flt(suNf_spinor_flt *r, suNf_spinor_flt *s);
double suNg_prod_im(suNg *u, suNg *v);
double suNf_prod_im(suNf *u, suNf *v);
double suNg_flt_prod_im(suNg_flt *u, suNg_flt *v);
double suNf_flt_prod_im(suNf_flt *u, suNf_flt *v);
double suNgc_prod_im(suNgc *u, suNgc *v);
double suNfc_prod_im(suNfc *u, suNfc *v);
double suNf_vector_prod_im(suNf_vector *r, suNf_vector *s);
double suNg_vector_prod_im(suNg_vector *r, suNg_vector *s);
double suNg_algebra_vector_prod_im(suNg_algebra_vector *r, suNg_algebra_vector *s);

#ifdef REPR_IS_REAL
#define __prod_im_suNfc(s1, s2) suNfc * : suNfc_prod_im((suNfc *)s1, (suNfc *)s2),
#else
#define __prod_im_suNfc(s1, s2)
#endif

#define prod_im(s1, s2)                                                                                           \
    _Generic((s2),                                                                                                \
        spinor_field *: prod_im_spinor_field((spinor_field *)s1, (spinor_field *)s2),                             \
        spinor_field_flt *: prod_im_spinor_field_flt((spinor_field_flt *)s1, (spinor_field_flt *)s2),             \
        scalar_field *: prod_im_scalar_field((scalar_field *)s1, (scalar_field *)s2),                             \
        suNg_field *: prod_im_suNg_field((suNg_field *)s1, (suNg_field *)s2),                                     \
        suNf_field *: prod_im_suNf_field((suNf_field *)s1, (suNf_field *)s2),                                     \
        suNfc_field *: prod_im_suNfc_field((suNfc_field *)s1, (suNfc_field *)s2),                                 \
        suNg_field_flt *: prod_im_suNg_field_flt((suNg_field_flt *)s1, (suNg_field_flt *)s2),                     \
        suNf_field_flt *: prod_im_suNf_field_flt((suNf_field_flt *)s1, (suNf_field_flt *)s2),                     \
        suNg_scalar_field *: prod_im_suNg_scalar_field((suNg_scalar_field *)s1, (suNg_scalar_field *)s2),         \
        suNg_av_field *: prod_im_suNg_av_field((suNg_av_field *)s1, (suNg_av_field *)s2),                         \
        gtransf *: prod_im_gtransf((gtransf *)s1, (gtransf *)s2),                                                 \
        clover_term *: prod_im_clover_term((clover_term *)s1, (clover_term *)s2),                                 \
        clover_force *: prod_im_clover_force((clover_force *)s1, (clover_force *)s2),                             \
        staple_field *: prod_im_staple_field((staple_field *)s1, (staple_field *)s2),                             \
        suNf_spinor *: spinor_prod_im_f((suNf_spinor *)s1, (suNf_spinor *)s2),                                    \
        suNf_spinor_flt *: spinor_prod_im_f_flt((suNf_spinor_flt *)s1, (suNf_spinor_flt *)s2),                    \
        __prod_im_suNfc(s1, s2) suNf *: suNf_prod_re((suNf *)s1, (suNf *)s2),                                     \
        suNg *: suNg_prod_re((suNg *)s1, (suNg *)s2),                                                             \
        suNf_flt *: suNf_flt_prod_re((suNf_flt *)s1, (suNf_flt *)s2),                                             \
        suNg_flt *: suNg_flt_prod_re((suNg_flt *)s1, (suNg_flt *)s2),                                             \
        suNf_vector *: suNf_vector_prod_im((suNf_vector *)s1, (suNf_vector *)s2),                                 \
        suNg_vector *: suNg_vector_prod_im((suNg_vector *)s1, (suNg_vector *)s2),                                 \
        suNg_algebra_vector *: suNg_algebra_vector_prod_im((suNg_algebra_vector *)s1, (suNg_algebra_vector *)s2), \
        double *: ({ 0.0; }),                                                                                     \
        float *: ({ 0.0; }))

hr_complex spinor_prod_f(suNf_spinor *r, suNf_spinor *s);
hr_complex spinor_prod_f_flt(suNf_spinor_flt *r, suNf_spinor_flt *s);
hr_complex suNg_prod(suNg *u, suNg *v);
hr_complex suNf_prod(suNf *u, suNf *v);
hr_complex suNg_flt_prod(suNg_flt *u, suNg_flt *v);
hr_complex suNf_flt_prod(suNf_flt *u, suNf_flt *v);
hr_complex suNgc_prod(suNgc *u, suNgc *v);
hr_complex suNfc_prod(suNfc *u, suNfc *v);
hr_complex suNf_vector_prod(suNf_vector *r, suNf_vector *s);
hr_complex suNg_vector_prod(suNg_vector *r, suNg_vector *s);
hr_complex suNg_algebra_vector_prod(suNg_algebra_vector *r, suNg_algebra_vector *s);

#ifdef REPR_IS_REAL
#define __prod_suNfc(s1, s2) suNfc * : suNfc_prod_re((suNfc *)s1, (suNfc *)s2),
#else
#define __prod_suNfc(s1, s2)
#endif

#define prod(s1, s2)                                                                                           \
    _Generic((s2),                                                                                             \
        spinor_field *: prod_spinor_field((spinor_field *)s1, (spinor_field *)s2),                             \
        spinor_field_flt *: ({                                                                                 \
            prod_re_spinor_field_flt((spinor_field_flt *)s1, (spinor_field_flt *)s2) +                         \
                I *prod_im_spinor_field_flt((spinor_field_flt *)s1, (spinor_field_flt *)s2);                   \
        }),                                                                                                    \
        scalar_field *: prod_scalar_field((scalar_field *)s1, (scalar_field *)s2),                             \
        suNg_field *: prod_suNg_field((suNg_field *)s1, (suNg_field *)s2),                                     \
        suNf_field *: prod_suNf_field((suNf_field *)s1, (suNf_field *)s2),                                     \
        suNfc_field *: prod_suNfc_field((suNfc_field *)s1, (suNfc_field *)s2),                                 \
        suNg_field_flt *: prod_suNg_field_flt((suNg_field_flt *)s1, (suNg_field_flt *)s2),                     \
        suNf_field_flt *: prod_suNf_field_flt((suNf_field_flt *)s1, (suNf_field_flt *)s2),                     \
        suNg_scalar_field *: prod_suNg_scalar_field((suNg_scalar_field *)s1, (suNg_scalar_field *)s2),         \
        suNg_av_field *: prod_suNg_av_field((suNg_av_field *)s1, (suNg_av_field *)s2),                         \
        gtransf *: prod_gtransf((gtransf *)s1, (gtransf *)s2),                                                 \
        clover_term *: prod_clover_term((clover_term *)s1, (clover_term *)s2),                                 \
        clover_force *: prod_clover_force((clover_force *)s1, (clover_force *)s2),                             \
        staple_field *: prod_staple_field((staple_field *)s1, (staple_field *)s2),                             \
        suNf_spinor *: spinor_prod_f((suNf_spinor *)s1, (suNf_spinor *)s2),                                    \
        suNf_spinor_flt *: spinor_prod_f_flt((suNf_spinor_flt *)s1, (suNf_spinor_flt *)s2),                    \
        __prod_suNfc(s1, s2) suNf *: suNf_prod((suNf *)s1, (suNf *)s2),                                        \
        suNg *: suNg_prod((suNg *)s1, (suNg *)s2),                                                             \
        suNf_flt *: suNf_flt_prod((suNf_flt *)s1, (suNf_flt *)s2),                                             \
        suNg_flt *: suNg_flt_prod((suNg_flt *)s1, (suNg_flt *)s2),                                             \
        suNf_vector *: suNf_vector_prod((suNf_vector *)s1, (suNf_vector *)s2),                                 \
        suNg_vector *: suNg_vector_prod((suNg_vector *)s1, (suNg_vector *)s2),                                 \
        suNg_algebra_vector *: suNg_algebra_vector_prod((suNg_algebra_vector *)s1, (suNg_algebra_vector *)s2), \
        double *: ({ (*(double *)s1) * (*(double *)s2); }),                                                    \
        float *: ({ (*(float *)s1) * (*(float *)s2); }))

double spinor_sqnorm_f(suNf_spinor *r);
double spinor_sqnorm_f_flt(suNf_spinor_flt *r);
double suNg_sqnorm(suNg *u);
double suNf_sqnorm(suNf *u);
double suNg_flt_sqnorm(suNg_flt *u);
double suNf_flt_sqnorm(suNf_flt *u);
double suNgc_sqnorm(suNgc *u);
double suNfc_sqnorm(suNfc *u);
double suNf_vector_sqnorm(suNf_vector *r);
double suNg_vector_sqnorm(suNg_vector *r);
double suNg_algebra_vector_sqnorm(suNg_algebra_vector *r);

#ifdef REPR_IS_REAL
#define __sqnorm_suNfc(s1) suNfc * : suNfc_sqnorm((suNfc *)s1),
#else
#define __sqnorm_suNfc(s1)
#endif

#define sqnorm(s1)                                                                    \
    _Generic((s1),                                                                    \
        spinor_field *: sqnorm_spinor_field((spinor_field *)s1),                      \
        spinor_field_flt *: sqnorm_spinor_field_flt((spinor_field_flt *)s1),          \
        scalar_field *: sqnorm_scalar_field((scalar_field *)s1),                      \
        suNg_field *: sqnorm_suNg_field((suNg_field *)s1),                            \
        suNf_field *: sqnorm_suNf_field((suNf_field *)s1),                            \
        suNfc_field *: sqnorm_suNfc_field((suNfc_field *)s1),                         \
        suNg_field_flt *: sqnorm_suNg_field_flt((suNg_field_flt *)s1),                \
        suNf_field_flt *: sqnorm_suNf_field_flt((suNf_field_flt *)s1),                \
        suNg_scalar_field *: sqnorm_suNg_scalar_field((suNg_scalar_field *)s1),       \
        suNg_av_field *: sqnorm_suNg_av_field((suNg_av_field *)s1),                   \
        gtransf *: sqnorm_gtransf((gtransf *)s1),                                     \
        clover_term *: sqnorm_clover_term((clover_term *)s1),                         \
        clover_force *: sqnorm_clover_force((clover_force *)s1),                      \
        staple_field *: sqnorm_staple_field((staple_field *)s1),                      \
        suNf_spinor *: spinor_sqnorm_f((suNf_spinor *)s1),                            \
        suNf_spinor_flt *: spinor_sqnorm_f_flt((suNf_spinor_flt *)s1),                \
        __sqnorm_suNfc(s1) suNf *: suNf_sqnorm((suNf *)s1),                           \
        suNg *: suNg_sqnorm((suNg *)s1),                                              \
        suNf_flt *: suNf_flt_sqnorm((suNf_flt *)s1),                                  \
        suNg_flt *: suNg_flt_sqnorm((suNg_flt *)s1),                                  \
        suNf_vector *: suNf_vector_sqnorm((suNf_vector *)s1),                         \
        suNg_vector *: suNg_vector_sqnorm((suNg_vector *)s1),                         \
        suNg_algebra_vector *: suNg_algebra_vector_sqnorm((suNg_algebra_vector *)s1), \
        double *: ({ (*(double *)s1) * (*(double *)s1); }),                           \
        float *: ({ (*(float *)s1) * (*(float *)s1); }))

double spinor_max_f(suNf_spinor *r);
double spinor_max_f_flt(suNf_spinor_flt *r);
double suNg_max(suNg *r);
double suNf_max(suNf *r);
double suNg_flt_max(suNg_flt *r);
double suNf_flt_max(suNf_flt *r);
double suNgc_max(suNgc *r);
double suNfc_max(suNfc *r);
double suNf_vector_max(suNf_vector *r);
double suNg_vector_max(suNg_vector *r);
double suNg_algebra_vector_max(suNg_algebra_vector *r);

#ifdef REPR_IS_REAL
#define __max_suNfc(s1) suNfc * : suNfc_max((suNfc *)s1),
#else
#define __max_suNfc(s1)
#endif

#define max(s1)                                                                                 \
    _Generic((s1),                                                                              \
        spinor_field *: max_spinor_field((spinor_field *)s1),                                   \
        spinor_field_flt *: max_spinor_field_flt((spinor_field_flt *)s1),                       \
        scalar_field *: max_scalar_field((scalar_field *)s1),                                   \
        suNg_field *: max_suNg_field((suNg_field *)s1),                                         \
        suNf_field *: max_suNf_field((suNf_field *)s1),                                         \
        suNfc_field *: max_suNfc_field((suNfc_field *)s1),                                      \
        suNg_field_flt *: max_suNg_field_flt((suNg_field_flt *)s1),                             \
        suNf_field_flt *: max_suNf_field_flt((suNf_field_flt *)s1),                             \
        suNg_scalar_field *: max_suNg_scalar_field((suNg_scalar_field *)s1),                    \
        suNg_av_field *: max_suNg_av_field((suNg_av_field *)s1),                                \
        gtransf *: max_gtransf((gtransf *)s1),                                                  \
        clover_term *: max_clover_term((clover_term *)s1),                                      \
        clover_force *: max_clover_force((clover_force *)s1),                                   \
        staple_field *: max_staple_field((staple_field *)s1),                                   \
        suNf_spinor *: spinor_max_f((suNf_spinor *)s1),                                         \
        suNf_spinor_flt *: spinor_max_f_flt((suNf_spinor_flt *)s1),                             \
        __max_suNfc(s1) suNf *: suNf_max((suNf *)s1),                                           \
        suNg *: suNg_max((suNg *)s1),                                                           \
        suNf_flt *: suNf_flt_max((suNf_flt *)s1),                                               \
        suNg_flt *: suNg_flt_max((suNg_flt *)s1),                                               \
        suNf_vector *: suNf_vector_max((suNf_vector *)s1),                                      \
        suNg_vector *: suNg_vector_max((suNg_vector *)s1),                                      \
        suNg_algebra_vector *: suNg_algebra_vector_max((suNg_algebra_vector *)s1),              \
        double *: ({ (*(double *)s1) > (*(double *)s1) ? (*(double *)s1) : (*(double *)s1); }), \
        float *: ({ (*(float *)s1) > (*(float *)s1) ? (*(float *)s1) : (*(float *)s1); }))

#define zero(s1)                                                              \
    _Generic((s1),                                                            \
        spinor_field *: zero_spinor_field((spinor_field *)s1),                \
        spinor_field_flt *: zero_spinor_field_flt((spinor_field_flt *)s1),    \
        scalar_field *: zero_scalar_field((scalar_field *)s1),                \
        suNg_field *: zero_suNg_field((suNg_field *)s1),                      \
        suNf_field *: zero_suNf_field((suNf_field *)s1),                      \
        suNfc_field *: zero_suNfc_field((suNfc_field *)s1),                   \
        suNg_field_flt *: zero_suNg_field_flt((suNg_field_flt *)s1),          \
        suNf_field_flt *: zero_suNf_field_flt((suNf_field_flt *)s1),          \
        suNg_scalar_field *: zero_suNg_scalar_field((suNg_scalar_field *)s1), \
        suNg_av_field *: zero_suNg_av_field((suNg_av_field *)s1),             \
        gtransf *: zero_gtransf((gtransf *)s1),                               \
        clover_term *: zero_clover_term((clover_term *)s1),                   \
        clover_force *: zero_clover_force((clover_force *)s1),                \
        staple_field *: zero_staple_field((staple_field *)s1))

#define copy(s1, s2)                                                                                   \
    _Generic((s2),                                                                                     \
        spinor_field *: copy_spinor_field((spinor_field *)s1, (spinor_field *)s2),                     \
        spinor_field_flt *: copy_spinor_field_flt((spinor_field_flt *)s1, (spinor_field_flt *)s2),     \
        scalar_field *: copy_scalar_field((scalar_field *)s1, (scalar_field *)s2),                     \
        suNg_field *: copy_suNg_field((suNg_field *)s1, (suNg_field *)s2),                             \
        suNf_field *: copy_suNf_field((suNf_field *)s1, (suNf_field *)s2),                             \
        suNfc_field *: copy_suNfc_field((suNfc_field *)s1, (suNfc_field *)s2),                         \
        suNg_field_flt *: copy_suNg_field_flt((suNg_field_flt *)s1, (suNg_field_flt *)s2),             \
        suNf_field_flt *: copy_suNf_field_flt((suNf_field_flt *)s1, (suNf_field_flt *)s2),             \
        suNg_scalar_field *: copy_suNg_scalar_field((suNg_scalar_field *)s1, (suNg_scalar_field *)s2), \
        suNg_av_field *: copy_suNg_av_field((suNg_av_field *)s1, (suNg_av_field *)s2),                 \
        gtransf *: copy_gtransf((gtransf *)s1, (gtransf *)s2),                                         \
        clover_term *: copy_clover_term((clover_term *)s1, (clover_term *)s2),                         \
        clover_force *: copy_clover_force((clover_force *)s1, (clover_force *)s2),                     \
        staple_field *: copy_staple_field((staple_field *)s1, (staple_field *)s2))

double spinor_g5_prod_re_f(suNf_spinor *r, suNf_spinor *s);
double spinor_g5_prod_re_f_flt(suNf_spinor_flt *r, suNf_spinor_flt *s);

#define g5_prod_re(s1, s2)                                                                               \
    _Generic((s2),                                                                                       \
        spinor_field *: g5_prod_re_spinor_field((spinor_field *)s1, (spinor_field *)s2),                 \
        spinor_field_flt *: g5_prod_re_spinor_field_flt((spinor_field_flt *)s1, (spinor_field_flt *)s2), \
        suNf_spinor *: spinor_g5_prod_re_f((suNf_spinor *)s1, (suNf_spinor *)s2),                        \
        suNf_spinor_flt *: spinor_g5_prod_re_f_flt((suNf_spinor_flt *)s1, (suNf_spinor_flt *)s2))

double spinor_g5_prod_im_f(suNf_spinor *r, suNf_spinor *s);
double spinor_g5_prod_im_f_flt(suNf_spinor_flt *r, suNf_spinor_flt *s);

#define g5_prod_im(s1, s2)                                                                               \
    _Generic((s2),                                                                                       \
        spinor_field *: g5_prod_im_spinor_field((spinor_field *)s1, (spinor_field *)s2),                 \
        spinor_field_flt *: g5_prod_im_spinor_field_flt((spinor_field_flt *)s1, (spinor_field_flt *)s2), \
        suNf_spinor *: spinor_g5_prod_im_f((suNf_spinor *)s1, (suNf_spinor *)s2),                        \
        suNf_spinor_flt *: spinor_g5_prod_im_f_flt((suNf_spinor_flt *)s1, (suNf_spinor_flt *)s2))

#define g5(s1, s2)                                                                               \
    _Generic((s2),                                                                               \
        spinor_field *: g5_spinor_field((spinor_field *)s1, (spinor_field *)s2),                 \
        spinor_field_flt *: g5_spinor_field_flt((spinor_field_flt *)s1, (spinor_field_flt *)s2), \
        suNf_spinor *: ({ _spinor_g5_f(*(suNf_spinor *)s1, *(suNf_spinor *)s2); }),              \
        suNf_spinor_flt *: ({ _spinor_g5_f(*(suNf_spinor_flt *)s1, *(suNf_spinor_flt *)s2); }))

#define g5_assign(s1)                                                           \
    _Generic((s1),                                                              \
        spinor_field *: g5_assign_spinor_field((spinor_field *)s1),             \
        spinor_field_flt *: g5_assign_spinor_field_flt((spinor_field_flt *)s1), \
        suNf_spinor *: ({ _spinor_g5_assign_f(*(suNf_spinor *)s1); }),          \
        suNf_spinor_flt *: ({ _spinor_g5_assign_f(*(suNf_spinor_flt *)s1); }))

#define g0(s1, s2)                                                                               \
    _Generic((s2),                                                                               \
        spinor_field *: g0_spinor_field((spinor_field *)s1, (spinor_field *)s2),                 \
        spinor_field_flt *: g0_spinor_field_flt((spinor_field_flt *)s1, (spinor_field_flt *)s2), \
        suNf_spinor *: ({ _spinor_g0_f(*(suNf_spinor *)s1, *(suNf_spinor *)s2); }),              \
        suNf_spinor_flt *: ({ _spinor_g0_f(*(suNf_spinor_flt *)s1, *(suNf_spinor_flt *)s2); }))

#define g1(s1, s2)                                                                               \
    _Generic((s2),                                                                               \
        spinor_field *: g1_spinor_field((spinor_field *)s1, (spinor_field *)s2),                 \
        spinor_field_flt *: g1_spinor_field_flt((spinor_field_flt *)s1, (spinor_field_flt *)s2), \
        suNf_spinor *: ({ _spinor_g1_f(*(suNf_spinor *)s1, *(suNf_spinor *)s2); }),              \
        suNf_spinor_flt *: ({ _spinor_g1_f(*(suNf_spinor_flt *)s1, *(suNf_spinor_flt *)s2); }))

#define g2(s1, s2)                                                                               \
    _Generic((s2),                                                                               \
        spinor_field *: g2_spinor_field((spinor_field *)s1, (spinor_field *)s2),                 \
        spinor_field_flt *: g2_spinor_field_flt((spinor_field_flt *)s1, (spinor_field_flt *)s2), \
        suNf_spinor *: ({ _spinor_g2_f(*(suNf_spinor *)s1, *(suNf_spinor *)s2); }),              \
        suNf_spinor_flt *: ({ _spinor_g2_f(*(suNf_spinor_flt *)s1, *(suNf_spinor_flt *)s2); }))

#define g3(s1, s2)                                                                               \
    _Generic((s2),                                                                               \
        spinor_field *: g3_spinor_field((spinor_field *)s1, (spinor_field *)s2),                 \
        spinor_field_flt *: g3_spinor_field_flt((spinor_field_flt *)s1, (spinor_field_flt *)s2), \
        suNf_spinor *: ({ _spinor_g3_f(*(suNf_spinor *)s1, *(suNf_spinor *)s2); }),              \
        suNf_spinor_flt *: ({ _spinor_g3_f(*(suNf_spinor_flt *)s1, *(suNf_spinor_flt *)s2); }))

#define lc(r, k1, s1, k2, s2)                                                                                  \
    _Generic((s2),                                                                                             \
        spinor_field *: lc_spinor_field((spinor_field *)r, k1, (spinor_field *)s1, k2, (spinor_field *)s2),    \
        spinor_field_flt *: lc_spinor_field_flt((spinor_field_flt *)r, k1, (spinor_field_flt *)s1, k2,         \
                                                (spinor_field_flt *)s2),                                       \
        suNf_spinor *: ({ _spinor_lc_f(*(suNf_spinor *)r, k1, *(suNf_spinor *)s1, k2, *(suNf_spinor *)s2); }), \
        suNf_spinor_flt *: ({ _spinor_lc_f(*(suNf_spinor_flt *)r, k1, *(suNf_spinor_flt *)s1, k2, *(suNf_spinor_flt *)s2); }))

#define lc_add_assign(r, k1, s1, k2, s2)                                                                                  \
    _Generic((s2),                                                                                                        \
        spinor_field *: lc_add_assign_spinor_field((spinor_field *)r, k1, (spinor_field *)s1, k2, (spinor_field *)s2),    \
        spinor_field_flt *: lc_add_assign_spinor_field_flt((spinor_field_flt *)r, k1, (spinor_field_flt *)s1, k2,         \
                                                           (spinor_field_flt *)s2),                                       \
        suNf_spinor *: ({ _spinor_lc_add_assign_f(*(suNf_spinor *)r, k1, *(suNf_spinor *)s1, k2, *(suNf_spinor *)s2); }), \
        suNf_spinor_flt *: ({                                                                                             \
            _spinor_lc_add_assign_f(*(suNf_spinor_flt *)r, k1, *(suNf_spinor_flt *)s1, k2, *(suNf_spinor_flt *)s2);       \
        }))

#define clc(r, k1, s1, k2, s2)                                                                                  \
    _Generic((s2),                                                                                              \
        spinor_field *: clc_spinor_field((spinor_field *)r, k1, (spinor_field *)s1, k2, (spinor_field *)s2),    \
        spinor_field_flt *: clc_spinor_field_flt((spinor_field_flt *)r, k1, (spinor_field_flt *)s1, k2,         \
                                                 (spinor_field_flt *)s2),                                       \
        suNf_spinor *: ({ _spinor_clc_f(*(suNf_spinor *)r, k1, *(suNf_spinor *)s1, k2, *(suNf_spinor *)s2); }), \
        suNf_spinor_flt *: ({                                                                                   \
            _spinor_clc_f(*(suNf_spinor_flt *)r, k1, *(suNf_spinor_flt *)s1, k2, *(suNf_spinor_flt *)s2);       \
        }))

#define clc_add_assign(r, k1, s1, k2, s2)                                                                                  \
    _Generic((s2),                                                                                                         \
        spinor_field *: clc_add_assign_spinor_field((spinor_field *)r, k1, (spinor_field *)s1, k2, (spinor_field *)s2),    \
        spinor_field_flt *: clc_add_assign_spinor_field_flt((spinor_field_flt *)r, k1, (spinor_field_flt *)s1, k2,         \
                                                            (spinor_field_flt *)s2),                                       \
        suNf_spinor *: ({ _spinor_clc_add_assign_f(*(suNf_spinor *)r, k1, *(suNf_spinor *)s1, k2, *(suNf_spinor *)s2); }), \
        suNf_spinor_flt *: ({                                                                                              \
            _spinor_clc_add_assign_f(*(suNf_spinor_flt *)r, k1, *(suNf_spinor_flt *)s1, k2, *(suNf_spinor_flt *)s2);       \
        }))

#ifdef WITH_GPU

#define mul_add_assign_cpu(s1, rho, s2)                                                                                   \
    _Generic((s2),                                                                                                        \
        spinor_field *: mul_add_assign_spinor_field_cpu((spinor_field *)s1, rho, (spinor_field *)s2),                     \
        spinor_field_flt *: mul_add_assign_spinor_field_flt_cpu((spinor_field_flt *)s1, rho, (spinor_field_flt *)s2),     \
        scalar_field *: mul_add_assign_scalar_field_cpu((scalar_field *)s1, rho, (scalar_field *)s2),                     \
        suNg_field *: mul_add_assign_suNg_field_cpu((suNg_field *)s1, rho, (suNg_field *)s2),                             \
        suNf_field *: mul_add_assign_suNf_field_cpu((suNf_field *)s1, rho, (suNf_field *)s2),                             \
        suNfc_field *: mul_add_assign_suNfc_field_cpu((suNfc_field *)s1, rho, (suNfc_field *)s2),                         \
        suNg_field_flt *: mul_add_assign_suNg_field_flt_cpu((suNg_field_flt *)s1, rho, (suNg_field_flt *)s2),             \
        suNf_field_flt *: mul_add_assign_suNf_field_flt_cpu((suNf_field_flt *)s1, rho, (suNf_field_flt *)s2),             \
        suNg_scalar_field *: mul_add_assign_suNg_scalar_field_cpu((suNg_scalar_field *)s1, rho, (suNg_scalar_field *)s2), \
        suNg_av_field *: mul_add_assign_suNg_av_field_cpu((suNg_av_field *)s1, rho, (suNg_av_field *)s2),                 \
        gtransf *: mul_add_assign_gtransf_cpu((gtransf *)s1, rho, (gtransf *)s2),                                         \
        clover_term *: mul_add_assign_clover_term_cpu((clover_term *)s1, rho, (clover_term *)s2),                         \
        clover_force *: mul_add_assign_clover_force_cpu((clover_force *)s1, rho, (clover_force *)s2),                     \
        staple_field *: mul_add_assign_staple_field_cpu((staple_field *)s1, rho, (staple_field *)s2))

#define mulc_add_assign_cpu(s1, rho, s2)                                                                                   \
    _Generic((s2),                                                                                                         \
        spinor_field *: mulc_add_assign_spinor_field_cpu((spinor_field *)s1, rho, (spinor_field *)s2),                     \
        spinor_field_flt *: mulc_add_assign_spinor_field_flt_cpu((spinor_field_flt *)s1, rho, (spinor_field_flt *)s2),     \
        scalar_field *: mulc_add_assign_scalar_field_cpu((scalar_field *)s1, rho, (scalar_field *)s2),                     \
        suNg_field *: mulc_add_assign_suNg_field_cpu((suNg_field *)s1, rho, (suNg_field *)s2),                             \
        suNf_field *: mulc_add_assign_suNf_field_cpu((suNf_field *)s1, rho, (suNf_field *)s2),                             \
        suNfc_field *: mulc_add_assign_suNfc_field_cpu((suNfc_field *)s1, rho, (suNfc_field *)s2),                         \
        suNg_field_flt *: mulc_add_assign_suNg_field_flt_cpu((suNg_field_flt *)s1, rho, (suNg_field_flt *)s2),             \
        suNf_field_flt *: mulc_add_assign_suNf_field_flt_cpu((suNf_field_flt *)s1, rho, (suNf_field_flt *)s2),             \
        suNg_scalar_field *: mulc_add_assign_suNg_scalar_field_cpu((suNg_scalar_field *)s1, rho, (suNg_scalar_field *)s2), \
        suNg_av_field *: mulc_add_assign_suNg_av_field_cpu((suNg_av_field *)s1, rho, (suNg_av_field *)s2),                 \
        gtransf *: mulc_add_assign_gtransf_cpu((gtransf *)s1, rho, (gtransf *)s2),                                         \
        clover_term *: mulc_add_assign_clover_term_cpu((clover_term *)s1, rho, (clover_term *)s2),                         \
        clover_force *: mulc_add_assign_clover_force_cpu((clover_force *)s1, rho, (clover_force *)s2),                     \
        staple_field *: mulc_add_assign_staple_field_cpu((staple_field *)s1, rho, (staple_field *)s2))

#define mul_cpu(s1, rho, s2)                                                                                   \
    _Generic((s2),                                                                                             \
        spinor_field *: mul_spinor_field_cpu((spinor_field *)s1, rho, (spinor_field *)s2),                     \
        spinor_field_flt *: mul_spinor_field_flt_cpu((spinor_field_flt *)s1, rho, (spinor_field_flt *)s2),     \
        scalar_field *: mul_scalar_field_cpu((scalar_field *)s1, rho, (scalar_field *)s2),                     \
        suNg_field *: mul_suNg_field_cpu((suNg_field *)s1, rho, (suNg_field *)s2),                             \
        suNf_field *: mul_suNf_field_cpu((suNf_field *)s1, rho, (suNf_field *)s2),                             \
        suNfc_field *: mul_suNfc_field_cpu((suNfc_field *)s1, rho, (suNfc_field *)s2),                         \
        suNg_field_flt *: mul_suNg_field_flt_cpu((suNg_field_flt *)s1, rho, (suNg_field_flt *)s2),             \
        suNf_field_flt *: mul_suNf_field_flt_cpu((suNf_field_flt *)s1, rho, (suNf_field_flt *)s2),             \
        suNg_scalar_field *: mul_suNg_scalar_field_cpu((suNg_scalar_field *)s1, rho, (suNg_scalar_field *)s2), \
        suNg_av_field *: mul_suNg_av_field_cpu((suNg_av_field *)s1, rho, (suNg_av_field *)s2),                 \
        gtransf *: mul_gtransf_cpu((gtransf *)s1, rho, (gtransf *)s2),                                         \
        clover_term *: mul_clover_term_cpu((clover_term *)s1, rho, (clover_term *)s2),                         \
        clover_force *: mul_clover_force_cpu((clover_force *)s1, rho, (clover_force *)s2),                     \
        staple_field *: mul_staple_field_cpu((staple_field *)s1, rho, (staple_field *)s2))

#define mulc_cpu(s1, rho, s2)                                                                                   \
    _Generic((s2),                                                                                              \
        spinor_field *: mulc_spinor_field_cpu((spinor_field *)s1, rho, (spinor_field *)s2),                     \
        spinor_field_flt *: mulc_spinor_field_flt_cpu((spinor_field_flt *)s1, rho, (spinor_field_flt *)s2),     \
        scalar_field *: mulc_scalar_field_cpu((scalar_field *)s1, rho, (scalar_field *)s2),                     \
        suNg_field *: mulc_suNg_field_cpu((suNg_field *)s1, rho, (suNg_field *)s2),                             \
        suNf_field *: mulc_suNf_field_cpu((suNf_field *)s1, rho, (suNf_field *)s2),                             \
        suNfc_field *: mulc_suNfc_field_cpu((suNfc_field *)s1, rho, (suNfc_field *)s2),                         \
        suNg_field_flt *: mulc_suNg_field_flt_cpu((suNg_field_flt *)s1, rho, (suNg_field_flt *)s2),             \
        suNf_field_flt *: mulc_suNf_field_flt_cpu((suNf_field_flt *)s1, rho, (suNf_field_flt *)s2),             \
        suNg_scalar_field *: mulc_suNg_scalar_field_cpu((suNg_scalar_field *)s1, rho, (suNg_scalar_field *)s2), \
        suNg_av_field *: mulc_suNg_av_field_cpu((suNg_av_field *)s1, rho, (suNg_av_field *)s2),                 \
        gtransf *: mulc_gtransf_cpu((gtransf *)s1, rho, (gtransf *)s2),                                         \
        clover_term *: mulc_clover_term_cpu((clover_term *)s1, rho, (clover_term *)s2),                         \
        clover_force *: mulc_clover_force_cpu((clover_force *)s1, rho, (clover_force *)s2),                     \
        staple_field *: mulc_staple_field_cpu((staple_field *)s1, rho, (staple_field *)s2))

#define add_cpu(r, s1, s2)                                                                                                   \
    _Generic((s2),                                                                                                           \
        spinor_field *: add_spinor_field_cpu((spinor_field *)r, (spinor_field *)s1, (spinor_field *)s2),                     \
        spinor_field_flt *: add_spinor_field_flt_cpu((spinor_field_flt *)r, (spinor_field_flt *)s1, (spinor_field_flt *)s2), \
        scalar_field *: add_scalar_field_cpu((scalar_field *)r, (scalar_field *)s1, (scalar_field *)s2),                     \
        suNg_field *: add_suNg_field_cpu((suNg_field *)r, (suNg_field *)s1, (suNg_field *)s2),                               \
        suNf_field *: add_suNf_field_cpu((suNf_field *)r, (suNf_field *)s1, (suNf_field *)s2),                               \
        suNfc_field *: add_suNfc_field_cpu((suNfc_field *)r, (suNfc_field *)s1, (suNfc_field *)s2),                          \
        suNg_field_flt *: add_suNg_field_flt_cpu((suNg_field_flt *)r, (suNg_field_flt *)s1, (suNg_field_flt *)s2),           \
        suNf_field_flt *: add_suNf_field_flt_cpu((suNf_field_flt *)r, (suNf_field_flt *)s1, (suNf_field_flt *)s2),           \
        suNg_scalar_field *: add_suNg_scalar_field_cpu((suNg_scalar_field *)r, (suNg_scalar_field *)s1,                      \
                                                       (suNg_scalar_field *)s2),                                             \
        suNg_av_field *: add_suNg_av_field_cpu((suNg_av_field *)r, (suNg_av_field *)s1, (suNg_av_field *)s2),                \
        gtransf *: add_gtransf_cpu((gtransf *)r, (gtransf *)s1, (gtransf *)s2),                                              \
        clover_term *: add_clover_term_cpu((clover_term *)r, (clover_term *)s1, (clover_term *)s2),                          \
        clover_force *: add_clover_force_cpu((clover_force *)r, (clover_force *)s1, (clover_force *)s2),                     \
        staple_field *: add_staple_field_cpu((staple_field *)r, (staple_field *)s1, (staple_field *)s2))

#define sub_cpu(r, s1, s2)                                                                                                   \
    _Generic((s2),                                                                                                           \
        spinor_field *: sub_spinor_field_cpu((spinor_field *)r, (spinor_field *)s1, (spinor_field *)s2),                     \
        spinor_field_flt *: sub_spinor_field_flt_cpu((spinor_field_flt *)r, (spinor_field_flt *)s1, (spinor_field_flt *)s2), \
        scalar_field *: sub_scalar_field_cpu((scalar_field *)r, (scalar_field *)s1, (scalar_field *)s2),                     \
        suNg_field *: sub_suNg_field_cpu((suNg_field *)r, (suNg_field *)s1, (suNg_field *)s2),                               \
        suNf_field *: sub_suNf_field_cpu((suNf_field *)r, (suNf_field *)s1, (suNf_field *)s2),                               \
        suNfc_field *: sub_suNfc_field_cpu((suNfc_field *)r, (suNfc_field *)s1, (suNfc_field *)s2),                          \
        suNg_field_flt *: sub_suNg_field_flt_cpu((suNg_field_flt *)r, (suNg_field_flt *)s1, (suNg_field_flt *)s2),           \
        suNf_field_flt *: sub_suNf_field_flt_cpu((suNf_field_flt *)r, (suNf_field_flt *)s1, (suNf_field_flt *)s2),           \
        suNg_scalar_field *: sub_suNg_scalar_field_cpu((suNg_scalar_field *)r, (suNg_scalar_field *)s1,                      \
                                                       (suNg_scalar_field *)s2),                                             \
        suNg_av_field *: sub_suNg_av_field_cpu((suNg_av_field *)r, (suNg_av_field *)s1, (suNg_av_field *)s2),                \
        gtransf *: sub_gtransf_cpu((gtransf *)r, (gtransf *)s1, (gtransf *)s2),                                              \
        clover_term *: sub_clover_term_cpu((clover_term *)r, (clover_term *)s1, (clover_term *)s2),                          \
        clover_force *: sub_clover_force_cpu((clover_force *)r, (clover_force *)s1, (clover_force *)s2),                     \
        staple_field *: sub_staple_field_cpu((staple_field *)r, (staple_field *)s1, (staple_field *)s2))

#define sub_assign_cpu(s1, s2)                                                                                   \
    _Generic((s2),                                                                                               \
        spinor_field *: sub_assign_spinor_field_cpu((spinor_field *)s1, (spinor_field *)s2),                     \
        spinor_field_flt *: sub_assign_spinor_field_flt_cpu((spinor_field_flt *)s1, (spinor_field_flt *)s2),     \
        scalar_field *: sub_assign_scalar_field_cpu((scalar_field *)s1, (scalar_field *)s2),                     \
        suNg_field *: sub_assign_suNg_field_cpu((suNg_field *)s1, (suNg_field *)s2),                             \
        suNf_field *: sub_assign_suNf_field_cpu((suNf_field *)s1, (suNf_field *)s2),                             \
        suNfc_field *: sub_assign_suNfc_field_cpu((suNfc_field *)s1, (suNfc_field *)s2),                         \
        suNg_field_flt *: sub_assign_suNg_field_flt_cpu((suNg_field_flt *)s1, (suNg_field_flt *)s2),             \
        suNf_field_flt *: sub_assign_suNf_field_flt_cpu((suNf_field_flt *)s1, (suNf_field_flt *)s2),             \
        suNg_scalar_field *: sub_assign_suNg_scalar_field_cpu((suNg_scalar_field *)s1, (suNg_scalar_field *)s2), \
        suNg_av_field *: sub_assign_suNg_av_field_cpu((suNg_av_field *)s1, (suNg_av_field *)s2),                 \
        gtransf *: sub_assign_gtransf_cpu((gtransf *)s1, (gtransf *)s2),                                         \
        clover_term *: sub_assign_clover_term_cpu((clover_term *)s1, (clover_term *)s2),                         \
        clover_force *: sub_assign_clover_force_cpu((clover_force *)s1, (clover_force *)s2),                     \
        staple_field *: sub_assign_staple_field_cpu((staple_field *)s1, (staple_field *)s2))

#define minus_cpu(s1, s2)                                                                                   \
    _Generic((s2),                                                                                          \
        spinor_field *: minus_spinor_field_cpu((spinor_field *)s1, (spinor_field *)s2),                     \
        spinor_field_flt *: minus_spinor_field_flt_cpu((spinor_field_flt *)s1, (spinor_field_flt *)s2),     \
        scalar_field *: minus_scalar_field_cpu((scalar_field *)s1, (scalar_field *)s2),                     \
        suNg_field *: minus_suNg_field_cpu((suNg_field *)s1, (suNg_field *)s2),                             \
        suNf_field *: minus_suNf_field_cpu((suNf_field *)s1, (suNf_field *)s2),                             \
        suNfc_field *: minus_suNfc_field_cpu((suNfc_field *)s1, (suNfc_field *)s2),                         \
        suNg_field_flt *: minus_suNg_field_flt_cpu((suNg_field_flt *)s1, (suNg_field_flt *)s2),             \
        suNf_field_flt *: minus_suNf_field_flt_cpu((suNf_field_flt *)s1, (suNf_field_flt *)s2),             \
        suNg_scalar_field *: minus_suNg_scalar_field_cpu((suNg_scalar_field *)s1, (suNg_scalar_field *)s2), \
        suNg_av_field *: minus_suNg_av_field_cpu((suNg_av_field *)s1, (suNg_av_field *)s2),                 \
        gtransf *: minus_gtransf_cpu((gtransf *)s1, (gtransf *)s2),                                         \
        clover_term *: minus_clover_term_cpu((clover_term *)s1, (clover_term *)s2),                         \
        clover_force *: minus_clover_force_cpu((clover_force *)s1, (clover_force *)s2),                     \
        staple_field *: minus_staple_field_cpu((staple_field *)s1, (staple_field *)s2))

#define add_assign_cpu(s1, s2)                                                                                   \
    _Generic((s2),                                                                                               \
        spinor_field *: add_assign_spinor_field_cpu((spinor_field *)s1, (spinor_field *)s2),                     \
        spinor_field_flt *: add_assign_spinor_field_flt_cpu((spinor_field_flt *)s1, (spinor_field_flt *)s2),     \
        scalar_field *: add_assign_scalar_field_cpu((scalar_field *)s1, (scalar_field *)s2),                     \
        suNg_field *: add_assign_suNg_field_cpu((suNg_field *)s1, (suNg_field *)s2),                             \
        suNf_field *: add_assign_suNf_field_cpu((suNf_field *)s1, (suNf_field *)s2),                             \
        suNfc_field *: add_assign_suNfc_field_cpu((suNfc_field *)s1, (suNfc_field *)s2),                         \
        suNg_field_flt *: add_assign_suNg_field_flt_cpu((suNg_field_flt *)s1, (suNg_field_flt *)s2),             \
        suNf_field_flt *: add_assign_suNf_field_flt_cpu((suNf_field_flt *)s1, (suNf_field_flt *)s2),             \
        suNg_scalar_field *: add_assign_suNg_scalar_field_cpu((suNg_scalar_field *)s1, (suNg_scalar_field *)s2), \
        suNg_av_field *: add_assign_suNg_av_field_cpu((suNg_av_field *)s1, (suNg_av_field *)s2),                 \
        gtransf *: add_assign_gtransf_cpu((gtransf *)s1, (gtransf *)s2),                                         \
        clover_term *: add_assign_clover_term_cpu((clover_term *)s1, (clover_term *)s2),                         \
        clover_force *: add_assign_clover_force_cpu((clover_force *)s1, (clover_force *)s2),                     \
        staple_field *: add_assign_staple_field_cpu((staple_field *)s1, (staple_field *)s2))

#define prod_re_cpu(s1, s2)                                                                                   \
    _Generic((s2),                                                                                            \
        spinor_field *: prod_re_spinor_field_cpu((spinor_field *)s1, (spinor_field *)s2),                     \
        spinor_field_flt *: prod_re_spinor_field_flt_cpu((spinor_field_flt *)s1, (spinor_field_flt *)s2),     \
        scalar_field *: prod_re_scalar_field_cpu((scalar_field *)s1, (scalar_field *)s2),                     \
        suNg_field *: prod_re_suNg_field_cpu((suNg_field *)s1, (suNg_field *)s2),                             \
        suNf_field *: prod_re_suNf_field_cpu((suNf_field *)s1, (suNf_field *)s2),                             \
        suNfc_field *: prod_re_suNfc_field_cpu((suNfc_field *)s1, (suNfc_field *)s2),                         \
        suNg_field_flt *: prod_re_suNg_field_flt_cpu((suNg_field_flt *)s1, (suNg_field_flt *)s2),             \
        suNf_field_flt *: prod_re_suNf_field_flt_cpu((suNf_field_flt *)s1, (suNf_field_flt *)s2),             \
        suNg_scalar_field *: prod_re_suNg_scalar_field_cpu((suNg_scalar_field *)s1, (suNg_scalar_field *)s2), \
        suNg_av_field *: prod_re_suNg_av_field_cpu((suNg_av_field *)s1, (suNg_av_field *)s2),                 \
        gtransf *: prod_re_gtransf_cpu((gtransf *)s1, (gtransf *)s2),                                         \
        clover_term *: prod_re_clover_term_cpu((clover_term *)s1, (clover_term *)s2),                         \
        clover_force *: prod_re_clover_force_cpu((clover_force *)s1, (clover_force *)s2),                     \
        staple_field *: prod_re_staple_field_cpu((staple_field *)s1, (staple_field *)s2))

#define prod_im_cpu(s1, s2)                                                                                   \
    _Generic((s2),                                                                                            \
        spinor_field *: prod_im_spinor_field_cpu((spinor_field *)s1, (spinor_field *)s2),                     \
        spinor_field_flt *: prod_im_spinor_field_flt_cpu((spinor_field_flt *)s1, (spinor_field_flt *)s2),     \
        scalar_field *: prod_im_scalar_field_cpu((scalar_field *)s1, (scalar_field *)s2),                     \
        suNg_field *: prod_im_suNg_field_cpu((suNg_field *)s1, (suNg_field *)s2),                             \
        suNf_field *: prod_im_suNf_field_cpu((suNf_field *)s1, (suNf_field *)s2),                             \
        suNfc_field *: prod_im_suNfc_field_cpu((suNfc_field *)s1, (suNfc_field *)s2),                         \
        suNg_field_flt *: prod_im_suNg_field_flt_cpu((suNg_field_flt *)s1, (suNg_field_flt *)s2),             \
        suNf_field_flt *: prod_im_suNf_field_flt_cpu((suNf_field_flt *)s1, (suNf_field_flt *)s2),             \
        suNg_scalar_field *: prod_im_suNg_scalar_field_cpu((suNg_scalar_field *)s1, (suNg_scalar_field *)s2), \
        suNg_av_field *: prod_im_suNg_av_field_cpu((suNg_av_field *)s1, (suNg_av_field *)s2),                 \
        gtransf *: prod_im_gtransf_cpu((gtransf *)s1, (gtransf *)s2),                                         \
        clover_term *: prod_im_clover_term_cpu((clover_term *)s1, (clover_term *)s2),                         \
        clover_force *: prod_im_clover_force_cpu((clover_force *)s1, (clover_force *)s2),                     \
        staple_field *: prod_im_staple_field_cpu((staple_field *)s1, (staple_field *)s2))

#define prod_cpu(s1, s2)                                                                                   \
    _Generic((s2),                                                                                         \
        spinor_field *: prod_spinor_field_cpu((spinor_field *)s1, (spinor_field *)s2),                     \
        spinor_field_flt *: prod_spinor_field_flt_cpu((spinor_field_flt *)s1, (spinor_field_flt *)s2),     \
        scalar_field *: prod_scalar_field_cpu((scalar_field *)s1, (scalar_field *)s2),                     \
        suNg_field *: prod_suNg_field_cpu((suNg_field *)s1, (suNg_field *)s2),                             \
        suNf_field *: prod_suNf_field_cpu((suNf_field *)s1, (suNf_field *)s2),                             \
        suNfc_field *: prod_suNfc_field_cpu((suNfc_field *)s1, (suNfc_field *)s2),                         \
        suNg_field_flt *: prod_suNg_field_flt_cpu((suNg_field_flt *)s1, (suNg_field_flt *)s2),             \
        suNf_field_flt *: prod_suNf_field_flt_cpu((suNf_field_flt *)s1, (suNf_field_flt *)s2),             \
        suNg_scalar_field *: prod_suNg_scalar_field_cpu((suNg_scalar_field *)s1, (suNg_scalar_field *)s2), \
        suNg_av_field *: prod_suNg_av_field_cpu((suNg_av_field *)s1, (suNg_av_field *)s2),                 \
        gtransf *: prod_gtransf_cpu((gtransf *)s1, (gtransf *)s2),                                         \
        clover_term *: prod_clover_term_cpu((clover_term *)s1, (clover_term *)s2),                         \
        clover_force *: prod_clover_force_cpu((clover_force *)s1, (clover_force *)s2),                     \
        staple_field *: prod_staple_field_cpu((staple_field *)s1, (staple_field *)s2),                     \
        suNf_spinor *: spinor_prod_f((suNf_spinor *)s1, (suNf_spinor *)s2))

#define sqnorm_cpu(s1)                                                              \
    _Generic((s1),                                                                  \
        spinor_field *: sqnorm_spinor_field_cpu((spinor_field *)s1),                \
        spinor_field_flt *: sqnorm_spinor_field_flt_cpu((spinor_field_flt *)s1),    \
        scalar_field *: sqnorm_scalar_field_cpu((scalar_field *)s1),                \
        suNg_field *: sqnorm_suNg_field_cpu((suNg_field *)s1),                      \
        suNf_field *: sqnorm_suNf_field_cpu((suNf_field *)s1),                      \
        suNfc_field *: sqnorm_suNfc_field_cpu((suNfc_field *)s1),                   \
        suNg_field_flt *: sqnorm_suNg_field_flt_cpu((suNg_field_flt *)s1),          \
        suNf_field_flt *: sqnorm_suNf_field_flt_cpu((suNf_field_flt *)s1),          \
        suNg_scalar_field *: sqnorm_suNg_scalar_field_cpu((suNg_scalar_field *)s1), \
        suNg_av_field *: sqnorm_suNg_av_field_cpu((suNg_av_field *)s1),             \
        gtransf *: sqnorm_gtransf_cpu((gtransf *)s1),                               \
        clover_term *: sqnorm_clover_term_cpu((clover_term *)s1),                   \
        clover_force *: sqnorm_clover_force_cpu((clover_force *)s1),                \
        staple_field *: sqnorm_staple_field_cpu((staple_field *)s1))

#define max_cpu(s1)                                                              \
    _Generic((s1),                                                               \
        spinor_field *: max_spinor_field_cpu((spinor_field *)s1),                \
        spinor_field_flt *: max_spinor_field_flt_cpu((spinor_field_flt *)s1),    \
        scalar_field *: max_scalar_field_cpu((scalar_field *)s1),                \
        suNg_field *: max_suNg_field_cpu((suNg_field *)s1),                      \
        suNf_field *: max_suNf_field_cpu((suNf_field *)s1),                      \
        suNfc_field *: max_suNfc_field_cpu((suNfc_field *)s1),                   \
        suNg_field_flt *: max_suNg_field_flt_cpu((suNg_field_flt *)s1),          \
        suNf_field_flt *: max_suNf_field_flt_cpu((suNf_field_flt *)s1),          \
        suNg_scalar_field *: max_suNg_scalar_field_cpu((suNg_scalar_field *)s1), \
        suNg_av_field *: max_suNg_av_field_cpu((suNg_av_field *)s1),             \
        gtransf *: max_gtransf_cpu((gtransf *)s1),                               \
        clover_term *: max_clover_term_cpu((clover_term *)s1),                   \
        clover_force *: max_clover_force_cpu((clover_force *)s1),                \
        staple_field *: max_staple_field_cpu((staple_field *)s1),                \
        suNf_spinor *: spinor_max_f((suNf_spinor *)s1),                          \
        suNf_spinor_flt *: spinor_max_f_flt((suNf_spinor_flt *)s1))

#define zero_cpu(s1)                                                              \
    _Generic((s1),                                                                \
        spinor_field *: zero_spinor_field_cpu((spinor_field *)s1),                \
        spinor_field_flt *: zero_spinor_field_flt_cpu((spinor_field_flt *)s1),    \
        scalar_field *: zero_scalar_field_cpu((scalar_field *)s1),                \
        suNg_field *: zero_suNg_field_cpu((suNg_field *)s1),                      \
        suNf_field *: zero_suNf_field_cpu((suNf_field *)s1),                      \
        suNfc_field *: zero_suNfc_field_cpu((suNfc_field *)s1),                   \
        suNg_field_flt *: zero_suNg_field_flt_cpu((suNg_field_flt *)s1),          \
        suNf_field_flt *: zero_suNf_field_flt_cpu((suNf_field_flt *)s1),          \
        suNg_scalar_field *: zero_suNg_scalar_field_cpu((suNg_scalar_field *)s1), \
        suNg_av_field *: zero_suNg_av_field_cpu((suNg_av_field *)s1),             \
        gtransf *: zero_gtransf_cpu((gtransf *)s1),                               \
        clover_term *: zero_clover_term_cpu((clover_term *)s1),                   \
        clover_force *: zero_clover_force_cpu((clover_force *)s1),                \
        staple_field *: zero_staple_field_cpu((staple_field *)s1))

#define copy_cpu(s1, s2)                                                                                   \
    _Generic((s2),                                                                                         \
        spinor_field *: copy_spinor_field_cpu((spinor_field *)s1, (spinor_field *)s2),                     \
        spinor_field_flt *: copy_spinor_field_flt_cpu((spinor_field_flt *)s1, (spinor_field_flt *)s2),     \
        scalar_field *: copy_scalar_field_cpu((scalar_field *)s1, (scalar_field *)s2),                     \
        suNg_field *: copy_suNg_field_cpu((suNg_field *)s1, (suNg_field *)s2),                             \
        suNf_field *: copy_suNf_field_cpu((suNf_field *)s1, (suNf_field *)s2),                             \
        suNfc_field *: copy_suNfc_field_cpu((suNfc_field *)s1, (suNfc_field *)s2),                         \
        suNg_field_flt *: copy_suNg_field_flt_cpu((suNg_field_flt *)s1, (suNg_field_flt *)s2),             \
        suNf_field_flt *: copy_suNf_field_flt_cpu((suNf_field_flt *)s1, (suNf_field_flt *)s2),             \
        suNg_scalar_field *: copy_suNg_scalar_field_cpu((suNg_scalar_field *)s1, (suNg_scalar_field *)s2), \
        suNg_av_field *: copy_suNg_av_field_cpu((suNg_av_field *)s1, (suNg_av_field *)s2),                 \
        gtransf *: copy_gtransf_cpu((gtransf *)s1, (gtransf *)s2),                                         \
        clover_term *: copy_clover_term_cpu((clover_term *)s1, (clover_term *)s2),                         \
        clover_force *: copy_clover_force_cpu((clover_force *)s1, (clover_force *)s2),                     \
        staple_field *: copy_staple_field_cpu((staple_field *)s1, (staple_field *)s2))

#define g5_prod_re_cpu(s1, s2)                                                               \
    _Generic((s2),                                                                           \
        spinor_field *: g5_prod_re_spinor_field_cpu((spinor_field *)s1, (spinor_field *)s2), \
        spinor_field_flt *: g5_prod_re_spinor_field_flt_cpu((spinor_field_flt *)s1, (spinor_field_flt *)s2))

#define g5_prod_im_cpu(s1, s2)                                                               \
    _Generic((s2),                                                                           \
        spinor_field *: g5_prod_im_spinor_field_cpu((spinor_field *)s1, (spinor_field *)s2), \
        spinor_field_flt *: g5_prod_im_spinor_field_flt_cpu((spinor_field_flt *)s1, (spinor_field_flt *)s2))

#define g5_cpu(s1, s2)                                                               \
    _Generic((s2),                                                                   \
        spinor_field *: g5_spinor_field_cpu((spinor_field *)s1, (spinor_field *)s2), \
        spinor_field_flt *: g5_spinor_field_flt_cpu((spinor_field_flt *)s1, (spinor_field_flt *)s2))

#define g5_assign_cpu(s1)                                               \
    _Generic((s1),                                                      \
        spinor_field *: g5_assign_spinor_field_cpu((spinor_field *)s1), \
        spinor_field_flt *: g5_assign_spinor_field_flt_cpu((spinor_field_flt *)s1))

#define g0_cpu(s1, s2)                                                               \
    _Generic((s2),                                                                   \
        spinor_field *: g0_spinor_field_cpu((spinor_field *)s1, (spinor_field *)s2), \
        spinor_field_flt *: g0_spinor_field_flt_cpu((spinor_field_flt *)s1, (spinor_field_flt *)s2))

#define g1_cpu(s1, s2)                                                               \
    _Generic((s2),                                                                   \
        spinor_field *: g1_spinor_field_cpu((spinor_field *)s1, (spinor_field *)s2), \
        spinor_field_flt *: g1_spinor_field_flt_cpu((spinor_field_flt *)s1, (spinor_field_flt *)s2))

#define g2_cpu(s1, s2)                                                               \
    _Generic((s2),                                                                   \
        spinor_field *: g2_spinor_field_cpu((spinor_field *)s1, (spinor_field *)s2), \
        spinor_field_flt *: g2_spinor_field_flt_cpu((spinor_field_flt *)s1, (spinor_field_flt *)s2))

#define g3_cpu(s1, s2)                                                               \
    _Generic((s2),                                                                   \
        spinor_field *: g3_spinor_field_cpu((spinor_field *)s1, (spinor_field *)s2), \
        spinor_field_flt *: g3_spinor_field_flt_cpu((spinor_field_flt *)s1, (spinor_field_flt *)s2))

#define lc_cpu(r, k1, s1, k2, s2)                                                                               \
    _Generic((s2),                                                                                              \
        spinor_field *: lc_spinor_field_cpu((spinor_field *)r, k1, (spinor_field *)s1, k2, (spinor_field *)s2), \
        spinor_field_flt *: lc_spinor_field_flt_cpu((spinor_field_flt *)r, k1, (spinor_field_flt *)s1, k2,      \
                                                    (spinor_field_flt *)s2))

#define lc_add_assign_cpu(r, k1, s1, k2, s2)                                                                               \
    _Generic((s2),                                                                                                         \
        spinor_field *: lc_add_assign_spinor_field_cpu((spinor_field *)r, k1, (spinor_field *)s1, k2, (spinor_field *)s2), \
        spinor_field_flt *: lc_add_assign_spinor_field_flt_cpu((spinor_field_flt *)r, k1, (spinor_field_flt *)s1, k2,      \
                                                               (spinor_field_flt *)s2))

#define clc_cpu(r, k1, s1, k2, s2)                                                                               \
    _Generic((s2),                                                                                               \
        spinor_field *: clc_spinor_field_cpu((spinor_field *)r, k1, (spinor_field *)s1, k2, (spinor_field *)s2), \
        spinor_field_flt *: clc_spinor_field_flt_cpu((spinor_field_flt *)r, k1, (spinor_field_flt *)s1, k2,      \
                                                     (spinor_field_flt *)s2))

#define clc_add_assign_cpu(r, k1, s1, k2, s2)                                                                               \
    _Generic((s2),                                                                                                          \
        spinor_field *: clc_add_assign_spinor_field_cpu((spinor_field *)r, k1, (spinor_field *)s1, k2, (spinor_field *)s2), \
        spinor_field_flt *: clc_add_assign_spinor_field_flt_cpu((spinor_field_flt *)r, k1, (spinor_field_flt *)s1, k2,      \
                                                                (spinor_field_flt *)s2))

#endif

#endif