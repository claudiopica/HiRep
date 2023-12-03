/***************************************************************************\
* Copyright (c) 2023, Sofie Martins                                         *
* All rights reserved.                                                      *
\***************************************************************************/

/**
 * @file random_generic.h
 * @brief Generic random functions
 */

#ifndef RANDOM_GENERIC_H
#define RANDOM_GENERIC_H

#include "random.h"
#include "utils.h"

#define random_field(s1)                                                            \
    _Generic((s1),                                                                  \
        spinor_field *: gaussian_spinor_field((spinor_field *)s1),                  \
        spinor_field_flt *: gaussian_spinor_field_flt((spinor_field_flt *)s1),      \
        scalar_field *: gaussian_scalar_field((scalar_field *)s1),                  \
        suNg_field *: random_u((suNg_field *)s1),                                   \
        suNf_field *: random_u_f((suNf_field *)s1),                                 \
        suNfc_field *: random_suNfc_field_cpu((suNfc_field *)s1),                   \
        suNg_field_flt *: random_gfield_flt_cpu((suNg_field_flt *)s1),              \
        suNf_field_flt *: random_gfield_f_flt_cpu((suNf_field_flt *)s1),            \
        suNg_scalar_field *: random_suNg_scalar_field_cpu((suNg_scalar_field *)s1), \
        suNg_av_field *: random_avfield_cpu((suNg_av_field *)s1),                   \
        gtransf *: random_gtransf_cpu((gtransf *)s1),                               \
        clover_term *: random_clover_term_cpu((clover_term *)s1),                   \
        clover_force *: random_clover_force_cpu((clover_force *)s1),                \
        staple_field *: random_staple_field_cpu((staple_field *)s1))

#endif