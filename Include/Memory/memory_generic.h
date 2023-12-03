/***************************************************************************\
* Copyright (c) 2023, Sofie Martins                                         *
* All rights reserved.                                                      *
\***************************************************************************/

/**
 * @file memory_generic.h
 * @brief Generic memory allocation/deallocation
 */

#ifndef MEMORY_GENERIC_H
#define MEMORY_GENERIC_H

#include "random.h"
#include "utils.h"

#define alloc(_s1, _n, _geom)                                      \
    _Generic((_s1),                                                \
        spinor_field *: alloc_spinor_field(_n, _geom),             \
        spinor_field_flt *: alloc_spinor_field_flt(_n, _geom),     \
        scalar_field *: alloc_scalar_field(_n, _geom),             \
        suNg_field *: alloc_n_suNg_field(_n, _geom),               \
        suNf_field *: alloc_n_suNf_field(_n, _geom),               \
        suNfc_field *: alloc_n_suNfc_field(_n, _geom),             \
        suNg_field_flt *: alloc_n_suNg_field_flt(_n, _geom),       \
        suNf_field_flt *: alloc_n_suNf_field_flt(_n, _geom),       \
        suNg_scalar_field *: alloc_n_suNg_scalar_field(_n, _geom), \
        suNg_av_field *: alloc_n_suNg_av_field(_n, _geom),         \
        gtransf *: alloc_n_gtransf(_n, _geom),                     \
        clover_term *: alloc_n_clover_term(_n, _geom),             \
        clover_force *: alloc_n_clover_force(_n, _geom),           \
        staple_field *: alloc_n_staple_field(_n, _geom))

#define free_field(s1)                                                        \
    _Generic((s1),                                                            \
        spinor_field *: free_spinor_field((spinor_field *)s1),                \
        spinor_field_flt *: free_spinor_field_flt((spinor_field_flt *)s1),    \
        scalar_field *: free_scalar_field((scalar_field *)s1),                \
        suNg_field *: free_suNg_field((suNg_field *)s1),                      \
        suNf_field *: free_suNf_field((suNf_field *)s1),                      \
        suNfc_field *: free_suNfc_field((suNfc_field *)s1),                   \
        suNg_field_flt *: free_suNg_field_flt((suNg_field_flt *)s1),          \
        suNf_field_flt *: free_suNf_field_flt((suNf_field_flt *)s1),          \
        suNg_scalar_field *: free_suNg_scalar_field((suNg_scalar_field *)s1), \
        suNg_av_field *: free_suNg_av_field((suNg_av_field *)s1),             \
        gtransf *: free_gtransf((gtransf *)s1),                               \
        clover_term *: free_clover_term((clover_term *)s1),                   \
        clover_force *: free_clover_force((clover_force *)s1),                \
        staple_field *: free_staple_field((staple_field *)s1))

#ifdef WITH_GPU
#define copy_to_gpu(s1)                                                              \
    _Generic((s1),                                                                   \
        spinor_field *: copy_to_gpu_spinor_field((spinor_field *)s1),                \
        spinor_field_flt *: copy_to_gpu_spinor_field_flt((spinor_field_flt *)s1),    \
        scalar_field *: copy_to_gpu_scalar_field((scalar_field *)s1),                \
        suNg_field *: copy_to_gpu_suNg_field((suNg_field *)s1),                      \
        suNf_field *: copy_to_gpu_suNf_field((suNf_field *)s1),                      \
        suNfc_field *: copy_to_gpu_suNfc_field((suNfc_field *)s1),                   \
        suNg_field_flt *: copy_to_gpu_suNg_field_flt((suNg_field_flt *)s1),          \
        suNf_field_flt *: copy_to_gpu_suNf_field_flt((suNf_field_flt *)s1),          \
        suNg_scalar_field *: copy_to_gpu_suNg_scalar_field((suNg_scalar_field *)s1), \
        suNg_av_field *: copy_to_gpu_suNg_av_field((suNg_av_field *)s1),             \
        gtransf *: copy_to_gpu_gtransf((gtransf *)s1),                               \
        clover_term *: copy_to_gpu_clover_term((clover_term *)s1),                   \
        clover_force *: copy_to_gpu_clover_force((clover_force *)s1),                \
        staple_field *: copy_to_gpu_staple_field((staple_field *)s1))

#define copy_from_gpu(s1)                                                              \
    _Generic((s1),                                                                     \
        spinor_field *: copy_from_gpu_spinor_field((spinor_field *)s1),                \
        spinor_field_flt *: copy_from_gpu_spinor_field_flt((spinor_field_flt *)s1),    \
        scalar_field *: copy_from_gpu_scalar_field((scalar_field *)s1),                \
        suNg_field *: copy_from_gpu_suNg_field((suNg_field *)s1),                      \
        suNf_field *: copy_from_gpu_suNf_field((suNf_field *)s1),                      \
        suNfc_field *: copy_from_gpu_suNfc_field((suNfc_field *)s1),                   \
        suNg_field_flt *: copy_from_gpu_suNg_field_flt((suNg_field_flt *)s1),          \
        suNf_field_flt *: copy_from_gpu_suNf_field_flt((suNf_field_flt *)s1),          \
        suNg_scalar_field *: copy_from_gpu_suNg_scalar_field((suNg_scalar_field *)s1), \
        suNg_av_field *: copy_from_gpu_suNg_av_field((suNg_av_field *)s1),             \
        gtransf *: copy_from_gpu_gtransf((gtransf *)s1),                               \
        clover_term *: copy_from_gpu_clover_term((clover_term *)s1),                   \
        clover_force *: copy_from_gpu_clover_force((clover_force *)s1),                \
        staple_field *: copy_from_gpu_staple_field((staple_field *)s1))
#endif

#endif