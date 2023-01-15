/***************************************************************************\
* Copyright (c) 2022, Sofie Martins                                         *   
* All rights reserved.                                                      * 
\***************************************************************************/

#ifdef WITH_GPU

#include "memory.h"
#include "gpu.h"

//TODO: inline comments
//TODO: doxygen docstrings


 /* Declare function to copy field from host to device */
    #define _DECLARE_COPY_TO(_name, _field_type, _site_type, _size, _is_spinorlike)                                 \
        void copy_to_gpu_##_name(_field_type *f)                                                                    \
        {                                                                                                           \
            _field_type *tmp = _ALLOCATE(_name);                                                                    \
            to_gpu_format_##_name(tmp, f);                                                                          \
            int field_size = _size * f->type->gsize_spinor * sizeof(*(f->gpu_ptr));                                 \
            cudaMemcpy(f->gpu_ptr, tmp->ptr, field_size, cudaMemcpyHostToDevice);                                   \
            free_##_name(tmp);                                                                                      \
        }

    /* Declare function to copy field from device to host */
    #define _DECLARE_COPY_FROM(_name, _field_type, _site_type, _size, _is_spinorlike)                               \
        void copy_from_gpu_##_name(_field_type *f)                                                                  \
        {                                                                                                           \
            _field_type *tmp = _ALLOCATE(_name);                                                                    \
            int field_size = _size * f->type->gsize_spinor * sizeof(*(f->gpu_ptr));                                 \
            cudaMemcpy(tmp->ptr, f->gpu_ptr, field_size, cudaMemcpyDeviceToHost);                                   \
            to_cpu_format_##_name(f, tmp);                                                                          \
            free_##_name(tmp);                                                                                      \
        }

#define _DECLARE_TRANSFER_FUNC(_name, _field_type, _site_type, _size, _geom)                                        \
    _DECLARE_COPY_TO(_name, _field_type, _site_type, _size, _geom)                                                  \
    _DECLARE_COPY_FROM(_name, _field_type, _site_type, _size, _geom)                                                \

#define _ALLOCATE(_name) alloc_##_name(1, f->type)
_DECLARE_TRANSFER_FUNC(spinor_field_f, spinor_field, suNf_spinor, 1, 1);
_DECLARE_TRANSFER_FUNC(spinor_field_f_flt, spinor_field_flt, suNf_spinor_flt, 1, 1);
_DECLARE_TRANSFER_FUNC(sfield, scalar_field, double, 1, 1);
#undef _ALLOCATE

#define _ALLOCATE(_name) alloc_##_name(f->type)
_DECLARE_TRANSFER_FUNC(gfield, suNg_field, suNg, 4, 0);
_DECLARE_TRANSFER_FUNC(gfield_flt, suNg_field_flt, suNg_flt, 4, 0);
_DECLARE_TRANSFER_FUNC(gfield_f, suNf_field, suNf, 4, 0);
_DECLARE_TRANSFER_FUNC(gfield_f_flt, suNf_field_flt, suNf_flt, 4, 0);
_DECLARE_TRANSFER_FUNC(suNg_scalar_field, suNg_scalar_field, suNg_vector, 1, 0);
_DECLARE_TRANSFER_FUNC(avfield, suNg_av_field, suNg_algebra_vector, 4, 0);
_DECLARE_TRANSFER_FUNC(gtransf, suNg_field, suNg, 1, 0);
_DECLARE_TRANSFER_FUNC(clover_ldl, ldl_field, ldl_t, 1, 0);
_DECLARE_TRANSFER_FUNC(clover_term, suNfc_field, suNfc, 4, 0);
_DECLARE_TRANSFER_FUNC(clover_force, suNf_field, suNf, 6, 0);
#undef _ALLOCATE

#undef _DECLARE_TRANSFER_FUNC
#undef _DECLARE_COPY_TO
#undef _DECLARE_COPY_FROM

#endif
