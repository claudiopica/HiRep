/***************************************************************************\
* Copyright (c) 2022, Sofie Martins                                         *   
* All rights reserved.                                                      * 
\***************************************************************************/

#define buffer_index_to_box_index \
        coord4 c = icoord[__idx_out_global]; \
        int __idx_in_global = ipt_ext_gpu(c.x[0], c.x[1], c.x[2], c.x[3]); 

#define _DECLARE_KERNEL(_name, _type, _size) \
    __global__ void box_to_buffer_kernel_##_name(kernel_field_input* input, coord4* icoord, int* ipt_gpu) \
    { \
        _KERNEL_FOR (input, _type, _size) \
        { \
            _find_index(buffer_index_to_box_index); \
            \
            _type site; \
            for (int comp = 0; comp < (_size); ++comp) { \
                _IN_FIELD_AT(site, _type, comp); \
                _WRITE_OUT_FIELD(site, _type, comp); \
            } \
        } \
    }

_DECLARE_KERNEL(spinor_field_f, suNf_spinor, 1);
_DECLARE_KERNEL(spinor_field_f_flt, suNf_spinor_flt, 1);
_DECLARE_KERNEL(sfield, double, 1);

_DECLARE_KERNEL(gfield, suNg, 4);
_DECLARE_KERNEL(gfield_f, suNf, 4);
_DECLARE_KERNEL(gfield_flt, suNg_flt, 4);
_DECLARE_KERNEL(gfield_f_flt, suNf_flt, 4);
_DECLARE_KERNEL(suNg_scalar_field, suNg_vector, 1);
_DECLARE_KERNEL(avfield, suNg_algebra_vector, 4);
_DECLARE_KERNEL(gtransf, suNg, 1);
_DECLARE_KERNEL(clover_ldl, ldl_t, 1);
_DECLARE_KERNEL(clover_term, suNfc, 4);
_DECLARE_KERNEL(clover_force, suNf, 6);
