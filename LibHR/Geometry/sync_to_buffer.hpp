//TODO this code is chaotic and needs better integration with the concepts in the
//     new geometry.

#include "./kernel_structure.h"

#define _DECLARE_KERNEL(_name, _type, _size) \
    __global__ void box_to_buffer_kernel_##_name(_type* field_in, int stride_in, int start_in, int master_shift_in, \
                                                _type* field_out, int stride_out, int start_out, int master_shift_out, \
                                                coord4* icoord, int* ipt_gpu) \
    { \
        _KERNEL_FOR (field_in, stride_in, start_in, master_shift_in, \
                    field_out, stride_out, start_out, master_shift_out, \
                    _type, _size) \
        { \
            _find_index(buffer_index_to_box_index); \
            _transfer_site(_type, _size); \
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
