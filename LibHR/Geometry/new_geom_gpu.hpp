//TODO this code is chaotic and needs better integration with the concepts in the
//     new geometry.

#define _DECLARE_KERNEL(_name, _type, _size) \
    __global__ void box_to_buffer_kernel_##_name(void *dst, _type *lattice_block, int base_index, int stride, coord4* icoord, \
                                                int* ipt_gpu, int vol, int block_start) \
    { \
        _type src_uncast; \
        int dix_loc = blockIdx.x*BLOCK_SIZE + threadIdx.x; \
        int dix = dix_loc + base_index; \
        if (dix_loc < vol) { \
            coord4 c = icoord[dix]; \
            int six = ipt_ext_gpu(c.x[0], c.x[1], c.x[2], c.x[3]); \
            int six_loc = six - block_start; \
            _type* dst_in = (_type*)dst; \
            _type* dst_uncast = _DFIELD_AT_PTR(dst_in, base_index, 0, 0, (_size));/* Add master shift here*/ \
            /*_type* dst_uncast = dst_in + (_size)*base_index;*/ \
            for (int comp = 0; comp < (_size); ++comp) { \
                read_gpu_##_type(stride, src_uncast, lattice_block, six_loc, comp); \
                write_gpu_##_type(vol, src_uncast, dst_uncast, dix_loc, comp); \
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
