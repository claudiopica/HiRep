/***************************************************************************\
* Copyright (c) 2022, Claudio Pica, Sofie Martins                           *   
* All rights reserved.                                                      * 
\***************************************************************************/
/**
 * @file gpu_geometry
 * @brief Implementation of geometry macros in (best possible) analogy to 
 *        the CPU geometry macros.
 */

#ifndef GPU_GEOMETRY_H
#define GPU_GEOMETRY_H

typedef struct _kernel_field_input {
    void* field_in;
    int stride_in;
    int start_in;
    int master_shift_in;
    void* field_out;
    int stride_out;
    int start_out;
    int master_shift_out;
} kernel_field_input;

#define _GPU_FIELD_BLK(s,i) (((s)->gpu_ptr) + ((s)->type->master_start[(i)] - (s)->type->master_shift))
#define _GPU_4FIELD_BLK(s,i) (((s)->gpu_ptr) + 4*((s)->type->master_start[(i)]))
#define _GPU_DFIELD_BLK(s,i,size) (((s)->gpu_ptr) + size*((s)->type->master_start[(i)] - (s)->type->master_shift))

#define _BUF_GPU_FIELD_BLK(s,i) (((s)->gpu_ptr) + ((s)->type->rbuf_start[(i)] - (s)->type->master_shift))
#define _BUF_GPU_4FIELD_BLK(s,i) (((s)->gpu_ptr) + 4*((s)->type->rbuf_start[(i)] - (s)->type->master_shift))
#define _BUF_GPU_DFIELD_BLK(s,i,size) (((s)->gpu_ptr) + size*((s)->type->rbuf_start[(i)] - (s)->type->master_shift))

#define _GPU_IDX_TO_LOCAL(_in, ix, ixp) ix - (_in)->type->master_start[(ixp)];

#define CONCAT(_name, _suffix) _name ## _suffix
#define _F_NAME(_name, _suffix) CONCAT(_name, _suffix)


// TODO: in spinor field -> read by comp, out spinor field -> read full spinor -> this needs to be clarified (SAM)
#define _IN_SPINOR_FIELD(_site, _comp) \
            read_gpu_suNf_vector(__stride_in, (_site), __in, __idx_in_local, _comp, 4);

#define _OUT_SPINOR_FIELD(_site) \
            read_gpu_suNf_spinor(__stride_out, (_site), __out, __idx_out_local, 0, 1);

#define _IN_GAUGE_FIELD(_site, _comp) \
            read_gpu_suNf(__stride_in, (_site), __in_gauge, __idx_in_local, _comp, 4);

#define _OUT_GAUGE_FIELD(_site, _comp) \
            read_gpu_suNf(__stride_out, (_site), __out_gauge, __idx_out_local, _comp, 4);

#define _WRITE_OUT_SPINOR_FIELD(_site) \
            write_gpu_suNf_spinor(__stride_out, (_site), __out, __idx_out_local, 0, 1);

#define _IN_FIELD_AT(_site, _type, _comp) \
            _F_NAME(read_gpu_,_type)(__stride_in, (_site), __in, __idx_in_local, _comp, 1); 

#define _OUT_FIELD_AT(_site, _comp) \
            _F_NAME(read_gpu_,_type)(__stride_out, (_site), __out, __idx_out_local, _comp, 1);

#define _IN_BLOCK_GAUGE_AT(_gauge_site, _comp) \
            _F_NAME(read_gpu_,_type)(__stride_in, (_site), __in_gauge, __idx_in_local, _comp, 4); 

#define _OUT_BLOCK_GAUGE_AT(_gauge_site, _comp) \
            _F_NAME(read_gpu_,_type)(__stride_in, (*_site), __out_gauge, __idx_out_local, _comp, 4); 

#define _WRITE_OUT_FIELD(_site, _type, _comp) \
            _F_NAME(write_gpu_,_type)(__stride_out, (_site), __out, __idx_out_local, _comp, 1);

#define _IN_DFIELD_AT(_site, _type, _comp, _dim) \
            _F_NAME(read_gpu_,_type)(__stride_in, (_site), __in, __idx_in_local, _comp, _dim);

#define _WRITE_OUT_DFIELD(_site, _type, _comp, _dim) \
            _F_NAME(write_gpu_,_type)(__stride_out, (_site), __out, __idx_out_local, _comp, _dim);

// TODO cast to void to avoid warnings when something is not needed (SAM)
#define _setup_striding(_input) \
    int __stride_in        = _input->stride_in; \
    int __start_in         = _input->start_in; \
    int __master_shift_in  = _input->master_shift_in; \
    int __stride_out       = _input->stride_out; \
    int __start_out        = _input->start_out; \
    int __master_shift_out = _input->master_shift_out; \
    int __idx_out_local    = blockIdx.x * BLOCK_SIZE + threadIdx.x; \
    int __idx_out_global   = __idx_out_local + __start_out; \

#define _find_index(__find_in_idx) \
    __find_in_idx; \
    int __idx_in_local   = __idx_in_global - __start_in;

#define _offset_fields(_input, _type, _size) \
    _type* __in  = _DFIELD_AT_PTR(((_type*)_input->field_in), __start_in, 0, __master_shift_in, (_size)); \
    _type* __out = _DFIELD_AT_PTR(((_type*)_input->field_out), __start_out, 0, __master_shift_out, (_size)); 

#define _offset_gauge_field(_gauge_type, _size, _gauge) \
    _gauge_type* __in_gauge = _DFIELD_AT_PTR((_gauge_type*)_gauge, __start_in, 0, 0, (_size)); \
    _gauge_type* __out_gauge = _DFIELD_AT_PTR((_gauge_type*)_gauge, __start_out, 0, 0, (_size));

#define _KERNEL_FOR(_input, _type, _size) \
    _setup_striding(_input) \
    _offset_fields(_input, _type, _size); \
    if (__idx_out_local < __stride_out)  

#define _KERNEL_FOR_WITH_GAUGE(_input, _gauge, _type, _gauge_type) \
    _setup_striding(_input) \
    _offset_fields(_input, _type, 1); \
    _offset_gauge_field(_gauge_type, 4, _gauge);  \
    if (__idx_out_local < __stride_out)


#undef _F_NAME
#undef CONCAT
#endif