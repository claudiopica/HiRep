#define buffer_index_to_box_index \
        coord4 c = icoord[__idx_out_global]; \
        int __idx_in_global = ipt_ext_gpu(c.x[0], c.x[1], c.x[2], c.x[3]); 


#define _transfer_site(_type, _size) \
    _type __site; \
    for (int comp = 0; comp < (_size); ++comp) { \
        read_gpu_##_type (__stride_in,  __site, __in,  __idx_in_local,  comp); \
        write_gpu_##_type(__stride_out, __site, __out, __idx_out_local, comp); \
    }


