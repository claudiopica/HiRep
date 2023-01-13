#define _IN_FIELD_AT(_site, _ix, _comp) \
            read_gpu_##_type(__stride_in, (*_site), __in, _ix, _comp); 

#define _OUT_FIELD_AT(_site, _ix, _comp) \
            read_gpu_##_type(__stride_out, (*site), __out, _ix, _comp);

#define _IN_BLOCK_GAUGE_AT(_gauge_site, _ix, _comp) \
            read_gpu_##_type(__stride_in, (*site), __out, _ix, _comp); 

#define _IN_BLOCK_GAUGE_AT(_gauge_site, _ix, _comp) \
            read_gpu_##_type(__stride_in, (*site), __out, _ix, _comp); 

#define _setup_striding(_stride_in, _start_in, _master_shift_in, \
                        _stride_out, _start_out, _master_shift_out) \
    int __stride_in        = _stride_in; \
    int __start_in         = _start_in; \
    int __master_shift_in  = _master_shift_in; \
    int __stride_out       = _stride_out; \
    int __start_out        = _start_out; \
    int __master_shift_out = _master_shift_out; \
    int __idx_out_local    = blockIdx.x * BLOCK_SIZE + threadIdx.x; \
    int __idx_out_global   = __idx_out_local + __start_out - __master_shift_out; \

#define _find_index(__find_in_idx) \
    __find_in_idx; \
    int __idx_in_local   = __idx_in_global - __start_in;

#define _offset_fields(_type, _size, _field_in, _field_out) \
    _type* __in  = _DFIELD_AT_PTR(((_type*)_field_in), __start_in, 0, __master_shift_in, (_size)); \
    _type* __out = _DFIELD_AT_PTR(((_type*)_field_out), __start_out, 0, __master_shift_out, (_size)); 

#define _offset_gauge_field(_gauge_type, _size, _gauge) \
    _gauge_type* __in_gauge = _DFIELD_AT_PTR((_gauge_type*)_gauge, __start_in, 0, __master_shift_in, (_size)); \
    _gauge_type* __out_gauge = _DFIELD_AT_PTR((_gauge_type*)_gauge, __start_out, 0, __master_shift_out, (_size));

#define _KERNEL_FOR(_field_in, _stride_in, _start_in, _master_shift_in, \
                    _field_out, _stride_out, _start_out, _master_shift_out, \
                    _type, _size) \
    _setup_striding(_stride_in, _start_in, _master_shift_in,  \
                    _stride_out, _start_out, _master_shift_out); \
    _offset_fields(_type, _size, _field_in, _field_out); \
    if (__idx_out_local < __stride_out)  

#define _KERNEL_FOR_WITH_GAUGE(_field_in, _stride_in, _start_in, _master_shift_in, \
                            _field_out, _stride_out, _start_out, _master_shift_out, \
                            _gauge, \
                            _type, _gauge_type, _size) \
        _setup_striding(_stride_in, _stride_in, _master_shift_in, \
                        _stride_out, _start_out, _master_shift_out); \
        _offset_fields(_type, _size, _field_in, _field_out); \
        _offset_gauge_field(_gauge_type, _size, _gauge); \
        if (__idx_out_local < __stride_out)
