
#include "suN.h"
#include "suN_types.h"
#include "spinor_field.h"
#include "gpu.h"
#include "global.h"
#include "new_geometry.h"
#ifdef WITH_GPU
#include "global_gpu.h"


#define _DECLARE_KERNEL(_name, _type) \
    __global__ void box_to_buffer_kernel_##_name(void *dst, _type *lattice, int base_index, int stride, coord4* icoord, int bytes_per_site, \
                                                int* ipt_gpu) \
    { \
        _type* src_uncast; \
        int dix = threadIdx.x*BLOCK_SIZE + blockIdx.x + base_index;\
        coord4 c = icoord[dix]; \
        int six = ipt_ext_gpu(c.x[0], c.x[1], c.x[2], c.x[3]); \
        read_gpu_##_type(stride, (*src_uncast), lattice, six, 0); \
        char *srcbuf = (char*)src_uncast;\
        char *dstbuf = ((char*)dst) + dix*bytes_per_site; \
        cudaMemcpyAsync(dstbuf, srcbuf, bytes_per_site, cudaMemcpyDeviceToDevice); \
    } \
    \
    __global__ void buffer_to_box_kernel_##_name(void *src, _type *lattice, int base_index, int stride, coord4* icoord, int bytes_per_site, \
                                            int* ipt_gpu) \
    { \
        int six = threadIdx.x*BLOCK_SIZE + blockIdx.x + base_index; \
        coord4 c = icoord[six]; \
        int dix = ipt_ext_gpu(c.x[0], c.x[1], c.x[2], c.x[3]); \
        \
        char* srcbuf = ((char*)src) + six * bytes_per_site; \
        _type* src_uncast = (_type*)srcbuf; \
        write_gpu_##_type(stride, (*src_uncast), lattice, dix, 0); \
    }

_DECLARE_KERNEL(gfield_f, suNf);
_DECLARE_KERNEL(spinor_field_f, suNf_spinor);

#ifdef __cplusplus
    extern "C" {
#endif

#define _DECLARE_SYNC_BOX(_name, _type, _size) \
    void sync_box_to_buffer_gpu_##_name(geometry_descriptor *gd, \
                                    box_t *src, \
                                    _type *lattice, \
                                    void *sendbuf) \
    { \
        const int bytes_per_site = (_size)*sizeof(*lattice); \
        _PIECE_FOR(gd, ixp) \
        { \
            int vol = 0; \
            int stride = 0; \
            if ((ixp%2)==0) vol = boxEvenVolume(src); \
            else vol = boxOddVolume(src); \
            coord4 *c = src->icoord; \
            coord4 *d_c; \
            /* These memory transfers might be a performance problem (SAM) */ \
            cudaMalloc((void**)&d_c, sizeof(c)); \
            cudaMemcpy(d_c, c, sizeof(c), cudaMemcpyHostToDevice); \
            int grid = (vol -1)/BLOCK_SIZE +1; \
            stride = gd->master_start[ixp] - gd->master_end[ixp] + 1; \
            box_to_buffer_kernel_##_name<<<grid, BLOCK_SIZE>>>(sendbuf, lattice, src->base_index, stride, d_c, bytes_per_site, ipt_gpu); \
        } \
    } \
    \
    void sync_buffer_to_box_gpu_##_name(geometry_descriptor *gd,  \
                                        box_t *src, \
                                        _type *lattice, \
                                        void *recvbuf) \
    { \
        const int bytes_per_site = (_size)*sizeof(*lattice); \
        _PIECE_FOR(gd, ixp) \
        { \
            int vol = 0; \
            int stride = 0; \
            if ((ixp % 2) == 0) vol = boxEvenVolume(src); \
            else vol = boxOddVolume(src); \
            coord4 *c = src->icoord; \
            coord4 *d_c; \
            cudaMalloc((void**)&d_c, sizeof(c)); \
            cudaMemcpy(d_c, c, sizeof(c), cudaMemcpyHostToDevice); \
            int grid = (vol - 1)/BLOCK_SIZE + 1; \
            stride  = gd->master_start[ixp] - gd->master_end[ixp] + 1;\
            buffer_to_box_kernel_##_name<<<grid, BLOCK_SIZE>>>(recvbuf, lattice, src->base_index, stride, d_c, bytes_per_site, ipt_gpu); \
        } \
    }

_DECLARE_SYNC_BOX(gfield_f, suNf, 4);
_DECLARE_SYNC_BOX(spinor_field_f, suNf_spinor, 1);

#ifdef __cplusplus
    }
#endif
#endif