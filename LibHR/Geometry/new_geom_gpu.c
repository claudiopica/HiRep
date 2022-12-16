
#include "suN.h"
#include "suN_types.h"
#include "spinor_field.h"
#include "gpu.h"
#include "global.h"
#include "new_geometry.h"
#ifdef WITH_GPU
#include "global_gpu.h"

//TODO the code in this file is broken

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


//__global__ void box_to_buffer_kernel_gfield_f(void *dst, suNf *lattice, int base_index, int stride, coord4* icoord, int bytes_per_site, \
                                                int* ipt_gpu, int vol) {}
_DECLARE_KERNEL(gfield_f, suNf, 4);
_DECLARE_KERNEL(spinor_field_f, suNf_spinor, 1);

#ifdef __cplusplus
    extern "C" {
#endif

#define _DECLARE_SYNC_BOX(_name, _type, _size) \
    void sync_box_to_buffer_gpu_##_name(geometry_descriptor *gd, \
                                    box_t *src, \
                                    _type *lattice, \
                                    void *sendbuf) \
    { \
        int full_vol = boxVolume(src); \
        coord4 *c = src->icoord; \
        coord4 *d_c; \
        /*TODO: This is incorrect, we need to calculate the volume right.*/ \
        cudaMalloc((void**)&d_c, 4*full_vol*sizeof(*c)); \
        cudaMemcpy(d_c, c, 4*full_vol*sizeof(*c), cudaMemcpyHostToDevice); \
        _PIECE_FOR(gd, ixp) \
        { \
            int vol = 0; \
            int stride = 0; \
            int base_idx = 0; \
            if ((ixp % 2)==0) { \
                vol = boxEvenVolume(src); \
                base_idx = src->base_index; \
            } else { \
                vol = boxOddVolume(src); \
                base_idx = src->base_index_odd; \
            } \
            if ((ixp%2)==0) vol = boxEvenVolume(src); \
            else vol = boxOddVolume(src); \
            int grid = (vol -1)/BLOCK_SIZE +1; \
            stride = gd->master_end[ixp] - gd->master_start[ixp] + 1; \
            int block_start = gd->master_start[ixp]; \
            _type* lattice_block = lattice + (_size)*block_start; \
            box_to_buffer_kernel_##_name<<<grid, BLOCK_SIZE>>>(sendbuf, lattice_block, base_idx, stride, d_c, ipt_gpu, vol, block_start); \
        } \
    } 

_DECLARE_SYNC_BOX(gfield_f, suNf, 4);
//void sync_box_to_buffer_gpu_gfield_f(geometry_descriptor *gd, box_t *src, suNf *lattice, void *sendbuf) {}

_DECLARE_SYNC_BOX(spinor_field_f, suNf_spinor, 1);

#ifdef __cplusplus
    }
#endif
#endif