#ifdef WITH_GPU
#include "suN.h"
#include "suN_types.h"
#include "spinor_field.h"
#include "gpu.h"
#include "global.h"
#include "new_geometry.h"
#ifdef WITH_GPU
#include "global_gpu.h"

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
_DECLARE_KERNEL(scalar_field, suNg_vector, 1);
_DECLARE_KERNEL(avfield, suNg_algebra_vector, 4);
_DECLARE_KERNEL(gtransf, suNg, 1);
_DECLARE_KERNEL(clover_ldl, ldl_t, 1);
_DECLARE_KERNEL(clover_term, suNfc, 4);
_DECLARE_KERNEL(clover_force, suNf, 6);

#ifdef __cplusplus
    extern "C" {
#endif

#define _DECLARE_SYNC_BOX(_name, _type, _size) \
    void sync_box_to_buffer_gpu_##_name(geometry_descriptor *gd, \
                                    box_t *src, \
                                    _type *lattice, \
                                    void *sendbuf) \
    { \
        /* icoord array gives coordinates of inner lattice given a buffer index */ \
        coord4 *c = src->icoord; \
        coord4 *d_c; \
        /*TODO: This is incorrect, we need to calculate the volume right.*/ \
        int full_vol = boxVolume(src); \
        cudaMalloc((void**)&d_c, gd->nbuffers_gauge*full_vol*sizeof(coord4)); \
        cudaMemcpy(d_c, c, gd->nbuffers_gauge*full_vol*sizeof(coord4), cudaMemcpyHostToDevice); \
        \
        /* Iterate over pieces and fill out buffers */\
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
            int grid = (vol - 1)/BLOCK_SIZE +1; \
            stride = gd->master_end[ixp] - gd->master_start[ixp] + 1; \
            int block_start = gd->master_start[ixp]; \
            _type* lattice_block = lattice + (_size)*block_start; \
            box_to_buffer_kernel_##_name<<<grid, BLOCK_SIZE>>>(sendbuf, lattice_block, base_idx, stride, d_c, ipt_gpu, vol, block_start); \
        } \
    } 

#define _DECLARE_SYNC_FIELD(_name, _type, _geom) \
    void sync_field_to_buffer_gpu_##_name(geometry_descriptor *gd, \
                        _type *lattice,  \
                        void *sendbuf) \
    { \
        error(geometryBoxes==NULL, 1, __func__, "geometryBoxes are not initialized.\n"); \
        \
        /* Query first buffer box in the list */ \
        box_t *L = geometryBoxes->next; \
        \
        /* Iterate over all boxes*/ \
        /* The i-counter is necessary, because spinor-like and gauge-like fields have */ \
        /* different numbers of buffers */ \
        int i = 0; \
        while (L && i < gd->nbuffers_##_geom) \
        { \
            sync_box_to_buffer_gpu_##_name(gd, L->sendBox, lattice, sendbuf); \
            L=L->next; i++; \
        } \
    } 

#define _DECLARE_SYNC_FUNCTIONS(_name, _type, _size, _geom) \
    _DECLARE_SYNC_BOX(_name, _type, _size) \
    _DECLARE_SYNC_FIELD(_name, _type, _geom)

_DECLARE_SYNC_FUNCTIONS(spinor_field_f, suNf_spinor, 1, spinor);
_DECLARE_SYNC_FUNCTIONS(spinor_field_f_flt, suNf_spinor_flt, 1, spinor);
_DECLARE_SYNC_FUNCTIONS(sfield, double, 1, spinor);

_DECLARE_SYNC_FUNCTIONS(gfield, suNg, 4, gauge);
_DECLARE_SYNC_FUNCTIONS(gfield_flt, suNg_flt, 4, gauge);
_DECLARE_SYNC_FUNCTIONS(gfield_f, suNf, 4, gauge);
_DECLARE_SYNC_FUNCTIONS(gfield_f_flt, suNf_flt, 4, gauge);
_DECLARE_SYNC_FUNCTIONS(scalar_field, suNg_vector, 1, gauge);
_DECLARE_SYNC_FUNCTIONS(avfield, suNg_algebra_vector, 4, gauge);
_DECLARE_SYNC_FUNCTIONS(gtransf, suNg, 1, gauge);
_DECLARE_SYNC_FUNCTIONS(clover_ldl, ldl_t, 1, gauge);
_DECLARE_SYNC_FUNCTIONS(clover_term, suNfc, 4, gauge);
_DECLARE_SYNC_FUNCTIONS(clover_force, suNf, 6, gauge);

#ifdef __cplusplus
    }
#endif

#endif
#endif