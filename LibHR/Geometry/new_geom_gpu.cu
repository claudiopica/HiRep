#include "geometry.h"
#include "libhr_core.h"

#include "./new_geom_gpu.hpp"

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
            int grid = (vol - 1)/BLOCK_SIZE +1; \
            stride = gd->master_end[ixp] - gd->master_start[ixp] + 1; \
            int block_start = gd->master_start[ixp]; \
            _type* lattice_block = lattice + (_size)*block_start; \
            box_to_buffer_kernel_##_name<<<grid, BLOCK_SIZE>>>(sendbuf, lattice_block, base_idx, stride, d_c, ipt_gpu, vol, block_start); \
        } \
    } 

#define _DECLARE_SYNC_FUNCTIONS(_name, _type, _size, _geom) \
    _DECLARE_SYNC_BOX(_name, _type, _size) 

_DECLARE_SYNC_FUNCTIONS(spinor_field_f, suNf_spinor, 1, spinor);
_DECLARE_SYNC_FUNCTIONS(spinor_field_f_flt, suNf_spinor_flt, 1, spinor);
_DECLARE_SYNC_FUNCTIONS(sfield, double, 1, spinor);

_DECLARE_SYNC_FUNCTIONS(gfield, suNg, 4, gauge);
_DECLARE_SYNC_FUNCTIONS(gfield_flt, suNg_flt, 4, gauge);
_DECLARE_SYNC_FUNCTIONS(gfield_f, suNf, 4, gauge);
_DECLARE_SYNC_FUNCTIONS(gfield_f_flt, suNf_flt, 4, gauge);
_DECLARE_SYNC_FUNCTIONS(suNg_scalar_field, suNg_vector, 1, gauge);
_DECLARE_SYNC_FUNCTIONS(avfield, suNg_algebra_vector, 4, gauge);
_DECLARE_SYNC_FUNCTIONS(gtransf, suNg, 1, gauge);
_DECLARE_SYNC_FUNCTIONS(clover_ldl, ldl_t, 1, gauge);
_DECLARE_SYNC_FUNCTIONS(clover_term, suNfc, 4, gauge);
_DECLARE_SYNC_FUNCTIONS(clover_force, suNf, 6, gauge);