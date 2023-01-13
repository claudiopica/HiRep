#include "geometry.h"
#include "libhr_core.h"

#include "./sync_to_buffer.hpp"

#define _DECLARE_SYNC_BOX(_name, _type, _size) \
    void sync_box_to_buffer_gpu_##_name(geometry_descriptor *gd, \
                                    box_t *src, \
                                    _type *lattice, \
                                    void *sendbuf) \
    { \
        /* TODO: we do not want to compare pointers */ \
        enum gd_type gd_t = GLOBAL; \
        if (gd == &glat_even) gd_t = EVEN; \
        if (gd == &glat_odd) gd_t = ODD; \
        \
        /* icoord array gives coordinates of inner lattice given a sendbuffer index */ \
        coord4 *c = src->icoord; \
        coord4 *d_c; \
        int full_vol = boxVolume(src); \
        cudaMalloc((void**)&d_c, glattice.nbuffers_gauge*full_vol*sizeof(coord4)); \
        cudaMemcpy(d_c, c, glattice.nbuffers_gauge*full_vol*sizeof(coord4), cudaMemcpyHostToDevice); \
        \
        if (gd_t & EVEN) { \
            const int ixp = 0; \
            const int buffer_stride = boxEvenVolume(src); \
            const int buffer_start = src->base_index; \
            const int grid = (buffer_stride - 1)/BLOCK_SIZE_SYNC + 1; \
            const int block_stride = gd->master_end[ixp] - gd->master_start[ixp] + 1; \
            const int block_start = gd->master_start[ixp]; \
            box_to_buffer_kernel_##_name<<<grid, BLOCK_SIZE_SYNC>>>((_type*)lattice, block_stride, block_start, gd->master_shift, \
                                                                    (_type*)sendbuf, buffer_stride, buffer_start, 0, d_c, ipt_gpu); \
        } \
        if (gd_t & ODD) { \
            const int ixp = (gd_t == GLOBAL) ? 1 : 0; \
            const int buffer_stride = boxOddVolume(src); \
            const int buffer_start = src->base_index_odd; \
            const int grid = (buffer_stride - 1)/BLOCK_SIZE_SYNC + 1; \
            const int block_stride = gd->master_end[ixp] - gd->master_start[ixp] + 1; \
            const int block_start = gd->master_start[ixp]; \
            box_to_buffer_kernel_##_name<<<grid, BLOCK_SIZE_SYNC>>>((_type*)lattice, block_stride, block_start, gd->master_shift, \
                                                                    (_type*)sendbuf, buffer_stride, buffer_start, 0, d_c, ipt_gpu); \
        }\
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