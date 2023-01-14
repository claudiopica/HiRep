#include "geometry.h"
#include "libhr_core.h"

#include "./sync_to_buffer.hpp"

static kernel_field_input* get_even_input(
    void* lattice, 
    void* sendbuf, 
    box_t* srcbox) {
        kernel_field_input* input;
        input = (kernel_field_input*)malloc(sizeof(kernel_field_input));

        input->field_in = lattice;
        input->start_in = geometryBoxes->base_index;
        input->stride_in = boxEvenVolume(geometryBoxes);
        input->master_shift_in = 0;

        input->field_out = sendbuf;
        input->stride_out = boxEvenVolume(srcbox);
        input->start_out = srcbox->base_index;
        input->master_shift_out = 0;

        kernel_field_input* d_input;
        cudaError_t err;
        cudaMalloc((void**)&d_input, sizeof(kernel_field_input));
        cudaMemcpy(d_input, input, sizeof(kernel_field_input), cudaMemcpyHostToDevice);
        return d_input;
}

static kernel_field_input* get_odd_input(
    int master_shift,
    void* lattice, 
    void* sendbuf, 
    box_t* srcbox) {
        kernel_field_input* input;
        input = (kernel_field_input*)malloc(sizeof(kernel_field_input));

        input->field_in = lattice;
        input->start_in = geometryBoxes->base_index_odd;
        input->stride_in = boxOddVolume(geometryBoxes);
        input->master_shift_in = master_shift;

        input->field_out = sendbuf;
        input->stride_out = boxOddVolume(srcbox);
        input->start_out = srcbox->base_index_odd;
        input->master_shift_out = 0;

        kernel_field_input* d_input;
        cudaMalloc((void**)&d_input, sizeof(kernel_field_input));
        cudaMemcpy(d_input, input, sizeof(kernel_field_input), cudaMemcpyHostToDevice);
        return d_input;
}

#define _DECLARE_SYNC_BOX(_name) \
    void sync_box_to_buffer_gpu_##_name(geometry_descriptor *gd, \
                                    box_t *src, \
                                    void *lattice, \
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
            kernel_field_input* input = get_even_input(lattice, sendbuf, src); \
            const int grid = (boxEvenVolume(src) - 1)/BLOCK_SIZE_SYNC + 1; \
            box_to_buffer_kernel_##_name<<<grid, BLOCK_SIZE_SYNC>>>(input, d_c, ipt_gpu); \
        } \
        if (gd_t & ODD) { \
            kernel_field_input* input = get_odd_input(gd->master_shift, lattice, sendbuf, src); \
            const int grid = (boxOddVolume(src) - 1)/BLOCK_SIZE_SYNC + 1; \
            box_to_buffer_kernel_##_name<<<grid, BLOCK_SIZE_SYNC>>>(input, d_c, ipt_gpu); \
        }\
    } 

_DECLARE_SYNC_BOX(spinor_field_f);
_DECLARE_SYNC_BOX(spinor_field_f_flt);
_DECLARE_SYNC_BOX(sfield);

_DECLARE_SYNC_BOX(gfield);
_DECLARE_SYNC_BOX(gfield_flt);
_DECLARE_SYNC_BOX(gfield_f);
_DECLARE_SYNC_BOX(gfield_f_flt);
_DECLARE_SYNC_BOX(suNg_scalar_field);
_DECLARE_SYNC_BOX(avfield);
_DECLARE_SYNC_BOX(gtransf);
_DECLARE_SYNC_BOX(clover_ldl);
_DECLARE_SYNC_BOX(clover_term);
_DECLARE_SYNC_BOX(clover_force);