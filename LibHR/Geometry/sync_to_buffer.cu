/***************************************************************************\
* Copyright (c) 2022, Sofie Martins                                         *   
* All rights reserved.                                                      * 
\***************************************************************************/

#include "geometry.h"
#include "libhr_core.h"

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

#define buffer_index_to_box_index \
        coord4 c = icoord[__idx_out_global]; \
        int __idx_in_global = ipt_ext_gpu(c.x[0], c.x[1], c.x[2], c.x[3]); 

#define _DECLARE_KERNEL(_name, _type, _size) \
    __global__ void box_to_buffer_kernel_##_name(kernel_field_input* input, coord4* icoord, int* ipt_gpu) \
    { \
        _KERNEL_FOR (input, _type, _size) \
        { \
            _find_index(buffer_index_to_box_index); \
            \
            _type site; \
            for (int comp = 0; comp < (_size); ++comp) { \
                _IN_FIELD_AT(site, _type, comp); \
                _WRITE_OUT_FIELD(site, _type, comp); \
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