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

/* Spinor-like fields */
#define _GEOM_TYPE spinor

#define _FIELD_NAME spinor_field_f
#define _SITE_TYPE suNf_spinor
#define _FIELD_DIM 1
#include "TMPL/sync_to_buffer.cu.tmpl"

#define _FIELD_NAME spinor_field_f_flt
#define _SITE_TYPE suNf_spinor_flt
#define _FIELD_DIM 1
#include "TMPL/sync_to_buffer.cu.tmpl"

#define _FIELD_NAME sfield
#define _SITE_TYPE double
#define _FIELD_DIM 1
#include "TMPL/sync_to_buffer.cu.tmpl"

#undef _GEOM_TYPE

/* Gauge fields */
#define _GEOM_TYPE gauge

#define _FIELD_NAME gfield
#define _SITE_TYPE suNg
#define _FIELD_DIM 4
#include "TMPL/sync_to_buffer.cu.tmpl"

#define _FIELD_NAME gfield_flt
#define _SITE_TYPE suNg_flt
#define _FIELD_DIM 4
#include "TMPL/sync_to_buffer.cu.tmpl"

#define _FIELD_NAME gfield_f
#define _SITE_TYPE suNf
#define _FIELD_DIM 4
#include "TMPL/sync_to_buffer.cu.tmpl"

#define _FIELD_NAME gfield_f_flt
#define _SITE_TYPE suNf_flt
#define _FIELD_DIM 4
#include "TMPL/sync_to_buffer.cu.tmpl"

#define _FIELD_NAME suNg_scalar_field
#define _SITE_TYPE suNg_vector
#define _FIELD_DIM 1
#include "TMPL/sync_to_buffer.cu.tmpl"

#define _FIELD_NAME avfield
#define _SITE_TYPE suNg_algebra_vector
#define _FIELD_DIM 4
#include "TMPL/sync_to_buffer.cu.tmpl"

#define _FIELD_NAME gtransf
#define _SITE_TYPE suNg
#define _FIELD_DIM 1
#include "TMPL/sync_to_buffer.cu.tmpl"

#define _FIELD_NAME clover_ldl
#define _SITE_TYPE ldl_t
#define _FIELD_DIM 1
#include "TMPL/sync_to_buffer.cu.tmpl"

#define _FIELD_NAME clover_term
#define _SITE_TYPE suNfc
#define _FIELD_DIM 4
#include "TMPL/sync_to_buffer.cu.tmpl"

#define _FIELD_NAME clover_force
#define _SITE_TYPE suNf
#define _FIELD_DIM 6
#include "TMPL/sync_to_buffer.cu.tmpl"

#define _FIELD_NAME staple_field
#define _SITE_TYPE suNg
#define _FIELD_DIM 3
#include "TMPL/sync_to_buffer.cu.tmpl"

#undef _GEOM_TYPE