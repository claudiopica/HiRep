/***************************************************************************\
* Copyright (c) 2022, Claudio Pica, Sofie Martins                           *   
* All rights reserved.                                                      * 
\***************************************************************************/
/**
 * @file gpu_geometry
 * @brief Implementation of geometry macros in (best possible) analogy to 
 *        the CPU geometry macros.
 */

#ifndef GPU_GEOMETRY_H
#define GPU_GEOMETRY_H

#include "new_geometry.h"
#include "Utils/generics.h"

#define _GPU_FIELD_BLK(s, i) (((s)->gpu_ptr) + ((s)->type->master_start[(i)] - (s)->type->master_shift))
#define _GPU_4FIELD_BLK(s, i) (((s)->gpu_ptr) + 4 * ((s)->type->master_start[(i)]))
#define _GPU_DFIELD_BLK(s, i, size) (((s)->gpu_ptr) + size * ((s)->type->master_start[(i)] - (s)->type->master_shift))

#define _BUF_GPU_FIELD_BLK(s, i) (((s)->gpu_ptr) + ((s)->type->rbuf_start[(i)] - (s)->type->master_shift))
#define _BUF_GPU_4FIELD_BLK(s, i) (((s)->gpu_ptr) + 4 * ((s)->type->rbuf_start[(i)] - (s)->type->master_shift))
#define _BUF_GPU_DFIELD_BLK(s, i, size) (((s)->gpu_ptr) + size * ((s)->type->rbuf_start[(i)] - (s)->type->master_shift))

// Kernel structure
#define _KERNEL_PIECE_FOR(_piece) for (int _piece = EVEN; _piece <= ODD; _piece++)

#define _IF_IN_BOX_IN(_input, _piece) \
    if (_input->gd_in & piece)        \
        if (blockIdx.x * BLOCK_SIZE + threadIdx.x < _input->vol_in[piece - 1])

#define _IF_IN_BOX_OUT(_input, _piece) \
    if (_input->gd_in & piece)         \
        if (blockIdx.x * BLOCK_SIZE + threadIdx.x < _input->vol_out[piece - 1])

typedef struct _kernel_field_input {
    void *field_in;
    size_t start_in_even;
    size_t start_in;
    size_t base_in[2];
    int stride_in;
    size_t vol_in[2];
    size_t master_shift_in;
    void *field_out;
    size_t stride_out;
    size_t base_out[2];
    size_t start_out;
    size_t vol_out[2];
    size_t master_shift_out;
    void *gauge;
    int *iup_gpu;
    int *idn_gpu;
    char *imask_gpu;
    enum gd_type gd_in;
} kernel_field_input;

#endif
#undef local_index