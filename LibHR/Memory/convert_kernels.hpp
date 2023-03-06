/***************************************************************************\
* Copyright (c) 2023, Sofie Martins                                         *
* All rights reserved.                                                      *
\***************************************************************************/

#ifndef CONVERT_KERNELS_HPP
#define CONVERT_KERNELS_HPP

#ifdef WITH_GPU

#include "libhr_core.h"
#include "geometry.h"

template <int FIELD_DIM, typename REAL, typename SITE_TYPE>
__global__ void to_cpu_format_kernel(SITE_TYPE *out, SITE_TYPE *in, int N, int block_size, int master_shift) {
    SITE_TYPE *target;
    for (int ix = blockIdx.x * blockDim.x + threadIdx.x; ix < N; ix += blockDim.x * gridDim.x) {
        ix += block_size;
        for (int comp = 0; comp < FIELD_DIM; ++comp) {
            target = _DFIELD_AT_PTR(out, ix, comp, master_shift, FIELD_DIM);
            read_gpu<REAL>(0, target, in, (ix - master_shift), comp, FIELD_DIM);
        }
    }
}

template <int FIELD_DIM, typename REAL, typename SITE_TYPE>
__global__ void to_gpu_format_kernel(SITE_TYPE *out, SITE_TYPE *in, int N, int block_size, int master_shift) {
    SITE_TYPE *s;
    for (int ix = blockIdx.x * blockDim.x + threadIdx.x; ix < N; ix += blockDim.x * gridDim.x) {
        ix += block_size;
        for (int comp = 0; comp < FIELD_DIM; ++comp) {
            s = _DFIELD_AT_PTR(in, ix, comp, master_shift, FIELD_DIM);
            write_gpu<REAL>(0, s, out, (ix - master_shift), comp, FIELD_DIM);
        }
    }
}

template <int FIELD_DIM, typename REAL, typename FIELD_TYPE> void to_cpu_format_convert(FIELD_TYPE *out, FIELD_TYPE *in) {
    _PIECE_FOR(in->type, ixp) {
        const int N = in->type->master_end[ixp] - in->type->master_start[ixp] + 1;
        const int block_start = in->type->master_start[ixp];
        const int grid = (N - 1) / BLOCK_SIZE + 1;
        to_cpu_format_kernel<FIELD_DIM, REAL>
            <<<grid, BLOCK_SIZE>>>(out->gpu_ptr, in->gpu_ptr, N, block_start, in->type->master_shift);
    }
}

template <int FIELD_DIM, typename REAL, typename FIELD_TYPE> void to_gpu_format_convert(FIELD_TYPE *out, FIELD_TYPE *in) {
    _PIECE_FOR(in->type, ixp) {
        const int N = in->type->master_end[ixp] - in->type->master_start[ixp] + 1;
        const int block_start = in->type->master_start[ixp];
        const int grid = (N - 1) / BLOCK_SIZE + 1;
        to_gpu_format_kernel<FIELD_DIM, REAL>
            <<<grid, BLOCK_SIZE>>>(out->gpu_ptr, in->gpu_ptr, N, block_start, in->type->master_shift);
    }
}

#endif
#endif
