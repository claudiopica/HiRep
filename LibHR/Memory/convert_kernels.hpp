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
__global__ void to_cpu_format_kernel(SITE_TYPE *out, SITE_TYPE *in, size_t N, size_t block_size, size_t master_shift) {
    for (size_t id = blockIdx.x * blockDim.x + threadIdx.x; id < N; id += blockDim.x * gridDim.x) {
        SITE_TYPE *target;
        size_t ix = id + block_size;
        for (int comp = 0; comp < FIELD_DIM; ++comp) {
            target = _DFIELD_AT_PTR(out, ix, comp, master_shift, FIELD_DIM);
            read_gpu<REAL>(0, target, in, (ix - master_shift), comp, FIELD_DIM);
        }
    }
}

template <int FIELD_DIM, typename REAL, typename SITE_TYPE>
__global__ void to_gpu_format_kernel(SITE_TYPE *out, SITE_TYPE *in, size_t N, size_t block_size, size_t master_shift) {
    for (size_t id = blockIdx.x * blockDim.x + threadIdx.x; id < N; id += blockDim.x * gridDim.x) {
        SITE_TYPE *s;
        size_t ix = id + block_size;
        for (int comp = 0; comp < FIELD_DIM; ++comp) {
            s = _DFIELD_AT_PTR(in, ix, comp, master_shift, FIELD_DIM);
            write_gpu<REAL>(0, s, out, (ix - master_shift), comp, FIELD_DIM);
        }
    }
}

template <int FIELD_DIM, typename REAL, typename FIELD_TYPE> void to_cpu_format_convert(FIELD_TYPE *out, FIELD_TYPE *in) {
    _PIECE_FOR(in->type, ixp) {
        const size_t N = in->type->master_end[ixp] - in->type->master_start[ixp] + 1;
        const size_t block_start = in->type->master_start[ixp];
        const int grid = (N - 1) / BLOCK_SIZE + 1;
        to_cpu_format_kernel<FIELD_DIM, REAL>
            <<<grid, BLOCK_SIZE>>>(out->gpu_ptr, in->gpu_ptr, N, block_start, in->type->master_shift);
        CudaCheckError();
    }
}

template <int FIELD_DIM, typename REAL, typename FIELD_TYPE> void to_gpu_format_convert(FIELD_TYPE *out, FIELD_TYPE *in) {
    _PIECE_FOR(in->type, ixp) {
        const size_t N = in->type->master_end[ixp] - in->type->master_start[ixp] + 1;
        const size_t block_start = in->type->master_start[ixp];
        const int grid = (N - 1) / BLOCK_SIZE + 1;
        to_gpu_format_kernel<FIELD_DIM, REAL>
            <<<grid, BLOCK_SIZE>>>(out->gpu_ptr, in->gpu_ptr, N, block_start, in->type->master_shift);
        CudaCheckError();
    }
}

#endif
#endif
