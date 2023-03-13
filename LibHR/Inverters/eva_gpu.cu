/***************************************************************************\
* Copyright (c) 2023, Sofie Martins                                         *
* All rights reserved.                                                      *
\***************************************************************************/

// Based on implementation in eva.c by Luigi Del Debbio, Claudio Pica and
// Agostino Patella

#ifdef WITH_GPU

#include "libhr_core.h"

__global__ void rotate_kernel(int n, suNf_spinor **pkk, hr_complex *v, int block_start, int block_size) {
    for (int id = blockDim.x * blockIdx.x + threadIdx.x; id < block_size * n; id += gridDim.x * blockDim.x) {
        suNf_spinor pk, pj;
        hr_complex *z;
        const int ix = id / n + block_start;
        const int k = id % n;

        read_gpu<double>(0, &pj, pkk[0], ix, 0, 1);
        z = &v[k];
        _spinor_mul_f(pk, *z, pj);
        for (int j = 1; j < n; j++) {
            read_gpu<double>(0, &pj, pkk[j], ix, 0, 1);
            z += n;
            _spinor_mul_add_assign_f(pk, *z, pj);
        }

        __syncthreads(); // All spinors read out, now we need to write
        write_gpu<double>(0, &pk, pkk[k], ix, 0, 1);
    }
}

extern "C" {

void rotate_gpu(int n, spinor_field *pkk, hr_complex v[]) {
    suNf_spinor **pkk_gpu = (suNf_spinor **)malloc(n * sizeof(suNf_spinor *));
    suNf_spinor **pkk_gpu_d;
    for (int i = 0; i < n; i++) {
        pkk_gpu[i] = (&pkk[i])->gpu_ptr;
    }

    cudaMalloc((void **)&pkk_gpu_d, n * sizeof(suNf_spinor *));
    cudaMemcpy(pkk_gpu_d, pkk_gpu, n * sizeof(suNf_spinor *), cudaMemcpyHostToDevice);

    _PIECE_FOR(pkk->type, ixp) {
        const size_t N = pkk->type->master_end[ixp] - pkk->type->master_start[ixp] + 1;
        const int start = pkk->type->master_start[ixp];
        const int grid = (N - 1) / BLOCK_SIZE + 1;
        rotate_kernel<<<grid, BLOCK_SIZE>>>(n, pkk_gpu_d, v, start, N);
    }

    cudaFree(pkk_gpu_d);
}
}

#endif
