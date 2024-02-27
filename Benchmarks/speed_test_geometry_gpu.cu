#include "geometry.h"
#include "libhr_core.h"
#include <stdio.h>

__global__ void read_speedtest_kernel(suNf_spinor *field, int N, int block_start) {
    suNf_spinor s;
    int ix = blockIdx.x * blockDim.x + threadIdx.x;
    if (ix < N) {
        ix += block_start;
        read_gpu<double>(0, &s, field, ix, 0, 1);
    }
}

__global__ void write_speedtest_kernel(suNf_spinor *field, int N, int block_start) {
    suNf_spinor s;
    _spinor_zero_f(s);
    int ix = blockIdx.x * blockDim.x + threadIdx.x;
    if (ix < N) {
        ix += block_start;
        write_gpu<double>(0, &s, field, ix, 0, 1);
    }
}

void simple_read_gpu(spinor_field *field) {
    _PIECE_FOR(field->type, ixp) {
        const int N = field->type->master_end[ixp] - field->type->master_start[ixp] + 1;
        const int block_start = field->type->master_start[ixp];
        const int grid = (N - 1) / BLOCK_SIZE + 1;
        read_speedtest_kernel<<<grid, BLOCK_SIZE, 0, 0>>>(field->gpu_ptr, N, block_start);
    }
}

void simple_write_gpu(spinor_field *field) {
    _PIECE_FOR(field->type, ixp) {
        const int N = field->type->master_end[ixp] - field->type->master_start[ixp] + 1;
        const int block_start = field->type->master_start[ixp];
        const int grid = (N - 1) / BLOCK_SIZE + 1;
        write_speedtest_kernel<<<grid, BLOCK_SIZE, 0, 0>>>(field->gpu_ptr, N, block_start);
    }
}
