#include "libhr_core.h"
#include "update.h"
#include "./staples_gpu.hpp"

__global__ void _force0_gpu(suNg *gauge, suNg_algebra_vector *force, double coeff, int *iup_gpu, int *idn_gpu,
                            double *plaq_weight, int N, int block_start) {
    for (int id = blockIdx.x * blockDim.x + threadIdx.x; id < N * 4; id += gridDim.x * blockDim.x) {
        const int ix = id / 4 + block_start;
        const int mu = id % 4;
        suNg s1, s2, u;
        suNg_algebra_vector f, f2;
        read_gpu<double>(0, &u, gauge, ix, mu, 4);
        staples_dev(ix, mu, &s1, gauge, iup_gpu, idn_gpu, plaq_weight);
        read_gpu<double>(0, &f2, force, ix, mu, 4);
        _suNg_times_suNg_dagger(s2, u, s1);
        _fund_algebra_project(f, s2);
        _algebra_vector_mul_add_assign_g(f2, coeff, f);
        write_gpu<double>(0, &f2, force, ix, mu, 4);
    }
}

void force0_kernel_gpu(suNg_av_field *force, double coeff) {
    _PIECE_FOR(&glattice, ixp) {
        const int N = glattice.master_end[ixp] - glattice.master_start[ixp] + 1;
        const int block_start = glattice.master_start[ixp];
        const int grid = (N * 4 - 1) / BLOCK_SIZE + 1;
        _force0_gpu<<<grid, BLOCK_SIZE, 0, 0>>>(u_gauge->gpu_ptr, force->gpu_ptr, coeff, iup_gpu, idn_gpu, plaq_weight_gpu, N,
                                                block_start);
        CudaCheckError();
    }
}