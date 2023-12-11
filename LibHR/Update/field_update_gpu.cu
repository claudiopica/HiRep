#include "geometry.h"
#include "libhr_core.h"
#include "update.h"
#include "utils.h"

#define _PROJ_BIT (1 << 4)

__global__ void field_update_kernel(suNg *suNg_field, suNg_algebra_vector *force, int N, int block_start, double dt) {
    suNg_algebra_vector f;
    suNg u;
    for (int id = blockIdx.x * blockDim.x + threadIdx.x; id < N; id += gridDim.x * blockDim.x) {
        const int ix = id + block_start;
        for (int comp = 0; comp < 4; comp++) {
            read_gpu<double>(0, &u, suNg_field, ix, comp, 4);
            read_gpu<double>(0, &f, force, ix, comp, 4);
            ExpX(dt, &f, &u);
            write_gpu<double>(0, &u, suNg_field, ix, comp, 4);
            write_gpu<double>(0, &f, force, ix, comp, 4);
        }
    }
}

void exec_field_update(suNg_field *suNg_field, suNg_av_field *force, double dt) {
    _PIECE_FOR(&glattice, ixp) {
        const int N = glattice.master_end[ixp] - glattice.master_start[ixp] + 1;
        const int block_start = glattice.master_start[ixp];
        const int grid = (N - 1) / BLOCK_SIZE + 1;
        field_update_kernel<<<grid, BLOCK_SIZE>>>(suNg_field->gpu_ptr, force->gpu_ptr, N, block_start, dt);
    }
}
