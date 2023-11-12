#include "libhr_core.h"
#include "update.h"
#include "utils.h"

__global__ void project_kernel(suNg* gauge, int N, int block_start) {
    suNg u;
    for (int id = blockDim.x * blockIdx.x + threadIdx.x; id < N; id += gridDim.x * blockDim.x) {
      const int ix = id + block_start;
      for (int comp = 0; comp < 4; comp++) {
            read_gpu<double>(0, &u, gauge, ix, comp, 4);
            project_to_suNg(&u);
            write_gpu<double>(0, &u, gauge, ix, comp, 4);
      }
    }
}

void exec_project(void) {
    _PIECE_FOR(&glattice, ixp) {
        const int N = glattice.master_end[ixp] - glattice.master_start[ixp] + 1;
        const int block_start = glattice.master_start[ixp];
        const int grid = (N-1)/BLOCK_SIZE + 1;
        project_kernel<<<grid, BLOCK_SIZE>>>(u_gauge->gpu_ptr, N, block_start);
    }
}