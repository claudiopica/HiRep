
#include "update.h"
#include "libhr_core.h"
#include <math.h>
#include "memory.h"
#include "geometry.h"
#include "utils.h"

#ifdef WITH_SMEARING
#define gauge_ptr u_gauge_s->gpu_ptr
#else
#define gauge_ptr u_gauge->gpu_ptr
#endif

#ifdef WITH_GPU

__global__ void represent_gauge_field_kernel(suNf *gauge_f, suNg *gauge, int N, int block_start) {
    for (int id = blockIdx.x * blockDim.x + threadIdx.x; id < N; id += gridDim.x * blockDim.x) {
        suNf Ru;
        suNg u;
        int ix = id + block_start;
        for (int mu = 0; mu < 4; mu++) {
            read_gpu<double>(0, &u, gauge, ix, mu, 4);
            _group_represent2(&Ru, &u);
            write_gpu<double>(0, &Ru, gauge_f, ix, mu, 4);
        }
    }
}

void represent_gauge_field_gpu() {
    _PIECE_FOR(u_gauge_f->type, ixp) {
        const int N = u_gauge_f->type->master_end[ixp] - u_gauge_f->type->master_start[ixp] + 1;
        const int block_start = u_gauge_f->type->master_start[ixp];
        const int grid = (N - 1) / BLOCK_SIZE + 1;
        represent_gauge_field_kernel<<<grid, BLOCK_SIZE>>>(u_gauge_f->gpu_ptr, gauge_ptr, N, block_start);
    }

#ifdef DPHI_FLT
    assign_ud2u_f();
#endif

#if defined(WITH_CLOVER) || defined(WITH_EXPCLOVER)
    compute_clover_term();
#endif
}

#endif