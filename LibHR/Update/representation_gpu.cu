
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
    box_t *L = geometryBoxes;
    while (L) {
        int N = boxEvenVolume(L);
        int block_start = L->base_index;
        int grid = (N - 1) / BLOCK_SIZE + 1;
        represent_gauge_field_kernel<<<grid, BLOCK_SIZE, 0, 0>>>(u_gauge_f->gpu_ptr, gauge_ptr, N, block_start);

        N = boxOddVolume(L);
        block_start = L->base_index_odd;
        grid = (N - 1) / BLOCK_SIZE + 1;
        represent_gauge_field_kernel<<<grid, BLOCK_SIZE, 0, 0>>>(u_gauge_f->gpu_ptr, gauge_ptr, N, block_start);

        if (L->type == INNER) { complete_sendrecv_suNg_field(u_gauge); }

        L = L->next;
    }
}

#endif