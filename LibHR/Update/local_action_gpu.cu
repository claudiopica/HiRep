#include "libhr_core.h"
#include "update.h"
#include "memory.h"
#include "utils.h"
#include "geometry.h"

#ifdef WITH_GPU

__global__ void _flip_sign(double *l, int N, int block_start) {
    for (int id = blockDim.x * blockIdx.x + threadIdx.x; id < N; id += gridDim.x * blockDim.x) {
        const int ix = id + block_start;
        double site;
        read_gpu<double>(0, &site, l, ix, 0, 1);
        site = -site;
        write_gpu<double>(0, &site, l, ix, 0, 1);
    }
}

__global__ void _loc_action_hmc(suNg_algebra_vector *momenta, suNg_vector *momenta_s, double *loc_action, int N,
                                int block_start) {
    for (int id = blockDim.x * blockIdx.x + threadIdx.x; id < N; id += gridDim.x * blockDim.x) {
        const int ix = id + block_start;
        double a = 0;
        double tmp;
        suNg_algebra_vector cmom;
        for (int j = 0; j < 4; ++j) {
            read_gpu<double>(0, &cmom, momenta, ix, j, 4);
            _algebra_vector_sqnorm_g(tmp, cmom);
            a += tmp;
        }
        double P2 = 0.0;
        if (momenta_s != NULL) {
            suNg_vector P;
            read_gpu<double>(0, &P, momenta_s, ix, 0, 1);
            _vector_prod_re_g(P2, P, P);
        }
        a *= 0.5 * _FUND_NORM2;
        double l;
        read_gpu<double>(0, &l, loc_action, ix, 0, 1);
        l += a + P2;
        write_gpu<double>(0, &l, loc_action, ix, 0, 1);
    }
}

__global__ void _pf_local_action(double *loc_action, suNf_spinor *pf, int N, int block_start) {
    for (int id = blockDim.x * blockIdx.x + threadIdx.x; id < N; id += gridDim.x * blockDim.x) {
        const int ix = id + block_start;
        double a = 0.;
        suNf_spinor s;
        double l;
        read_gpu<double>(0, &s, pf, ix, 0, 1);
        _spinor_prod_re_f(a, s, s);
        read_gpu<double>(0, &l, loc_action, ix, 0, 1);
        l += a;
        write_gpu<double>(0, &l, loc_action, ix, 0, 1);
    }
}

void local_hmc_action_gpu(local_action_type type, scalar_field *loc_action, suNg_av_field *momenta,
                          suNg_scalar_field *momenta_s) {
    _TWO_SPINORS_MATCHING(u_gauge, loc_action); /* check that action is defined on the global lattice */
    _TWO_SPINORS_MATCHING(loc_action, momenta);

    switch (type) {
    case NEW:
        cudaMemset(loc_action->gpu_ptr, 0, loc_action->type->gsize_spinor * sizeof(double));
        break;
    case DELTA:
        _PIECE_FOR(&glattice, ixp) {
            const int N = glattice.master_end[ixp] - glattice.master_start[ixp] + 1;
            const int block_start = glattice.master_start[ixp];
            const int grid = (N - 1) / BLOCK_SIZE + 1;
            _flip_sign<<<grid, BLOCK_SIZE>>>(loc_action->gpu_ptr, N, block_start);
        }
        break;
    default:
        error(1, 1, "local_hmc_action", "Invalid type");
    }

    _PIECE_FOR(&glattice, ixp) {
        const int N = glattice.master_end[ixp] - glattice.master_start[ixp] + 1;
        const int block_start = glattice.master_start[ixp];
        const int grid = (N - 1) / BLOCK_SIZE + 1;

        suNg_vector *scalar_momenta = NULL;
        if (u_scalar != NULL) { scalar_momenta = momenta_s->gpu_ptr; }
        _loc_action_hmc<<<grid, BLOCK_SIZE>>>(momenta->gpu_ptr, scalar_momenta, loc_action->gpu_ptr, N, block_start);
    }

    int nmon = num_mon();
    for (int i = 0; i < nmon; ++i) {
        monomial const *m = mon_n(i);
        m->add_local_action(m, loc_action);
    }
}

void pf_local_action_gpu(scalar_field *loc_action, spinor_field *pf) {
    _PIECE_FOR(pf->type, ixp) {
        const int N = glattice.master_end[ixp] - glattice.master_start[ixp] + 1;
        const int block_start = glattice.master_start[ixp];
        const int grid = (N - 1) / BLOCK_SIZE + 1;
        _pf_local_action<<<grid, BLOCK_SIZE>>>(loc_action->gpu_ptr, pf->gpu_ptr, N, block_start);
    }
}

void (*local_hmc_action)(local_action_type, scalar_field *, suNg_av_field *, suNg_scalar_field *) = local_hmc_action_gpu;
void (*pf_local_action)(scalar_field *, spinor_field *) = pf_local_action_gpu;

#endif
