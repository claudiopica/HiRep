/***************************************************************************\
 * Copyright (c) 2024, Sofie Martins                                       *
 * All rights reserved.                                                    *
\***************************************************************************/

// Port of the file luscherweisz.c by Agostino Patella and Martin Hansen

#include "update.h"
#include "libhr_core.h"
#include "memory.h"
#include "io.h"
#include <string.h>
#include "utils.h"
#include "geometry.h"
#include "inverters.h"

#ifdef WITH_GPU

#define COMM (1 == 1)
#define NOCOMM (1 == 0)

static staple_field *stfld[8] = { NULL };
static suNg *stflt_gpu_ptr[8];
static suNg **stflt_gpu_ptr_d;

// These are from the linear algebra
template <class T> T global_sum_gpu(T *vector, int size);

__host__ __device__ double _lw_action_density(suNg **stfld, suNg *gauge, int ix, double beta, double c0, double c1) {
    double plaqs = 0;
    double rects = 0;
    double p;
    suNg w1;
    suNg u, s1, s2;

    for (int nu = 0; nu < 3; nu++) {
        for (int mu = nu + 1; mu < 4; mu++) {
            int i = mu - nu - 1;
            read_gpu<double>(0, &s1, stfld[2 * nu + 1], ix, i, 3);
            read_gpu<double>(0, &u, gauge, ix, mu, 4);
            _suNg_times_suNg_dagger(w1, s1, u);
            _suNg_trace_re(p, w1);
#ifdef PLAQ_WEIGHTS
            p *= plaq_weight[16 * ix + 4 * mu + nu];
#endif
            plaqs -= p;
        }
    }

    for (int nu = 0; nu < 4; nu++) {
        for (int i = 0; i < 3; i++) {
            int ixpnu = ix;
            read_gpu<double>(0, &s1, stfld[2 * nu + 1], ixpnu, i, 3);
            read_gpu<double>(0, &s2, stfld[2 * nu + 0], ixpnu, i, 3);
            _suNg_times_suNg_dagger(w1, s1, s2);
            _suNg_trace_re(p, w1);
#ifdef PLAQ_WEIGHTS
            int mu = (nu + i + 1) & 0x3;
            ixpnu = idn(ix, nu);
            p *= rect_weight[16 * ixpnu + 4 * mu + nu];
#endif
            rects -= p;
        }
    }

    return (beta / NG) * (c0 * plaqs + c1 * rects);
}

__global__ void _calculate_stfld(suNg **stfld, suNg *gauge, int N, int block_start, int *iup_gpu, int *idn_gpu) {
    suNg u, u1, s, wu1;
    _suNg_zero(s);
    int ix = blockIdx.x * blockDim.x + threadIdx.x;
    if (ix < N) {
        ix += block_start;
        for (int nu = 0; nu < 4; nu++) {
            int ixpnu = iup_gpu[4 * ix + nu];
            int ixmnu = idn_gpu[4 * ix + nu];

            for (int i = 0; i < 3; i++) {
                int mu = (nu + i + 1) & 0x3;
                int ixpmu = iup_gpu[4 * ix + mu];
                int ixpmumnu = idn_gpu[4 * ixpmu + nu];

                read_gpu<double>(0, &u, gauge, ixmnu, mu, 4);
                read_gpu<double>(0, &u1, gauge, ixpmumnu, nu, 4);
                _suNg_times_suNg(wu1, u, u1);
                read_gpu<double>(0, &u, gauge, ixmnu, nu, 4);
                _suNg_dagger_times_suNg(s, u, wu1);
                write_gpu<double>(0, &s, stfld[2 * nu + 0], ix, i, 3);

                read_gpu<double>(0, &u, gauge, ix, nu, 4);
                read_gpu<double>(0, &u1, gauge, ixpnu, mu, 4);
                _suNg_times_suNg(wu1, u, u1);
                read_gpu<double>(0, &u, gauge, ixpmu, nu, 4);
                _suNg_times_suNg_dagger(s, wu1, u);
                write_gpu<double>(0, &s, stfld[2 * nu + 1], ix, i, 3);
            }
        }
    }
}

__global__ void _lw_force(suNg **stfld, suNg *gauge, suNg_algebra_vector *force, int N, int block_start, double dt, double beta,
                          double c0, double c1, int *iup_gpu, int *idn_gpu) {
    suNg ws[4], wu1, wu2;
    suNg s;
    suNg u, u1, u0;
    suNg_algebra_vector f, wf1;

    int ix = blockIdx.x * blockDim.x + threadIdx.x;
    if (ix < N) {
        ix += block_start;

        for (int mu = 0; mu < 4; mu++) {
            _suNg_zero(ws[mu]);
        }

        for (int nu = 0; nu < 4; nu++) {
            for (int i = 0; i < 3; i++) {
                int mu = (nu + i + 1) & 0x3;
                read_gpu<double>(0, &s, stfld[2 * nu + 0], ix, i, 3);
                _suNg_add_assign(ws[mu], s);
                read_gpu<double>(0, &s, stfld[2 * nu + 1], ix, i, 3);
                _suNg_add_assign(ws[mu], s);
            }
        }

        for (int mu = 0; mu < 4; mu++) {
            read_gpu<double>(0, &u, gauge, ix, mu, 4);
            read_gpu<double>(0, &f, force, ix, mu, 4);
            _suNg_times_suNg_dagger(wu1, u, ws[mu]);
            _fund_algebra_project(wf1, wu1);
            _algebra_vector_mul_add_assign_g(f, dt * (-beta * c0 / NG), wf1);
            write_gpu<double>(0, &f, force, ix, mu, 4);
            _suNg_zero(ws[mu]);
        }

        for (int nu = 0; nu < 4; nu++) {
            int ixpnu = iup_gpu[4 * ix + nu];
            int ixmnu = idn_gpu[4 * ix + nu];

            for (int i = 0; i < 3; i++) {
                int mu = (nu + i + 1) & 0x3;
                int ixpmu = iup_gpu[4 * ix + mu];
                int ixpmunnu = idn_gpu[4 * ixpmu + nu];

                read_gpu<double>(0, &u, gauge, ixmnu, nu, 4);
                read_gpu<double>(0, &s, stfld[2 * nu + 0], ixmnu, i, 3);
                _suNg_dagger_times_suNg(wu1, u, s);
                read_gpu<double>(0, &u, gauge, ixpmunnu, nu, 4);
                _suNg_times_suNg(wu2, wu1, u);
                _suNg_add_assign(ws[mu], wu2);

                read_gpu<double>(0, &u, gauge, ix, nu, 4);
                read_gpu<double>(0, &s, stfld[2 * nu + 1], ixpnu, i, 3);
                _suNg_times_suNg(wu1, u, s);
                read_gpu<double>(0, &u, gauge, ixpmu, nu, 4);
                _suNg_times_suNg_dagger(wu2, wu1, u);
                _suNg_add_assign(ws[mu], wu2);
            }
        }

        for (int mu = 0; mu < 4; mu++) {
            int ixpmu = iup_gpu[4 * ix + mu];

            for (int i = 0; i < 3; i++) {
                int nu = (mu + i + 1) & 0x3;
                int ixpnu = iup_gpu[4 * ix + nu];
                int ixmnu = idn_gpu[4 * ix + nu];
                int ixmnupmu = iup_gpu[4 * ixmnu + mu];

                // one
                read_gpu<double>(0, &u, gauge, ixmnu, nu, 4);
                read_gpu<double>(0, &u1, gauge, ixmnu, mu, 4);
                _suNg_dagger_times_suNg(wu1, u, u1);
                read_gpu<double>(0, &s, stfld[2 * mu + 1], ixmnupmu, i, 3);
                _suNg_times_suNg(wu2, wu1, s);
                _suNg_add_assign(ws[mu], wu2);

                // two
                read_gpu<double>(0, &u, gauge, ix, nu, 4);
                read_gpu<double>(0, &u1, gauge, ixpnu, mu, 4);
                _suNg_times_suNg(wu1, u, u1);
                read_gpu<double>(0, &s, stfld[2 * mu + 1], ixpmu, i, 3);
                _suNg_times_suNg_dagger(wu2, wu1, s);
                _suNg_add_assign(ws[mu], wu2);

                // three
                read_gpu<double>(0, &u, gauge, ixmnu, mu, 4);
                read_gpu<double>(0, &s, stfld[2 * mu + 0], ixmnu, i, 3);
                _suNg_dagger_times_suNg(wu1, s, u);
                read_gpu<double>(0, &u, gauge, ixmnupmu, nu, 4);
                _suNg_times_suNg(wu2, wu1, u);
                _suNg_add_assign(ws[mu], wu2);

                // four
                read_gpu<double>(0, &s, stfld[2 * mu + 0], ix, i, 3);
                read_gpu<double>(0, &u, gauge, ixpnu, mu, 4);
                _suNg_times_suNg(wu1, s, u);
                read_gpu<double>(0, &u, gauge, ixpmu, nu, 4);
                _suNg_times_suNg_dagger(wu2, wu1, u);
                _suNg_add_assign(ws[mu], wu2);
            }
        }

        for (int mu = 0; mu < 4; mu++) {
            read_gpu<double>(0, &u, gauge, ix, mu, 4);
            _suNg_times_suNg_dagger(wu1, u, ws[mu]);
            _fund_algebra_project(wf1, wu1);
            read_gpu<double>(0, &f, force, ix, mu, 4);
            _algebra_vector_mul_add_assign_g(f, dt * (-beta * c1 / NG), wf1);
            write_gpu<double>(0, &f, force, ix, mu, 4);
        }
    }
}

__global__ void _lw_action(suNg **stfld, suNg *gauge, double beta, double c0, double c1, double *s, int N, int block_start) {
    for (int id = blockDim.x * blockIdx.x + threadIdx.x; id < N; id += blockDim.x * gridDim.x) {
        const int ix = id + block_start;
        s[id] = _lw_action_density(stfld, gauge, ix, beta, c0, c1);
    }
}

__global__ void _add_assign_loc_action(double *loc_action, double val, int iy) {
    if (int ix = blockDim.x * blockIdx.x + threadIdx.x < 1) {
        double l;
        read_gpu<double>(0, &l, loc_action, iy, 0, 1);
        l += val;
        write_gpu<double>(0, &l, loc_action, iy, 0, 1);
    }
}

void calculate_stfld_gpu(int comm) {
    if (stfld[0] == NULL) {
        for (int k = 0; k < 8; k++) {
            stfld[k] = alloc_staple_field(&glattice);
            //stfld[k]->comm_type = (comm_t)ALL_COMMS;
        }

        for (int k = 0; k < 8; k++) {
            stflt_gpu_ptr[k] = stfld[k]->gpu_ptr;
        }
        cudaMalloc((void **)&stflt_gpu_ptr_d, 8 * sizeof(suNg *));
        cudaMemcpy(stflt_gpu_ptr_d, stflt_gpu_ptr, 8 * sizeof(suNg *), cudaMemcpyHostToDevice);
    }

    for (int k = 0; k < 8; k++) {
        cudaMemset(stfld[k]->gpu_ptr, 0, glattice.gsize_gauge * sizeof(suNg) * 3);
    }

    _PIECE_FOR(&glattice, ixp) {
        const int N = glattice.master_end[ixp] - glattice.master_start[ixp] + 1;
        const int block_start = glattice.master_start[ixp];
        const int grid = (N - 1) / BLOCK_SIZE + 1;
        _calculate_stfld<<<grid, BLOCK_SIZE>>>(stflt_gpu_ptr_d, u_gauge->gpu_ptr, N, block_start, iup_gpu, idn_gpu);
    }

    if (comm) {
        for (int k = 0; k < 8; k++) {
            start_sendrecv_staple_field(stfld[k]);
            complete_sendrecv_staple_field(stfld[k]);
        }
    }
}

void lw_force_gpu(double dt, void *vpar) {
    force_gauge_par *par = (force_gauge_par *)vpar;
    suNg_av_field *force = *par->momenta;
    double beta = par->beta;
    double c0 = par->c0;
    double c1 = par->c1;

    start_sendrecv_suNg_field(u_gauge);
    complete_sendrecv_suNg_field(u_gauge);
    start_sendrecv_suNg_av_field(force);
    complete_sendrecv_suNg_av_field(force);

    calculate_stfld_gpu(COMM);

    for (int k = 0; k < 8; k++) {
        start_sendrecv_staple_field(stfld[k]);
        complete_sendrecv_staple_field(stfld[k]);
    }

    _PIECE_FOR(&glattice, ixp) {
        const int N = glattice.master_end[ixp] - glattice.master_start[ixp] + 1;
        const int block_start = glattice.master_start[ixp];
        const int grid = (N - 1) / BLOCK_SIZE + 1;
        _lw_force<<<grid, BLOCK_SIZE>>>(stflt_gpu_ptr_d, u_gauge->gpu_ptr, force->gpu_ptr, N, block_start, dt, beta, c0, c1,
                                        iup_gpu, idn_gpu);
    }

    apply_BCs_on_momentum_field(force);
}

double lw_action_gpu(double beta, double c0, double c1) {
    double *resPiece;
    double res = 0.0;
    calculate_stfld_gpu(NOCOMM);

    _PIECE_FOR(&glattice, ixp) {
        const int N = glattice.master_end[ixp] - glattice.master_start[ixp] + 1;
        resPiece = alloc_double_sum_field(N);
        const int block_start = glattice.master_start[ixp];
        const int grid = (N - 1) / BLOCK_SIZE + 1;
        _lw_action<<<grid, BLOCK_SIZE>>>(stflt_gpu_ptr_d, u_gauge->gpu_ptr, beta, c0, c1, resPiece, N, block_start);
        res += global_sum_gpu(resPiece, N);
    }

#ifdef WITH_MPI
    global_sum(&res, 1);
#endif
    return res;
}

void lw_local_action_gpu(scalar_field *loc_action, double beta, double c0, double c1) {
    double *resPiece;
    double res = 0.0;
    calculate_stfld_gpu(COMM);

    _PIECE_FOR(&glattice, ixp) {
        const int N = glattice.master_end[ixp] - glattice.master_start[ixp] + 1;
        const int block_start = glattice.master_start[ixp];
        const int grid = (N - 1) / BLOCK_SIZE + 1;
        resPiece = alloc_double_sum_field(N);
        _lw_action<<<grid, BLOCK_SIZE>>>(stflt_gpu_ptr_d, u_gauge->gpu_ptr, beta, c0, c1, resPiece, N, block_start);
        res += global_sum_gpu(resPiece, N);
    }

    int iy = ipt(2, 0, 0, 0);
    _add_assign_loc_action<<<1, 1>>>(loc_action->gpu_ptr, res, iy);
}

void (*calculate_stfld)(int comm) = calculate_stfld_gpu;
void (*lw_force)(double dt, void *vpar) = lw_force_gpu;
double (*lw_action)(double beta, double c0, double c1) = lw_action_gpu;
void (*lw_local_action)(scalar_field *loc_action, double beta, double c0, double c1) = lw_local_action_gpu;

#undef COMM
#undef NOCOMM

#endif
