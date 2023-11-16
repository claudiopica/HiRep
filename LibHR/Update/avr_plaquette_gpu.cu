/***************************************************************************\
* Copyright (c) 2023, Sofie Martins                                         *
* All rights reserved.                                                      *
\***************************************************************************/

/*******************************************************************************
*
* File plaquette.c
*
* Routines for the average plaquette
*
*******************************************************************************/

// With patterns from linear_algebra_gpu.cu + templates/kernels
// based on avr_plaquette.c

#include "libhr_core.h"
#include "update.h"
#include "inverters.h"

#ifdef WITH_GPU

static int *timeslices;
static int init = 0;

#ifdef PLAQ_WEIGHTS
#define PLAQ_WEIGHT_ARG , plaq_weight
#define PLAQ_WEIGHT_ARG_DEF , double *plaq_weight
#else
#define PLAQ_WEIGHT_ARG
#define PLAQ_WEIGHT_ARG_DEF
#endif

template <class T> T global_sum_gpu(T *vector, int size);

// Based on the implementation in avr_plaquette.c
// specify gauge as last argument if you want to
// evaluate this not on u_gauge
__device__ static double plaq_dev(int ix, int mu, int nu, suNg *gauge, int *iup_gpu PLAQ_WEIGHT_ARG_DEF) {
    int iy, iz;
    double p;
    suNg v1, v2, v3, v4, w1, w2, w3;

    iy = iup_gpu[4 * ix + mu];
    iz = iup_gpu[4 * ix + nu];

    read_gpu<double>(0, &v1, gauge, ix, mu, 4);
    read_gpu<double>(0, &v2, gauge, iy, nu, 4);
    read_gpu<double>(0, &v3, gauge, iz, mu, 4);
    read_gpu<double>(0, &v4, gauge, ix, nu, 4);

    _suNg_times_suNg(w1, v1, v2);
    _suNg_times_suNg(w2, v4, v3);
    _suNg_times_suNg_dagger(w3, w1, w2);

    _suNg_trace_re(p, w3);

#ifdef PLAQ_WEIGHTS
    if (plaq_weight == NULL) { return p; }
    return plaq_weight[ix * 16 + mu * 4 + nu] * p;
#else
    return p;
#endif
}

__device__ static void cplaq_dev(hr_complex *res, int ix, int mu, int nu, suNg *gauge, int *iup_gpu PLAQ_WEIGHT_ARG_DEF) {
    int iy, iz;
    suNg v1, v2, v3, v4, w1, w2, w3;
    double tmpre = 0.;
    double tmpim = 0.;

    iy = iup_gpu[4 * ix + mu];
    iz = iup_gpu[4 * ix + nu];

    read_gpu<double>(0, &v1, gauge, ix, mu, 4);
    read_gpu<double>(0, &v2, gauge, iy, nu, 4);
    read_gpu<double>(0, &v3, gauge, iz, mu, 4);
    read_gpu<double>(0, &v4, gauge, ix, nu, 4);

    _suNg_times_suNg(w1, v1, v2);
    _suNg_times_suNg(w2, v4, v3);
    _suNg_times_suNg_dagger(w3, w1, w2);

    _suNg_trace_re(tmpre, w3);

#ifndef GAUGE_SON
    _suNg_trace_im(tmpim, w3);
#endif

    double *t = (double *)res;
    t[0] = tmpre;
    t[1] = tmpim;

#ifdef PLAQ_WEIGHTS
    if (plaq_weight != NULL) {
        t[0] *= plaq_weight[ix * 16 + mu * 4 + nu];
        t[1] *= plaq_weight[ix * 16 + mu * 4 + nu];
    }
#endif
}

__device__ static double local_plaq_dev(int ix, suNg *gauge, int *iup_gpu PLAQ_WEIGHT_ARG_DEF) {
    double pa;

    pa = plaq_dev(ix, 1, 0, gauge, iup_gpu PLAQ_WEIGHT_ARG);
    pa += plaq_dev(ix, 2, 0, gauge, iup_gpu PLAQ_WEIGHT_ARG);
    pa += plaq_dev(ix, 2, 1, gauge, iup_gpu PLAQ_WEIGHT_ARG);
    pa += plaq_dev(ix, 3, 0, gauge, iup_gpu PLAQ_WEIGHT_ARG);
    pa += plaq_dev(ix, 3, 1, gauge, iup_gpu PLAQ_WEIGHT_ARG);
    pa += plaq_dev(ix, 3, 2, gauge, iup_gpu PLAQ_WEIGHT_ARG);

    return pa;
}

// TODO: put somewhere central for inverters
#define _CUDA_FOR(s, ixp, body)                                                        \
    do {                                                                               \
        _PIECE_FOR((s)->type, (ixp)) {                                                 \
            int N = (s)->type->master_end[(ixp)] - (s)->type->master_start[(ixp)] + 1; \
            unsigned int grid_size = (N - 1) / BLOCK_SIZE_LINEAR_ALGEBRA + 1;          \
            body;                                                                      \
            CudaCheckError();                                                          \
        }                                                                              \
    } while (0)

__global__ void _avr_plaquette(suNg *u, double *resField, int *iup_gpu, int N, int block_start PLAQ_WEIGHT_ARG_DEF) {
    for (int id = blockIdx.x * blockDim.x + threadIdx.x; id < N; id += blockDim.x * gridDim.x) {
        const int ix = id + block_start;
        resField[id] = plaq_dev(ix, 1, 0, u, iup_gpu PLAQ_WEIGHT_ARG);
        resField[id] += plaq_dev(ix, 2, 0, u, iup_gpu PLAQ_WEIGHT_ARG);
        resField[id] += plaq_dev(ix, 2, 1, u, iup_gpu PLAQ_WEIGHT_ARG);
        resField[id] += plaq_dev(ix, 3, 0, u, iup_gpu PLAQ_WEIGHT_ARG);
        resField[id] += plaq_dev(ix, 3, 1, u, iup_gpu PLAQ_WEIGHT_ARG);
        resField[id] += plaq_dev(ix, 3, 2, u, iup_gpu PLAQ_WEIGHT_ARG);
    }
}

double avr_plaquette_gpu() {
    double res = 0.0;
    double *resPiece;

    complete_sendrecv_gfield(u_gauge);

    _CUDA_FOR(u_gauge, ixp, resPiece = alloc_double_sum_field(N); (_avr_plaquette<<<grid_size, BLOCK_SIZE_LINEAR_ALGEBRA>>>(
        u_gauge->gpu_ptr, resPiece, iup_gpu, N, u_gauge->type->master_start[ixp] PLAQ_WEIGHT_ARG));
              res += global_sum_gpu(resPiece, N););

#ifdef WITH_MPI
    global_sum(&res, 1);
#endif

#ifdef BC_T_OPEN
    res /= 6.0 * NG * GLB_VOL3 * (GLB_T - 1);
#elif BC_T_SF
    res /= 6.0 * NG * GLB_VOL3 * (GLB_T - 2);
#else
    res /= 6.0 * NG * GLB_VOLUME;
#endif

    return res;
}

__global__ void _avr_plaquette_time(suNg *g, double *resPiece, int zero, int global_T, int *iup_gpu, int *timeslices, int N,
                                    int T, int block_start PLAQ_WEIGHT_ARG_DEF) {
    for (int id = blockIdx.x * blockDim.x + threadIdx.x; id < N; id += blockDim.x * gridDim.x) {
        const int ix = id + block_start;
        int nt = timeslices[ix];
        const int tc = (zero + nt + global_T) % global_T;

        resPiece[id + N * tc] += plaq_dev(ix, 1, 0, g, iup_gpu PLAQ_WEIGHT_ARG);
        resPiece[id + N * tc] += plaq_dev(ix, 2, 0, g, iup_gpu PLAQ_WEIGHT_ARG);
        resPiece[id + N * tc] += plaq_dev(ix, 3, 0, g, iup_gpu PLAQ_WEIGHT_ARG);

        resPiece[id + N * tc + N * global_T] += plaq_dev(ix, 2, 1, g, iup_gpu PLAQ_WEIGHT_ARG);
        resPiece[id + N * tc + N * global_T] += plaq_dev(ix, 3, 1, g, iup_gpu PLAQ_WEIGHT_ARG);
        resPiece[id + N * tc + N * global_T] += plaq_dev(ix, 3, 2, g, iup_gpu PLAQ_WEIGHT_ARG);
    }
}

void init_avr_plaquette_time() {
    int nts[GLB_VOLUME];

    for (int nt = 0; nt < T; nt++) {
        for (int nx = 0; nx < X; nx++) {
            for (int ny = 0; ny < Y; ny++) {
                for (int nz = 0; nz < Z; nz++) {
                    int ix = ipt(nt, nx, ny, nz);
                    nts[ix] = nt;
                }
            }
        }
    }

    cudaMalloc((void **)&timeslices, GLB_VOLUME * sizeof(int));
    cudaMemcpy(timeslices, nts, GLB_VOLUME * sizeof(int), cudaMemcpyHostToDevice);
}

void avr_plaquette_time_gpu(double *plaqt, double *plaqs) {
    double *resPiece;

    for (int nt = 0; nt < GLB_T; nt++) {
        plaqt[nt] = plaqs[nt] = 0.0;
    }

    if (!init) {
        init_avr_plaquette_time();
        init = 1;
    }

    start_sendrecv_gfield(u_gauge);
    complete_sendrecv_gfield(u_gauge);

    _CUDA_FOR(
        u_gauge, ixp, resPiece = alloc_double_sum_field(N * GLB_T * 2); cudaMemset(resPiece, 0, N * GLB_T * 2 * sizeof(double));
        (_avr_plaquette_time<<<grid_size, BLOCK_SIZE_LINEAR_ALGEBRA>>>(u_gauge->gpu_ptr, resPiece, zerocoord[0], GLB_T, iup_gpu,
                                                                       timeslices, N, T,
                                                                       u_gauge->type->master_start[ixp] PLAQ_WEIGHT_ARG));
        for (int nt = 0; nt < GLB_T; nt++) {
            int tc = (zerocoord[0] + nt + GLB_T) % GLB_T;
            plaqt[tc] += global_sum_gpu(resPiece + N * tc, N) / 3.0 / NG / GLB_VOL3;
        } for (int nt = 0; nt < GLB_T; nt++) {
            int tc = (zerocoord[0] + nt + GLB_T) % GLB_T;
            plaqs[tc] += global_sum_gpu(resPiece + N * tc + N * GLB_T, N) / 3.0 / NG / GLB_VOL3;
        });

    for (int nt = 0; nt < GLB_T; nt++) {
        global_sum(&plaqt[nt], 1);
        global_sum(&plaqs[nt], 1);
    }
}

__global__ void _full_plaquette(suNg *u, hr_complex *resPiece, int *iup_gpu, int N, int block_start PLAQ_WEIGHT_ARG_DEF) {
    for (int id = blockIdx.x * blockDim.x + threadIdx.x; id < N; id += blockDim.x * gridDim.x) {
        const int ix = id + block_start;
        hr_complex tmp;

        cplaq_dev(&tmp, ix, 1, 0, u, iup_gpu PLAQ_WEIGHT_ARG);
        resPiece[id] = tmp;

        cplaq_dev(&tmp, ix, 2, 0, u, iup_gpu PLAQ_WEIGHT_ARG);
        resPiece[id + N] = tmp;

        cplaq_dev(&tmp, ix, 2, 1, u, iup_gpu PLAQ_WEIGHT_ARG);
        resPiece[id + 2 * N] = tmp;

        cplaq_dev(&tmp, ix, 3, 0, u, iup_gpu PLAQ_WEIGHT_ARG);
        resPiece[id + 3 * N] = tmp;

        cplaq_dev(&tmp, ix, 3, 1, u, iup_gpu PLAQ_WEIGHT_ARG);
        resPiece[id + 4 * N] = tmp;

        cplaq_dev(&tmp, ix, 3, 2, u, iup_gpu PLAQ_WEIGHT_ARG);
        resPiece[id + 5 * N] = tmp;
    }
}

void full_plaquette_gpu(void) {
    start_sendrecv_gfield(u_gauge);
    complete_sendrecv_gfield(u_gauge);

    hr_complex pa[6];

    hr_complex r0 = 0.0;
    hr_complex r1 = 0.0;
    hr_complex r2 = 0.0;
    hr_complex r3 = 0.0;
    hr_complex r4 = 0.0;
    hr_complex r5 = 0.0;

    hr_complex *resPiece;

    _CUDA_FOR(u_gauge, ixp, resPiece = alloc_complex_sum_field(N * 6); cudaMemset(resPiece, 0, N * 6 * sizeof(hr_complex));
              (_full_plaquette<<<grid_size, BLOCK_SIZE_LINEAR_ALGEBRA>>>(u_gauge->gpu_ptr, resPiece, iup_gpu, N,
                                                                         u_gauge->type->master_start[ixp] PLAQ_WEIGHT_ARG));
              r0 += global_sum_gpu(resPiece, N); r1 += global_sum_gpu(resPiece + N, N);
              r2 += global_sum_gpu(resPiece + 2 * N, N); r3 += global_sum_gpu(resPiece + 3 * N, N);
              r4 += global_sum_gpu(resPiece + 4 * N, N); r5 += global_sum_gpu(resPiece + 5 * N, N););

    pa[0] = r0;
    pa[1] = r1;
    pa[2] = r2;
    pa[3] = r3;
    pa[4] = r4;
    pa[5] = r5;

    global_sum((double *)pa, 12);
    double *pad = (double *)pa;
    for (int k = 0; k < 12; k++) {
#ifdef BC_T_OPEN
        pad[k] /= NG * GLB_VOLUME * (GLB_T - 1) / GLB_T;
#else
        pad[k] /= NG * GLB_VOLUME;
#endif
    }

    lprintf("PLAQ", 0, "Plaq(%d,%d) = ( %f , %f )\n", 1, 0, creal(pa[0]), cimag(pa[0]));
    lprintf("PLAQ", 0, "Plaq(%d,%d) = ( %f , %f )\n", 2, 0, creal(pa[1]), cimag(pa[1]));
    lprintf("PLAQ", 0, "Plaq(%d,%d) = ( %f , %f )\n", 2, 1, creal(pa[2]), cimag(pa[2]));
    lprintf("PLAQ", 0, "Plaq(%d,%d) = ( %f , %f )\n", 3, 0, creal(pa[3]), cimag(pa[3]));
    lprintf("PLAQ", 0, "Plaq(%d,%d) = ( %f , %f )\n", 3, 1, creal(pa[4]), cimag(pa[4]));
    lprintf("PLAQ", 0, "Plaq(%d,%d) = ( %f , %f )\n", 3, 2, creal(pa[5]), cimag(pa[5]));
}

double (*avr_plaquette)(void) = avr_plaquette_gpu;
void (*full_plaquette)(void) = full_plaquette_gpu;
void (*avr_plaquette_time)(double *plaqt, double *plaqs) = avr_plaquette_time_gpu;
#endif

#undef PLAQ_WEIGHT_ARG
#undef PLAQ_WEIGHT_ARG_DEF
