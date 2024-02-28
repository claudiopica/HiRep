/***************************************************************************
 * Copyright (c) 2023, Sofie Martins                                        *
 * All rights reserved.                                                     *
 ***************************************************************************/

#include "libhr_core.h"
#include "geometry.h"
#include "utils.h"
#include "memory.h"
#include "io.h"

#ifdef WITH_GPU
#if defined(WITH_CLOVER) || defined(WITH_EXPCLOVER)
static double sigma;
static double csw_value;

#define c_idx(i, j) ((i) * ((i) + 1) / 2 + (j))

#define _CUDA_CALL(_type, _grid, _N, _block_start, _ixp, _call)                  \
    _PIECE_FOR(_type, _ixp) {                                                    \
        const int N = (_type)->master_end[ixp] - (_type)->master_start[ixp] + 1; \
        const int grid = (N - 1) / BLOCK_SIZE + 1;                               \
        const int block_start = (_type)->master_start[ixp];                      \
        _call;                                                                   \
        CudaCheckError();                                                        \
    }

__device__ static void clover_loop(int ix, int mu, int nu, suNf *u, suNf *gauge, int *iup_gpu, int *idn_gpu) {
    int o1, o2, o3;
    suNf s1, s2, s3;
    suNf u1, u2;

    // Leaf 1
    o1 = iup_gpu[4 * ix + mu];
    o2 = iup_gpu[4 * ix + nu];
    read_gpu<double>(0, &u1, gauge, ix, mu, 4);
    read_gpu<double>(0, &u2, gauge, o1, nu, 4);
    _suNf_times_suNf(s1, u1, u2);
    read_gpu<double>(0, &u1, gauge, ix, nu, 4);
    read_gpu<double>(0, &u2, gauge, o2, mu, 4);
    _suNf_times_suNf(s2, u1, u2);
    _suNf_times_suNf_dagger(*u, s1, s2);

    // Leaf 2
    o1 = idn_gpu[4 * ix + mu];
    o2 = iup_gpu[4 * o1 + nu];

    read_gpu<double>(0, &u1, gauge, o1, nu, 4);
    read_gpu<double>(0, &u2, gauge, o2, mu, 4);
    _suNf_times_suNf(s1, u1, u2);
    read_gpu<double>(0, &u1, gauge, ix, nu, 4);
    _suNf_times_suNf_dagger(s2, u1, s1);
    read_gpu<double>(0, &u1, gauge, o1, mu, 4);
    _suNf_times_suNf(s3, s2, u1);
    _suNf_add_assign(*u, s3);

    // Leaf 3
    o1 = idn_gpu[4 * ix + mu];
    o2 = idn_gpu[4 * ix + nu];
    o3 = idn_gpu[4 * o1 + nu];
    read_gpu<double>(0, &u1, gauge, o3, nu, 4);
    read_gpu<double>(0, &u2, gauge, o1, mu, 4);
    _suNf_times_suNf(s1, u1, u2);
    read_gpu<double>(0, &u1, gauge, o3, mu, 4);
    _suNf_dagger_times_suNf(s2, s1, u1);
    read_gpu<double>(0, &u1, gauge, o2, nu, 4);
    _suNf_times_suNf(s3, s2, u1);
    _suNf_add_assign(*u, s3);

    // Leaf 4
    o1 = idn_gpu[4 * ix + nu];
    o2 = iup_gpu[4 * o1 + mu];
    read_gpu<double>(0, &u1, gauge, o1, nu, 4);
    read_gpu<double>(0, &u2, gauge, o1, mu, 4);
    _suNf_dagger_times_suNf(s1, u1, u2);
    read_gpu<double>(0, &u1, gauge, o2, nu, 4);
    read_gpu<double>(0, &u2, gauge, ix, mu, 4);
    _suNf_times_suNf_dagger(s2, u1, u2);
    _suNf_times_suNf(s3, s1, s2);
    _suNf_add_assign(*u, s3);
}

__global__ static void _ldl(ldl_t *cl_ldl_gpu, int dir, int N, int block_start) {
    for (int id = blockIdx.x * blockDim.x + threadIdx.x; id < N; id += gridDim.x * blockDim.x) {
        const int ix = id + block_start;
        for (int i = 0; i < 2 * NF; i++) {
            hr_complex c1;
            read_clover(&c1, cl_ldl_gpu, ix, dir, c_idx(i, i));

            for (int k = 0; k < i; k++) {
                hr_complex c2, c3;
                read_clover(&c2, cl_ldl_gpu, ix, dir, c_idx(i, k));
                read_clover(&c3, cl_ldl_gpu, ix, dir, c_idx(k, k));
                c1 -= _complex_prod(c2, c2) * creal(c3);
            }
            write_clover(&c1, cl_ldl_gpu, ix, dir, c_idx(i, i));
            __syncthreads();

            for (int j = i + 1; j < 2 * NF; j++) {
                hr_complex c5;
                read_clover(&c5, cl_ldl_gpu, ix, dir, c_idx(j, i));
                for (int k = 0; k < i; k++) {
                    hr_complex c2, c3, c4;
                    read_clover(&c2, cl_ldl_gpu, ix, dir, c_idx(i, k));
                    read_clover(&c3, cl_ldl_gpu, ix, dir, c_idx(j, k));
                    read_clover(&c4, cl_ldl_gpu, ix, dir, c_idx(k, k));
                    c5 -= _complex_prod(c2, c3) * creal(c4);
                }
                _complex_div(c5, c5, c1);
                write_clover(&c5, cl_ldl_gpu, ix, dir, c_idx(j, i));
            }
        }
    }
}

__global__ static void _compute_ldl_decomp(suNfc *clover_term_gpu, ldl_t *cl_ldl_gpu, double sigma, int N, int block_start) {
    for (int id = blockDim.x * blockIdx.x + threadIdx.x; id < N * NF * NF; id += gridDim.x * blockDim.x) {
        const int ix = id / NF / NF + block_start;
        const int i = (id % (NF * NF)) / NF;
        const int j = id % NF;

        hr_complex c;

        const int m = i + NF;
        const int n = j + NF;

        read_clover_term_comp(&c, clover_term_gpu, ix, 1, j, i);
        c = conj(c);
        write_clover(&c, cl_ldl_gpu, ix, UP, c_idx(m, j));

        read_clover_term_comp(&c, clover_term_gpu, ix, 3, j, i);
        c = conj(c);
        write_clover(&c, cl_ldl_gpu, ix, DOWN, c_idx(m, j));

        if (i >= j) {
            read_clover_term_comp(&c, clover_term_gpu, ix, 0, i, j);
            if (i == j) { _complex_add_assign(c, sigma); }
            write_clover(&c, cl_ldl_gpu, ix, UP, c_idx(i, j));

            read_clover_term_comp(&c, clover_term_gpu, ix, 0, i, j);
            _complex_minus(c, c);
            if (i == j) { _complex_add_assign(c, sigma); }
            write_clover(&c, cl_ldl_gpu, ix, UP, c_idx(m, n));

            read_clover_term_comp(&c, clover_term_gpu, ix, 2, i, j);
            if (i == j) { _complex_add_assign(c, sigma); }
            write_clover(&c, cl_ldl_gpu, ix, DOWN, c_idx(i, j));

            read_clover_term_comp(&c, clover_term_gpu, ix, 2, i, j);
            _complex_minus(c, c);
            if (i == j) { _complex_add_assign(c, sigma); }
            write_clover(&c, cl_ldl_gpu, ix, DOWN, c_idx(m, n));
        }
    }
}

__global__ static void _compute_clover_force(ldl_t *cl_ldl_gpu, suNf *cl_force_gpu, double coeff, int N, int block_start) {
    for (int id = blockDim.x * blockIdx.x + threadIdx.x; id < N; id += gridDim.x * blockDim.x) {
        int ix = id + block_start;
        hr_complex A[2 * NF][2 * NF];
        hr_complex B[2 * NF][2 * NF];
        memset(A, 0, sizeof(A));
        memset(B, 0, sizeof(B));

        hr_complex a11;
        hr_complex a12;
        hr_complex a21;
        hr_complex a22;
        hr_complex a33;
        hr_complex a34;
        hr_complex a43;
        hr_complex a44;

        hr_complex U, L, c;

        for (int n = 0; n < 2 * NF; n++) {
            A[n][n] = coeff;
            B[n][n] = coeff;

            for (int i = n; i < 2 * NF; i++) {
                for (int k = 0; k < i; k++) {
                    read_clover(&U, cl_ldl_gpu, ix, UP, c_idx(i, k));
                    A[i][n] -= U * A[k][n];
                    read_clover(&L, cl_ldl_gpu, ix, DOWN, c_idx(i, k));
                    B[i][n] -= L * B[k][n];
                }
            }
            for (int i = 2 * NF - 1; i >= n; i--) {
                read_clover(&U, cl_ldl_gpu, ix, UP, c_idx(i, i));
                A[i][n] /= creal(U);

                read_clover(&L, cl_ldl_gpu, ix, DOWN, c_idx(i, i));
                B[i][n] /= creal(L);

                for (int k = i + 1; k < 2 * NF; k++) {
                    read_clover(&U, cl_ldl_gpu, ix, UP, c_idx(k, i));
                    U = conj(U);
                    A[i][n] -= U * A[k][n];

                    read_clover(&L, cl_ldl_gpu, ix, DOWN, c_idx(k, i));
                    L = conj(L);
                    B[i][n] -= L * B[k][n];
                }
            }
        }

        // Construct force matrices
        for (int i = 0; i < NF; i++) {
            for (int j = 0; j < NF; j++) {
                a21 = A[i + NF][j];
                a12 = conj(A[j + NF][i]);
                a43 = B[i + NF][j];
                a34 = conj(B[j + NF][i]);

                if (i < j) {
                    a11 = conj(A[j][i]);
                    a22 = conj(A[j + NF][i + NF]);
                    a33 = conj(B[j][i]);
                    a44 = conj(B[j + NF][i + NF]);
                } else {
                    a11 = A[i][j];
                    a22 = A[i + NF][j + NF];
                    a33 = B[i][j];
                    a44 = B[i + NF][j + NF];
                }

#ifdef REPR_IS_REAL
                read_force(&c, cl_force_gpu, ix, 0, i, j);
                c += cimag(a12) + cimag(a21) - cimag(a34) - cimag(a43); // X_01
                write_force(&c, cl_force_gpu, ix, 0, i, j);

                read_force(&c, cl_force_gpu, ix, 1, i, j);
                c += creal(a12) - creal(a21) + creal(a43) - creal(a34); // X_02
                write_force(&c, cl_force_gpu, ix, 1, i, j);

                read_force(&c, cl_force_gpu, ix, 2, i, j);
                c += cimag(a22) - cimag(a11) + cimag(a44) - cimag(a33); // X_12
                write_force(&c, cl_force_gpu, ix, 2, i, j);

                read_force(&c, cl_force_gpu, ix, 3, i, j);
                c += cimag(a11) - cimag(a22) + cimag(a44) - cimag(a33); // X_03
                write_force(&c, cl_force_gpu, ix, 3, i, j);

                read_force(&c, cl_force_gpu, ix, 4, i, j);
                c += creal(a12) - creal(a21) + creal(a34) - creal(a43); // X_13
                write_force(&c, cl_force_gpu, ix, 4, i, j);

                read_force(&c, cl_force_gpu, ix, 5, i, j);
                c -= cimag(a12) + cimag(a21) + cimag(a34) + cimag(a43); // X_23*/
                write_force(&c, cl_force_gpu, ix, 5, i, j);
#else
                read_force(&c, cl_force_gpu, ix, 0, i, j);
                c -= I * (a12 + a21 - a34 - a43); // X_01
                write_force(&c, cl_force_gpu, ix, 0, i, j);

                read_force(&c, cl_force_gpu, ix, 1, i, j);
                c += a12 - a21 + a43 - a34; // X_02
                write_force(&c, cl_force_gpu, ix, 1, i, j);

                read_force(&c, cl_force_gpu, ix, 2, i, j);
                c -= I * (a22 - a11 + a44 - a33); // X_12
                write_force(&c, cl_force_gpu, ix, 2, i, j);

                read_force(&c, cl_force_gpu, ix, 3, i, j);
                c -= I * (a11 - a22 + a44 - a33); // X_03
                write_force(&c, cl_force_gpu, ix, 3, i, j);

                read_force(&c, cl_force_gpu, ix, 4, i, j);
                c += a12 - a21 + a34 - a43; // X_13
                write_force(&c, cl_force_gpu, ix, 4, i, j);

                read_force(&c, cl_force_gpu, ix, 5, i, j);
                c += I * (a12 + a21 + a34 + a43); // X_23*/
                write_force(&c, cl_force_gpu, ix, 5, i, j);
#endif
            }
        }
    }
}

__global__ static void _clover_la_logdet(double nf, double *la_gpu, ldl_t *cl_ldl_gpu, int N, int block_start) {
    for (int id = blockDim.x * blockIdx.x + threadIdx.x; id < N; id += gridDim.x * blockDim.x) {
        const int ix = id + block_start;
        double prod = 1;
        hr_complex up, dn;

        for (int n = 0; n < 2 * NF; n++) {
            read_clover(&up, cl_ldl_gpu, ix, UP, c_idx(n, n));
            read_clover(&dn, cl_ldl_gpu, ix, DOWN, c_idx(n, n));
            prod *= creal(up) * creal(dn);
        }

        double val;
        read_gpu<double>(0, &val, la_gpu, ix, 0, 1);
        val -= nf * log(prod);
        write_gpu<double>(0, &val, la_gpu, ix, 0, 1);
    }
}

__global__ static void _compute_clover_term(suNfc *cl_term_gpu, double csw_value, suNf *gauge, int *iup_gpu, int *idn_gpu,
                                            int N, int block_start) {
    for (int id = blockDim.x * blockIdx.x + threadIdx.x; id < N; id += gridDim.x * blockDim.x) {
        const int ix = id + block_start;
        suNf tmp[6];
        double csw;

        hr_complex atmp;
        hr_complex btmp;
        hr_complex ctmp;
        hr_complex dtmp;
        hr_complex c;

        csw = csw_value;
        csw = -csw / 16.0;

        clover_loop(ix, 0, 1, &tmp[0], gauge, iup_gpu, idn_gpu);
        clover_loop(ix, 0, 2, &tmp[1], gauge, iup_gpu, idn_gpu);
        clover_loop(ix, 0, 3, &tmp[2], gauge, iup_gpu, idn_gpu);
        clover_loop(ix, 1, 2, &tmp[3], gauge, iup_gpu, idn_gpu);
        clover_loop(ix, 1, 3, &tmp[4], gauge, iup_gpu, idn_gpu);
        clover_loop(ix, 2, 3, &tmp[5], gauge, iup_gpu, idn_gpu);

        for (int i = 0; i < NF; i++) {
            for (int j = 0; j < NF; j++) {
                int ij = i * NF + j;
                int ji = j * NF + i;

#ifdef REPR_IS_REAL
                atmp = I * (tmp[2].c[ji] - tmp[2].c[ij]);
                btmp = I * (tmp[3].c[ij] - tmp[3].c[ji]);
                ctmp = tmp[1].c[ji] - tmp[1].c[ij] + I * (tmp[0].c[ji] - tmp[0].c[ij]);
                dtmp = tmp[4].c[ij] - tmp[4].c[ji] + I * (tmp[5].c[ji] - tmp[5].c[ij]);

#else
                atmp = I * (conj(tmp[2].c[ji]) - tmp[2].c[ij]);
                btmp = I * (conj(tmp[3].c[ij]) - tmp[3].c[ji]);
                ctmp = I * (conj(tmp[0].c[ji]) - tmp[0].c[ij]) - tmp[1].c[ij] + conj(tmp[1].c[ji]);
                dtmp = tmp[4].c[ij] - conj(tmp[4].c[ji]) + I * (conj(tmp[5].c[ji]) - tmp[5].c[ij]);
#endif
                c = csw * (atmp - conj(btmp));
                write_clover_term_comp(&c, cl_term_gpu, ix, 0, i, j);
                c = csw * (ctmp - dtmp);
                write_clover_term_comp(&c, cl_term_gpu, ix, 1, i, j);
                c = -csw * (atmp + conj(btmp));
                write_clover_term_comp(&c, cl_term_gpu, ix, 2, i, j);
                c = -csw * (ctmp + dtmp);
                write_clover_term_comp(&c, cl_term_gpu, ix, 3, i, j);
            }
        }
    }
}

double get_csw_gpu() {
    return csw_value;
}

void compute_ldl_decomp_gpu(double sigma0) {
    if (sigma == sigma0) {
        return;
    } else {
        sigma = sigma0;
    }

    _CUDA_CALL((&glattice), grid, (N * NF * NF), block_start, ixp,
               (_compute_ldl_decomp<<<grid, BLOCK_SIZE, 0, 0>>>(cl_term->gpu_ptr, cl_ldl->gpu_ptr, sigma, N, block_start)));

    _CUDA_CALL((&glattice), grid, N, block_start, ixp, (_ldl<<<grid, BLOCK_SIZE, 0, 0>>>(cl_ldl->gpu_ptr, 0, N, block_start)));

    _CUDA_CALL((&glattice), grid, N, block_start, ixp, (_ldl<<<grid, BLOCK_SIZE, 0, 0>>>(cl_ldl->gpu_ptr, 1, N, block_start)));
}

void compute_clover_term_gpu() {
    sigma = 0xF00F;
    start_sendrecv_suNf_field(u_gauge_f);
    complete_sendrecv_suNf_field(u_gauge_f);
    _CUDA_CALL((&glattice), grid, N, block_start, ixp,
               (_compute_clover_term<<<grid, BLOCK_SIZE, 0, 0>>>(cl_term->gpu_ptr, csw_value, u_gauge_f->gpu_ptr, iup_gpu,
                                                                 idn_gpu, N, block_start)););

    apply_BCs_on_clover_term(cl_term);
}

void clover_la_logdet_gpu(double nf, double mass, scalar_field *la) {
    compute_ldl_decomp_gpu(4.0 + mass);

    _CUDA_CALL((&glat_odd), grid, N, block_start, ixp,
               (_clover_la_logdet<<<grid, BLOCK_SIZE, 0, 0>>>(nf, la->gpu_ptr, cl_ldl->gpu_ptr, N, block_start)));
}

void compute_force_logdet_gpu(double mass, double coeff) {
    compute_ldl_decomp_gpu(4.0 + mass);

    _CUDA_CALL((&glat_odd), grid, N, block_start, ixp,
               (_compute_clover_force<<<grid, BLOCK_SIZE, 0, 0>>>(cl_ldl->gpu_ptr, cl_force->gpu_ptr, coeff, N, block_start)));
}

void clover_init_gpu(double csw) {
    cl_term = alloc_clover_term(&glattice);
    cl_ldl = alloc_ldl_field(&glattice);
    cl_force = alloc_clover_force(&glattice);

    cudaMemset(cl_term->gpu_ptr, 0, 4 * sizeof(suNfc) * glattice.gsize_gauge);
    cudaMemset(cl_ldl->gpu_ptr, 0, sizeof(ldl_t) * glattice.gsize_gauge);
    cudaMemset(cl_force->gpu_ptr, 0, 6 * sizeof(suNf) * glattice.gsize_gauge);

    sigma = 0xF00F;
    csw_value = csw;
    lprintf("CLOVER", 10, "Initial Coefficient: csw = %1.6f\n", csw_value);
}

void set_csw_gpu(double *csw) {
    csw_value = *csw;
    lprintf("CLOVER", 10, "Coefficient: reset to csw = %1.6f\n", csw_value);
}

double (*get_csw)(void) = get_csw_gpu;
void (*compute_ldl_decomp)(double) = compute_ldl_decomp_gpu;
void (*compute_clover_term)(void) = compute_clover_term_gpu;
void (*clover_la_logdet)(double, double, scalar_field *) = clover_la_logdet_gpu;
void (*compute_force_logdet)(double, double) = compute_force_logdet_gpu;
void (*clover_init)(double) = clover_init_gpu;
void (*set_csw)(double *) = set_csw_gpu;

#endif
#endif
