/***************************************************************************\
* Copyright (c) 2023, Sofie Martins                                         *   
* All rights reserved.                                                      * 
\***************************************************************************/

#ifndef DPHI_GPU_KERNELS_HPP
#define DPHI_GPU_KERNELS_HPP

#include "./Dphi_gpu_twisted_bc.h"
#include "geometry.h"
#include "libhr_core.h"
#include "utils.h"

__device__ __constant__ char UP_MASK = T_UP_MASK | X_UP_MASK | Y_UP_MASK | Z_UP_MASK;
__device__ __constant__ char DN_MASK = T_DN_MASK | X_DN_MASK | Y_DN_MASK | Z_DN_MASK;
__device__ __constant__ char T_MASK = T_UP_MASK | T_DN_MASK;
__device__ __constant__ char X_MASK = X_UP_MASK | X_DN_MASK;
__device__ __constant__ char Y_MASK = Y_UP_MASK | Y_DN_MASK;
__device__ __constant__ char Z_MASK = Z_UP_MASK | Z_DN_MASK;

#define find_neighbor(input, _ix, _dir, _mu) ((_dir == UP) ? input->iup_gpu[4 * (_ix) + _mu] : input->idn_gpu[4 * (_ix) + _mu])

#define _FIND_BUFFER_DIRECTION(_ix, _iy, _mu, _dir, _piece, _input)            \
    _iy = blockIdx.x * BLOCK_SIZE + threadIdx.x + _input->base_in[_piece - 1]; \
    const char DIR_MASK = _input->imask_gpu[iy];                               \
    _mu = _MU(DIR_MASK);                                                       \
    const int dir_inverted = _DIR(DIR_MASK);                                   \
    _ix = find_neighbor(_input, _iy, dir_inverted, _mu);                       \
    _dir = !dir_inverted;

#define _LOOP_DIRECTIONS(_ix, _iy, _mu, _dir, _input, body)  \
    for (_mu = 0; _mu < 4; ++_mu) {                          \
        for (_dir = UP; _dir <= DOWN; ++_dir) {              \
            const char DIR_MASK = MASK(_mu, _dir);           \
            if (_input->imask_gpu[_ix] & DIR_MASK) {         \
                _iy = find_neighbor(_input, _ix, _dir, _mu); \
                body;                                        \
            }                                                \
        }                                                    \
    }

#define iup_on_gpu(_dir) int __idx_in_global = iup_d[4 * (__idx_out_global) + _dir]
#define idn_on_gpu(_dir) int __idx_in_global = idn_d[4 * (__idx_out_global) + _dir]
#define MASK(_mu, _dir) (1u << (2 * _mu + _dir));

#define _DIR(MASK) ((MASK & UP_MASK) ? UP : DOWN)
#define _MU(MASK) ((MASK & T_MASK) ? 0 : (MASK & X_MASK) ? 1 : (MASK & Y_MASK) ? 2 : 3)

template <typename HSPINOR_TYPE, class REAL, typename GAUGE_TYPE, typename SITE_TYPE>
__device__ void evaluate_direction(SITE_TYPE *r, SITE_TYPE *in, GAUGE_TYPE *gauge, int ix, int iy, int mu, int dir,
                                   int master_shift) {
    HSPINOR_TYPE sn;
    GAUGE_TYPE u;
    if (mu == 0) {
        if (dir == UP) {
            in_spinor_field<REAL>(&(sn.c[0]), in, iy - master_shift, 0);
            in_spinor_field<REAL>(&(sn.c[1]), in, iy - master_shift, 2);
            in_gauge_field<REAL>(&u, gauge, ix, iy, 0, dir);

            _vector_add_assign_f(sn.c[0], sn.c[1]);

            _suNf_theta_T_multiply(sn.c[1], u, sn.c[0]);

            _vector_add_assign_f((*r).c[0], sn.c[1]);
            _vector_add_assign_f((*r).c[2], sn.c[1]);

            in_spinor_field<REAL>(&(sn.c[0]), in, iy - master_shift, 1);
            in_spinor_field<REAL>(&(sn.c[1]), in, iy - master_shift, 3);
            _vector_add_assign_f(sn.c[0], sn.c[1]);

            _suNf_theta_T_multiply(sn.c[1], u, sn.c[0]);

            _vector_add_assign_f((*r).c[1], sn.c[1]);
            _vector_add_assign_f((*r).c[3], sn.c[1]);
        } else {
            in_spinor_field<REAL>(&(sn.c[0]), in, iy - master_shift, 0);
            in_spinor_field<REAL>(&(sn.c[1]), in, iy - master_shift, 2);
            in_gauge_field<REAL>(&u, gauge, ix, iy, 0, dir);

            _vector_sub_assign_f(sn.c[0], sn.c[1]);

            _suNf_theta_T_inverse_multiply(sn.c[1], u, sn.c[0]);

            _vector_add_assign_f((*r).c[0], sn.c[1]);
            _vector_sub_assign_f((*r).c[2], sn.c[1]);

            in_spinor_field<REAL>(&(sn.c[0]), in, iy - master_shift, 1);
            in_spinor_field<REAL>(&(sn.c[1]), in, iy - master_shift, 3);
            _vector_sub_assign_f(sn.c[0], sn.c[1]);

            _suNf_theta_T_inverse_multiply(sn.c[1], u, sn.c[0]);

            _vector_add_assign_f((*r).c[1], sn.c[1]);
            _vector_sub_assign_f((*r).c[3], sn.c[1]);
        }
    } else if (mu == 1) {
        if (dir == UP) {
            in_spinor_field<REAL>(&(sn.c[0]), in, iy - master_shift, 0);
            in_spinor_field<REAL>(&(sn.c[1]), in, iy - master_shift, 3);
            in_gauge_field<REAL>(&u, gauge, ix, iy, 1, dir);

            _vector_i_add_assign_f(sn.c[0], sn.c[1]);
            _suNf_theta_X_multiply(sn.c[1], u, sn.c[0]);

            _vector_add_assign_f((*r).c[0], sn.c[1]);
            _vector_i_sub_assign_f((*r).c[3], sn.c[1]);

            in_spinor_field<REAL>(&(sn.c[0]), in, iy - master_shift, 1);
            in_spinor_field<REAL>(&(sn.c[1]), in, iy - master_shift, 2);
            _vector_i_add_assign_f(sn.c[0], sn.c[1]);

            _suNf_theta_X_multiply(sn.c[1], u, sn.c[0]);

            _vector_add_assign_f((*r).c[1], sn.c[1]);
            _vector_i_sub_assign_f((*r).c[2], sn.c[1]);
        } else {
            in_spinor_field<REAL>(&(sn.c[0]), in, iy - master_shift, 0);
            in_spinor_field<REAL>(&(sn.c[1]), in, iy - master_shift, 3);
            in_gauge_field<REAL>(&u, gauge, ix, iy, 1, dir);

            _vector_i_sub_assign_f(sn.c[0], sn.c[1]);
            _suNf_theta_X_inverse_multiply(sn.c[1], u, sn.c[0]);

            _vector_add_assign_f((*r).c[0], sn.c[1]);
            _vector_i_add_assign_f((*r).c[3], sn.c[1]);

            in_spinor_field<REAL>(&(sn.c[0]), in, iy - master_shift, 1);
            in_spinor_field<REAL>(&(sn.c[1]), in, iy - master_shift, 2);
            _vector_i_sub_assign_f(sn.c[0], sn.c[1]);

            _suNf_theta_X_inverse_multiply(sn.c[1], u, sn.c[0]);

            _vector_add_assign_f((*r).c[1], sn.c[1]);
            _vector_i_add_assign_f((*r).c[2], sn.c[1]);
        }
    } else if (mu == 2) {
        if (dir == UP) {
            in_spinor_field<REAL>(&(sn.c[0]), in, iy - master_shift, 0);
            in_spinor_field<REAL>(&(sn.c[1]), in, iy - master_shift, 3);
            _vector_add_assign_f(sn.c[0], sn.c[1]);

            in_gauge_field<REAL>(&u, gauge, ix, iy, 2, dir);
            _suNf_theta_Y_multiply(sn.c[1], u, sn.c[0]);

            _vector_add_assign_f((*r).c[0], sn.c[1]);
            _vector_add_assign_f((*r).c[3], sn.c[1]);

            in_spinor_field<REAL>(&(sn.c[0]), in, iy - master_shift, 1);
            in_spinor_field<REAL>(&(sn.c[1]), in, iy - master_shift, 2);
            _vector_sub_assign_f(sn.c[0], sn.c[1]);

            _suNf_theta_Y_multiply(sn.c[1], u, sn.c[0]);

            _vector_add_assign_f((*r).c[1], sn.c[1]);
            _vector_sub_assign_f((*r).c[2], sn.c[1]);
        } else {
            in_spinor_field<REAL>(&(sn.c[0]), in, iy - master_shift, 0);
            in_spinor_field<REAL>(&(sn.c[1]), in, iy - master_shift, 3);
            _vector_sub_assign_f(sn.c[0], sn.c[1]);

            in_gauge_field<REAL>(&u, gauge, ix, iy, 2, dir);
            _suNf_theta_Y_inverse_multiply(sn.c[1], u, sn.c[0]);

            _vector_add_assign_f((*r).c[0], sn.c[1]);
            _vector_sub_assign_f((*r).c[3], sn.c[1]);

            in_spinor_field<REAL>(&(sn.c[0]), in, iy - master_shift, 1);
            in_spinor_field<REAL>(&(sn.c[1]), in, iy - master_shift, 2);
            _vector_add_assign_f(sn.c[0], sn.c[1]);

            _suNf_theta_Y_inverse_multiply(sn.c[1], u, sn.c[0]);

            _vector_add_assign_f((*r).c[1], sn.c[1]);
            _vector_add_assign_f((*r).c[2], sn.c[1]);
        }
    } else if (mu == 3) {
        if (dir == UP) {
            in_spinor_field<REAL>(&(sn.c[0]), in, iy - master_shift, 0);
            in_spinor_field<REAL>(&(sn.c[1]), in, iy - master_shift, 2);
            _vector_i_add_assign_f(sn.c[0], sn.c[1]);

            in_gauge_field<REAL>(&u, gauge, ix, iy, 3, dir);
            _suNf_theta_Z_multiply(sn.c[1], u, sn.c[0]);

            _vector_add_assign_f((*r).c[0], sn.c[1]);
            _vector_i_sub_assign_f((*r).c[2], sn.c[1]);

            in_spinor_field<REAL>(&(sn.c[0]), in, iy - master_shift, 1);
            in_spinor_field<REAL>(&(sn.c[1]), in, iy - master_shift, 3);
            _vector_i_sub_assign_f(sn.c[0], sn.c[1]);

            _suNf_theta_Z_multiply(sn.c[1], u, sn.c[0]);

            _vector_add_assign_f((*r).c[1], sn.c[1]);
            _vector_i_add_assign_f((*r).c[3], sn.c[1]);
        } else {
            in_spinor_field<REAL>(&(sn.c[0]), in, iy - master_shift, 0);
            in_spinor_field<REAL>(&(sn.c[1]), in, iy - master_shift, 2);
            _vector_i_sub_assign_f(sn.c[0], sn.c[1]);

            in_gauge_field<REAL>(&u, gauge, ix, iy, 3, dir);
            _suNf_theta_Z_inverse_multiply(sn.c[1], u, sn.c[0]);

            _vector_add_assign_f((*r).c[0], sn.c[1]);
            _vector_i_add_assign_f((*r).c[2], sn.c[1]);

            in_spinor_field<REAL>(&(sn.c[0]), in, iy - master_shift, 1);
            in_spinor_field<REAL>(&(sn.c[1]), in, iy - master_shift, 3);
            _vector_i_add_assign_f(sn.c[0], sn.c[1]);

            _suNf_theta_Z_inverse_multiply(sn.c[1], u, sn.c[0]);

            _vector_add_assign_f((*r).c[1], sn.c[1]);
            _vector_i_sub_assign_f((*r).c[3], sn.c[1]);
        }
    }
}

template <typename HSPINOR_TYPE, class REAL, typename SITE_TYPE, typename GAUGE_TYPE>
__global__ void Dphi_gpu_inner_kernel(kernel_field_input *input) {
    SITE_TYPE r;
    SITE_TYPE *field_out = (SITE_TYPE *)input->field_out;
    SITE_TYPE *field_in = (SITE_TYPE *)input->field_in;
    GAUGE_TYPE *gauge = (GAUGE_TYPE *)input->gauge;

    int iy, mu, dir;
    _KERNEL_PIECE_FOR(piece) {
        _IF_IN_BOX_OUT(input, piece) {
            _spinor_zero_f(r);
            const int ix = blockIdx.x * BLOCK_SIZE + threadIdx.x + input->base_out[piece - 1];

            _LOOP_DIRECTIONS(ix, iy, mu, dir, input,
                             (evaluate_direction<HSPINOR_TYPE, REAL, GAUGE_TYPE, SITE_TYPE>(&r, field_in, gauge, ix, iy, mu,
                                                                                            dir, input->master_shift_in));)

            _spinor_mul_f(r, -0.5, r);
            int ix_spinor = ix - input->master_shift_out;
            write_out_spinor_field<REAL>(&r, field_out, ix_spinor);
        }
    }
}

// Cannot run two boundary kernels at the same time -> race condition
template <typename HSPINOR_TYPE, class REAL, typename SITE_TYPE, typename GAUGE_TYPE>
__global__ void Dphi_gpu_boundary_kernel(kernel_field_input *input) {
    SITE_TYPE r;
    SITE_TYPE res;
    SITE_TYPE *field_out = (SITE_TYPE *)input->field_out;
    SITE_TYPE *field_in = (SITE_TYPE *)input->field_in;
    GAUGE_TYPE *gauge = (GAUGE_TYPE *)input->gauge;

    int iy, ix, mu, dir;
    _KERNEL_PIECE_FOR(piece) {
        _IF_IN_BOX_IN(input, piece) {
            _FIND_BUFFER_DIRECTION(ix, iy, mu, dir, piece, input);
            _spinor_zero_f(r);
            evaluate_direction<HSPINOR_TYPE, REAL, GAUGE_TYPE, SITE_TYPE>(&r, field_in, gauge, ix, iy, mu, dir,
                                                                          input->master_shift_in);

            const int ix_spinor = ix - input->master_shift_out;
            in_spinor_field<REAL>(&res, field_out, ix_spinor, 0);
            _spinor_mul_add_assign_f(res, -0.5, r);
            write_out_spinor_field<REAL>(&res, field_out, ix_spinor);
        }
    }
}

#ifdef WITH_CLOVER
template <typename VECTOR_TYPE, class REAL, typename SITE_TYPE>
__global__ void Cphi_gpu_kernel_(SITE_TYPE *dptr, SITE_TYPE *sptr, suNfc *cl_term, double mass, int assign, int N,
                                 int block_start) {
    for (int ix = blockIdx.x * blockDim.x + threadIdx.x; ix < N; ix += gridDim.x * blockDim.x) {
        VECTOR_TYPE v1, v2;
        VECTOR_TYPE s0, s1, out;
        suNfc cl_0, cl_1;

        read_gpu<REAL>(0, &s0, sptr, ix, 0, 1);
        read_gpu<REAL>(0, &s1, sptr, ix, 1, 1);
        read_gpu<double>(0, &cl_0, cl_term, ix, 0, 4);
        read_gpu<double>(0, &cl_1, cl_term, ix, 1, 4);
        _suNfc_multiply(v1, cl_0, s0);
        _suNfc_multiply(v2, cl_1, s1);
        _vector_add_f(v1, v1, v2);
        _vector_mul_add_assign_f(v1, mass, s0);
        if (assign) {
            read_gpu<REAL>(0, &out, dptr, ix, 0, 1);
            _vector_add_assign_f(v1, out);
        }
        write_gpu<REAL>(0, &v1, dptr, ix, 0, 1);

        _suNfc_inverse_multiply(v1, cl_1, s0);
        _suNfc_multiply(v2, cl_0, s1);
        _vector_sub_f(v1, v1, v2);
        _vector_mul_add_assign_f(v1, mass, s1);
        if (assign) {
            read_gpu<REAL>(0, &out, dptr, ix, 1, 1);
            _vector_add_assign_f(v1, out);
        }
        write_gpu<REAL>(0, &v1, dptr, ix, 1, 1);

        read_gpu<REAL>(0, &s0, sptr, ix, 2, 1);
        read_gpu<REAL>(0, &s1, sptr, ix, 3, 1);
        read_gpu<double>(0, &cl_0, cl_term, ix, 2, 4);
        read_gpu<double>(0, &cl_1, cl_term, ix, 3, 4);
        _suNfc_multiply(v1, cl_0, s0);
        _suNfc_multiply(v2, cl_1, s1);
        _vector_add_f(v1, v1, v2);
        _vector_mul_add_assign_f(v1, mass, s0);
        if (assign) {
            read_gpu<REAL>(0, &out, dptr, ix, 2, 1);
            _vector_add_assign_f(v1, out);
        }
        write_gpu<REAL>(0, &v1, dptr, ix, 2, 1);

        _suNfc_inverse_multiply(v1, cl_1, s0);
        _suNfc_multiply(v2, cl_0, s1);
        _vector_sub_f(v1, v1, v2);
        _vector_mul_add_assign_f(v1, mass, s1);
        if (assign) {
            read_gpu<REAL>(0, &out, dptr, ix, 3, 1);
            _vector_add_assign_f(v1, out);
        }
        write_gpu<REAL>(0, &v1, dptr, ix, 3, 1);
    }
}

template <typename VECTOR_TYPE, typename COMPLEX, class REAL, typename SITE_TYPE>
__global__ void Cphi_inv_kernel_(SITE_TYPE *dptr, SITE_TYPE *sptr, ldl_t *ldl_gpu, int mass, int assign, int N,
                                 int block_start) {
    for (int ix = blockIdx.x * blockDim.x + threadIdx.x; ix < N; ix += gridDim.x * blockDim.x) {
        ldl_t site_ldl;
        hr_complex *up, *dn, c;
        SITE_TYPE out, in;
        int n;

        read_gpu<double>(0, &site_ldl, ldl_gpu, ix, 0, 1);
        read_gpu<REAL>(0, &in, sptr, ix, 0, 1);

        up = site_ldl.up;
        dn = site_ldl.dn;

        COMPLEX *x;
        x = (COMPLEX *)&in;
        for (int i = 0; i < 2 * NF; i++) {
            for (int k = 0; k < i; k++) {
                n = i * (i + 1) / 2 + k;
                _complex_mul_sub_assign(x[i], (COMPLEX)up[n], x[k]);
                _complex_mul_sub_assign(x[i + 2 * NF], (COMPLEX)dn[n], x[k + 2 * NF]);
            }
        }

        for (int i = 2 * NF - 1; i >= 0; i--) {
            n = i * (i + 1) / 2 + i;
            _complex_mulr(x[i], 1. / ((REAL)creal(up[n])), x[i]);
            _complex_mulr(x[i + 2 * NF], 1. / ((REAL)creal(dn[n])), x[i + 2 * NF]);
            for (int k = i + 1; k < 2 * NF; k++) {
                n = k * (k + 1) / 2 + i;

                c = (COMPLEX)conj(up[n]);
                _complex_mul_sub_assign(x[i], c, x[k]);
                c = (COMPLEX)conj(dn[n]);
                _complex_mul_sub_assign(x[i + 2 * NF], c, x[k + 2 * NF]);
            }
        }

        if (assign) {
            read_gpu<REAL>(0, &out, dptr, ix, 0, 1);
            _spinor_add_assign_f(in, out);
        }

        write_gpu<REAL>(0, &in, dptr, ix, 0, 1);
    }
}

#endif

#ifdef WITH_EXPCLOVER

template <typename VECTOR_TYPE, class REAL, typename SITE_TYPE>
__global__ void Cphi_gpu_kernel_(SITE_TYPE *dptr, SITE_TYPE *sptr, suNfc *cl_term, double mass, double invexpmass, int assign,
                                 int N, int block_start, int NN_loc, int NNexp_loc) {
    for (int ix = blockIdx.x * blockDim.x + threadIdx.x; ix < N; ix += gridDim.x * blockDim.x) {
        suNfc Aplus[4];
        suNfc Aminus[4];

        suNfc expAplus[4];
        suNfc expAminus[4];
        VECTOR_TYPE v1, v2;
        SITE_TYPE out, in, tmp;
        suNfc s0, s1, s2, s3;

        read_gpu<REAL>(0, &in, sptr, ix, 0, 1);
        read_gpu<double>(0, &s0, cl_term, ix, 0, 4);
        read_gpu<double>(0, &s1, cl_term, ix, 1, 4);
        read_gpu<double>(0, &s2, cl_term, ix, 2, 4);
        read_gpu<double>(0, &s3, cl_term, ix, 3, 4);

        _suNfc_mul(Aplus[0], invexpmass, s0);
        _suNfc_mul(Aplus[1], invexpmass, s1);
        _suNfc_dagger(Aplus[2], Aplus[1]);
        _suNfc_mul(Aplus[3], -invexpmass, s0);

        _suNfc_mul(Aminus[0], invexpmass, s2);
        _suNfc_mul(Aminus[1], invexpmass, s3);
        _suNfc_dagger(Aminus[2], Aminus[1]);
        _suNfc_mul(Aminus[3], -invexpmass, s2);

        clover_exp(Aplus, expAplus, NN_loc, NNexp_loc);
        clover_exp(Aminus, expAminus, NN_loc, NNexp_loc);

        _suNfc_mul_assign(expAplus[0], mass);
        _suNfc_mul_assign(expAplus[1], mass);
        _suNfc_mul_assign(expAplus[2], mass);
        _suNfc_mul_assign(expAplus[3], mass);
        _suNfc_mul_assign(expAminus[0], mass);
        _suNfc_mul_assign(expAminus[1], mass);
        _suNfc_mul_assign(expAminus[2], mass);
        _suNfc_mul_assign(expAminus[3], mass);

        // Comp 0
        _suNfc_multiply(v1, expAplus[0], in.c[0]);
        _suNfc_multiply(v2, expAplus[1], in.c[1]);
        _vector_add_f(tmp.c[0], v1, v2);

        // Comp 1
        _suNfc_multiply(v1, expAplus[2], in.c[0]);
        _suNfc_multiply(v2, expAplus[3], in.c[1]);
        _vector_add_f(tmp.c[1], v1, v2);

        // Comp 2
        _suNfc_multiply(v1, expAminus[0], in.c[2]);
        _suNfc_multiply(v2, expAminus[1], in.c[3]);
        _vector_add_f(tmp.c[2], v1, v2);

        // Comp 3
        _suNfc_multiply(v1, expAminus[2], in.c[2]);
        _suNfc_multiply(v2, expAminus[3], in.c[3]);
        _vector_add_f(tmp.c[3], v1, v2);

        if (assign) {
            read_gpu<REAL>(0, &out, dptr, ix, 0, 1);
            _spinor_add_assign_f(tmp, out);
        }

        write_gpu<REAL>(0, &tmp, dptr, ix, 0, 1);
    }
}

#endif
#endif