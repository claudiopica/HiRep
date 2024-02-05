/***************************************************************************\
* Copyright (c) 2023, Sofie Martins                                         *   
* All rights reserved.                                                      * 
\***************************************************************************/

#ifndef DPHI_GPU_KERNELS_HPP
#define DPHI_GPU_KERNELS_HPP

#include "geometry.h"
#include "libhr_core.h"
#include "update.h"
#include "utils.h"

#define DPHI_T_UP_GPU(ix, iy, in, gauge, r, sn, u)        \
    in_spinor_field<REAL>(&((sn).c[0]), (in), (iy), 0);   \
    in_spinor_field<REAL>(&((sn).c[1]), (in), (iy), 2);   \
    in_gauge_field<REAL>(&u, (gauge), (ix), (iy), 0, UP); \
                                                          \
    _vector_add_assign_f((sn).c[0], (sn).c[1]);           \
    _suNf_theta_T_multiply((sn).c[1], u, (sn).c[0]);      \
    _vector_mul_add_assign_f((r).c[0], -0.5, (sn).c[1]);  \
    _vector_mul_add_assign_f((r).c[2], -0.5, (sn).c[1]);  \
                                                          \
    in_spinor_field<REAL>(&((sn).c[0]), (in), (iy), 1);   \
    in_spinor_field<REAL>(&((sn).c[1]), (in), (iy), 3);   \
                                                          \
    _vector_add_assign_f((sn).c[0], (sn).c[1]);           \
    _suNf_theta_T_multiply((sn).c[1], u, (sn).c[0]);      \
    _vector_mul_add_assign_f((r).c[1], -0.5, (sn).c[1]);  \
    _vector_mul_add_assign_f((r).c[3], -0.5, (sn).c[1]);

#define DPHI_T_DN_GPU(ix, iy, in, gauge, r, sn, u)             \
    in_spinor_field<REAL>(&((sn).c[0]), (in), (iy), 0);        \
    in_spinor_field<REAL>(&((sn).c[1]), (in), (iy), 2);        \
    in_gauge_field<REAL>(&(u), (gauge), (ix), (iy), 0, DOWN);  \
                                                               \
    _vector_sub_assign_f((sn).c[0], (sn).c[1]);                \
    _suNf_theta_T_inverse_multiply((sn).c[1], (u), (sn).c[0]); \
    _vector_mul_add_assign_f((r).c[0], -0.5, (sn).c[1]);       \
    _vector_mul_sub_assign_f((r).c[2], -0.5, (sn).c[1]);       \
                                                               \
    in_spinor_field<REAL>(&((sn).c[0]), (in), (iy), 1);        \
    in_spinor_field<REAL>(&((sn).c[1]), (in), (iy), 3);        \
                                                               \
    _vector_sub_assign_f((sn).c[0], (sn).c[1]);                \
    _suNf_theta_T_inverse_multiply((sn).c[1], (u), (sn).c[0]); \
    _vector_mul_add_assign_f((r).c[1], -0.5, (sn).c[1]);       \
    _vector_mul_sub_assign_f((r).c[3], -0.5, (sn).c[1]);

#define DPHI_X_UP_GPU(ix, iy, in, gauge, r, sn, u)          \
    in_spinor_field<REAL>(&((sn).c[0]), (in), (iy), 0);     \
    in_spinor_field<REAL>(&((sn).c[1]), (in), (iy), 3);     \
    in_gauge_field<REAL>(&(u), (gauge), (ix), (iy), 1, UP); \
                                                            \
    _vector_i_add_assign_f((sn).c[0], (sn).c[1]);           \
    _suNf_theta_X_multiply((sn).c[1], (u), (sn).c[0]);      \
    _vector_mul_add_assign_f((r).c[0], -0.5, (sn).c[1]);    \
    _vector_i_mul_sub_assign_f((r).c[3], -0.5, (sn).c[1]);  \
                                                            \
    in_spinor_field<REAL>(&((sn).c[0]), (in), (iy), 1);     \
    in_spinor_field<REAL>(&((sn).c[1]), (in), (iy), 2);     \
                                                            \
    _vector_i_add_assign_f((sn).c[0], (sn).c[1]);           \
    _suNf_theta_X_multiply((sn).c[1], (u), (sn).c[0]);      \
    _vector_mul_add_assign_f((r).c[1], -0.5, (sn).c[1]);    \
    _vector_i_mul_sub_assign_f((r).c[2], -0.5, (sn).c[1]);

#define DPHI_X_DN_GPU(ix, iy, in, gauge, r, sn, u)             \
    in_spinor_field<REAL>(&((sn).c[0]), (in), (iy), 0);        \
    in_spinor_field<REAL>(&((sn).c[1]), (in), (iy), 3);        \
    in_gauge_field<REAL>(&(u), (gauge), (ix), (iy), 1, DOWN);  \
                                                               \
    _vector_i_sub_assign_f((sn).c[0], (sn).c[1]);              \
    _suNf_theta_X_inverse_multiply((sn).c[1], (u), (sn).c[0]); \
    _vector_mul_add_assign_f((r).c[0], -0.5, (sn).c[1]);       \
    _vector_i_mul_add_assign_f((r).c[3], -0.5, (sn).c[1]);     \
                                                               \
    in_spinor_field<REAL>(&((sn).c[0]), (in), (iy), 1);        \
    in_spinor_field<REAL>(&((sn).c[1]), (in), (iy), 2);        \
                                                               \
    _vector_i_sub_assign_f((sn).c[0], (sn).c[1]);              \
    _suNf_theta_X_inverse_multiply((sn).c[1], (u), (sn).c[0]); \
    _vector_mul_add_assign_f((r).c[1], -0.5, (sn).c[1]);       \
    _vector_i_mul_add_assign_f((r).c[2], -0.5, (sn).c[1]);

#define DPHI_Y_UP_GPU(ix, iy, in, gauge, r, sn, u)          \
    in_spinor_field<REAL>(&((sn).c[0]), (in), (iy), 0);     \
    in_spinor_field<REAL>(&((sn).c[1]), (in), (iy), 3);     \
    in_gauge_field<REAL>(&(u), (gauge), (ix), (iy), 2, UP); \
                                                            \
    _vector_add_assign_f((sn).c[0], (sn).c[1]);             \
    _suNf_theta_Y_multiply((sn).c[1], (u), (sn).c[0]);      \
    _vector_mul_add_assign_f((r).c[0], -0.5, (sn).c[1]);    \
    _vector_mul_add_assign_f((r).c[3], -0.5, (sn).c[1]);    \
                                                            \
    in_spinor_field<REAL>(&((sn).c[0]), (in), (iy), 1);     \
    in_spinor_field<REAL>(&((sn).c[1]), (in), (iy), 2);     \
                                                            \
    _vector_sub_assign_f((sn).c[0], (sn).c[1]);             \
    _suNf_theta_Y_multiply((sn).c[1], (u), (sn).c[0]);      \
    _vector_mul_add_assign_f((r).c[1], -0.5, (sn).c[1]);    \
    _vector_mul_sub_assign_f((r).c[2], -0.5, (sn).c[1]);

#define DPHI_Y_DN_GPU(ix, iy, in, gauge, r, sn, u)             \
    in_spinor_field<REAL>(&((sn).c[0]), (in), (iy), 0);        \
    in_spinor_field<REAL>(&((sn).c[1]), (in), (iy), 3);        \
    in_gauge_field<REAL>(&(u), (gauge), (ix), (iy), 2, DOWN);  \
                                                               \
    _vector_sub_assign_f((sn).c[0], (sn).c[1]);                \
    _suNf_theta_Y_inverse_multiply((sn).c[1], (u), (sn).c[0]); \
    _vector_mul_add_assign_f((r).c[0], -0.5, (sn).c[1]);       \
    _vector_mul_sub_assign_f((r).c[3], -0.5, (sn).c[1]);       \
                                                               \
    in_spinor_field<REAL>(&((sn).c[0]), (in), (iy), 1);        \
    in_spinor_field<REAL>(&((sn).c[1]), (in), (iy), 2);        \
                                                               \
    _vector_add_assign_f((sn).c[0], (sn).c[1]);                \
    _suNf_theta_Y_inverse_multiply((sn).c[1], (u), (sn).c[0]); \
    _vector_mul_add_assign_f((r).c[1], -0.5, (sn).c[1]);       \
    _vector_mul_add_assign_f((r).c[2], -0.5, (sn).c[1]);

#define DPHI_Z_UP_GPU(ix, iy, in, gauge, r, sn, u)          \
    in_spinor_field<REAL>(&((sn).c[0]), (in), (iy), 0);     \
    in_spinor_field<REAL>(&((sn).c[1]), (in), (iy), 2);     \
    in_gauge_field<REAL>(&(u), (gauge), (ix), (iy), 3, UP); \
                                                            \
    _vector_i_add_assign_f((sn).c[0], (sn).c[1]);           \
    _suNf_theta_Z_multiply((sn).c[1], (u), (sn).c[0]);      \
    _vector_mul_add_assign_f((r).c[0], -0.5, (sn).c[1]);    \
    _vector_i_mul_sub_assign_f((r).c[2], -0.5, (sn).c[1]);  \
                                                            \
    in_spinor_field<REAL>(&((sn).c[0]), (in), (iy), 1);     \
    in_spinor_field<REAL>(&((sn).c[1]), (in), (iy), 3);     \
                                                            \
    _vector_i_sub_assign_f((sn).c[0], (sn).c[1]);           \
    _suNf_theta_Z_multiply((sn).c[1], (u), (sn).c[0]);      \
    _vector_mul_add_assign_f((r).c[1], -0.5, (sn).c[1]);    \
    _vector_i_mul_add_assign_f((r).c[3], -0.5, (sn).c[1]);

#define DPHI_Z_DN_GPU(ix, iy, in, gauge, r, sn, u)             \
    in_spinor_field<REAL>(&((sn).c[0]), (in), (iy), 0);        \
    in_spinor_field<REAL>(&((sn).c[1]), (in), (iy), 2);        \
    in_gauge_field<REAL>(&(u), (gauge), (ix), (iy), 3, DOWN);  \
                                                               \
    _vector_i_sub_assign_f((sn).c[0], (sn).c[1]);              \
    _suNf_theta_Z_inverse_multiply((sn).c[1], (u), (sn).c[0]); \
    _vector_mul_add_assign_f((r).c[0], -0.5, (sn).c[1]);       \
    _vector_i_mul_add_assign_f((r).c[2], -0.5, (sn).c[1]);     \
                                                               \
    in_spinor_field<REAL>(&((sn).c[0]), (in), (iy), 1);        \
    in_spinor_field<REAL>(&((sn).c[1]), (in), (iy), 3);        \
                                                               \
    _vector_i_add_assign_f((sn).c[0], (sn).c[1]);              \
    _suNf_theta_Z_inverse_multiply((sn).c[1], (u), (sn).c[0]); \
    _vector_mul_add_assign_f((r).c[1], -0.5, (sn).c[1]);       \
    _vector_i_mul_sub_assign_f((r).c[3], -0.5, (sn).c[1]);

#define read_reduced(iy, in, sn, piece)                                        \
    do {                                                                       \
        const int block_offset = input->base_in[(piece)-1];                    \
        const HSPINOR_TYPE *in_offset = (HSPINOR_TYPE *)((in) + block_offset); \
        const int iy_loc = (iy)-block_offset;                                  \
        read_gpu<REAL>(0, &(sn), in_offset, iy_loc, 0, 1);                     \
    } while (0)

#define DPHI_RED_T_UP_GPU(r, sn)                     \
    _vector_mul_add_assign_f(r.c[0], -0.5, sn.c[0]); \
    _vector_mul_add_assign_f(r.c[2], -0.5, sn.c[0]); \
    _vector_mul_add_assign_f(r.c[1], -0.5, sn.c[1]); \
    _vector_mul_add_assign_f(r.c[3], -0.5, sn.c[1]);

#define DPHI_RED_T_DN_GPU(r, sn)                     \
    _vector_mul_add_assign_f(r.c[0], -0.5, sn.c[0]); \
    _vector_mul_sub_assign_f(r.c[2], -0.5, sn.c[0]); \
    _vector_mul_add_assign_f(r.c[1], -0.5, sn.c[1]); \
    _vector_mul_sub_assign_f(r.c[3], -0.5, sn.c[1]);

#define DPHI_RED_X_UP_GPU(r, sn)                       \
    _vector_mul_add_assign_f(r.c[0], -0.5, sn.c[0]);   \
    _vector_i_mul_sub_assign_f(r.c[3], -0.5, sn.c[0]); \
    _vector_mul_add_assign_f(r.c[1], -0.5, sn.c[1]);   \
    _vector_i_mul_sub_assign_f(r.c[2], -0.5, sn.c[1]);

#define DPHI_RED_X_DN_GPU(r, sn)                       \
    _vector_mul_add_assign_f(r.c[0], -0.5, sn.c[0]);   \
    _vector_i_mul_add_assign_f(r.c[3], -0.5, sn.c[0]); \
    _vector_mul_add_assign_f(r.c[1], -0.5, sn.c[1]);   \
    _vector_i_mul_add_assign_f(r.c[2], -0.5, sn.c[1]);

#define DPHI_RED_Y_UP_GPU(r, sn)                     \
    _vector_mul_add_assign_f(r.c[0], -0.5, sn.c[0]); \
    _vector_mul_add_assign_f(r.c[3], -0.5, sn.c[0]); \
    _vector_mul_add_assign_f(r.c[1], -0.5, sn.c[1]); \
    _vector_mul_sub_assign_f(r.c[2], -0.5, sn.c[1]);

#define DPHI_RED_Y_DN_GPU(r, sn)                     \
    _vector_mul_add_assign_f(r.c[0], -0.5, sn.c[0]); \
    _vector_mul_sub_assign_f(r.c[3], -0.5, sn.c[0]); \
    _vector_mul_add_assign_f(r.c[1], -0.5, sn.c[1]); \
    _vector_mul_add_assign_f(r.c[2], -0.5, sn.c[1]);

#define DPHI_RED_Z_UP_GPU(r, sn)                       \
    _vector_mul_add_assign_f(r.c[0], -0.5, sn.c[0]);   \
    _vector_i_mul_sub_assign_f(r.c[2], -0.5, sn.c[0]); \
    _vector_mul_add_assign_f(r.c[1], -0.5, sn.c[1]);   \
    _vector_i_mul_add_assign_f(r.c[3], -0.5, sn.c[1]);

#define DPHI_RED_Z_DN_GPU(r, sn)                       \
    _vector_mul_add_assign_f(r.c[0], -0.5, sn.c[0]);   \
    _vector_i_mul_add_assign_f(r.c[2], -0.5, sn.c[0]); \
    _vector_mul_add_assign_f(r.c[1], -0.5, sn.c[1]);   \
    _vector_i_mul_sub_assign_f(r.c[3], -0.5, sn.c[1]);

template <typename HSPINOR_TYPE, class REAL, typename COMPLEX, typename SITE_TYPE, typename GAUGE_TYPE>
__global__ void Dphi_gpu_inner_kernel(kernel_field_input *input) {
    _KERNEL_PIECE_FOR(piece) {
        if (input->gd_in & piece) {
            for (int id = blockIdx.x * blockDim.x + threadIdx.x; id < input->vol_out[piece - 1]; id += gridDim.x * blockDim.x) {
                int ix = id + input->base_out[piece - 1];
                SITE_TYPE *out = (SITE_TYPE *)input->field_out;
                SITE_TYPE *in = ((SITE_TYPE *)input->field_in);
                GAUGE_TYPE *gauge = (GAUGE_TYPE *)input->gauge;

                SITE_TYPE r;
                HSPINOR_TYPE sn;
                GAUGE_TYPE u;

                _spinor_zero_f(r);

                /******************************* direction +0 *********************************/
                if (input->imask_gpu[ix] & T_UP_MASK) {
                    const int iy = input->iup_gpu[4*ix];
                    DPHI_T_UP_GPU(ix, iy, in, gauge, r, sn, u);
                }

                /******************************* direction -0 *********************************/
                if (input->imask_gpu[ix] & T_DN_MASK) {
                    const int iy = input->idn_gpu[4*ix];
                    DPHI_T_DN_GPU(ix, iy, in, gauge, r, sn, u);
                }

                /******************************* direction +1 *********************************/
                if (input->imask_gpu[ix] & X_UP_MASK) {
                    const int iy = input->iup_gpu[4*ix+1];
                    DPHI_X_UP_GPU(ix, iy, in, gauge, r, sn, u);
                }

                /******************************* direction -1 *********************************/
                if (input->imask_gpu[ix] & X_DN_MASK) {
                    const int iy = input->idn_gpu[4*ix+1];
                    DPHI_X_DN_GPU(ix, iy, in, gauge, r, sn, u);
                }

                /******************************* direction +2 *********************************/
                if (input->imask_gpu[ix] & Y_UP_MASK) {
                    const int iy = input->iup_gpu[4*ix+2];
                    DPHI_Y_UP_GPU(ix, iy, in, gauge, r, sn, u);
                }

                /******************************* direction -2 *********************************/
                if (input->imask_gpu[ix] & Y_DN_MASK) {
                    const int iy = input->idn_gpu[4*ix+2];
                    DPHI_Y_DN_GPU(ix, iy, in, gauge, r, sn, u);
                }

                /******************************* direction +3 *********************************/
                if (input->imask_gpu[ix] & Z_UP_MASK) {
                    const int iy = input->iup_gpu[4*ix+3];
                    DPHI_Z_UP_GPU(ix, iy, in, gauge, r, sn, u);
                }

                /******************************* direction -3 *********************************/
                if (input->imask_gpu[ix] & Z_DN_MASK) {
                    const int iy = input->idn_gpu[4*ix+3];
                    DPHI_Z_DN_GPU(ix, iy, in, gauge, r, sn, u);
                }

                write_out_spinor_field<REAL>(&r, out, ix);
            }
        }
    }
}

// Cannot run two boundary kernels at the same time -> race condition
template <typename HSPINOR_TYPE, class REAL, typename SITE_TYPE, typename GAUGE_TYPE>
__global__ void Dphi_gpu_boundary_kernel(kernel_field_input *input) {
    _KERNEL_PIECE_FOR(piece) {
        if (input->gd_in & piece) {
            for (int id = blockIdx.x * blockDim.x + threadIdx.x; id < input->vol_in[piece - 1]; id += gridDim.x * blockDim.x) {
                int ix = 0;
                int iy = id + input->base_in[piece - 1];

                SITE_TYPE *out = (SITE_TYPE *)input->field_out;
                SITE_TYPE *in = (SITE_TYPE *)input->field_in;

                SITE_TYPE r;
                HSPINOR_TYPE sn;
                GAUGE_TYPE u;

                _spinor_zero_f(r);

                /******************************* direction +0 *********************************/
                if (input->mask & T_UP_MASK) {
                    ix = find_neighbor(input, iy, DOWN, 0);
                    read_reduced(iy, in, sn, piece);
                    DPHI_RED_T_UP_GPU(r, sn);
                }

                /******************************* direction -0 *********************************/
                if (input->mask & T_DN_MASK) {
                    ix = find_neighbor(input, iy, UP, 0);
                    read_reduced(iy, in, sn, piece);
                    DPHI_RED_T_DN_GPU(r, sn);
                }

                /******************************* direction +1 *********************************/
                if (input->mask & X_UP_MASK) {
                    ix = find_neighbor(input, iy, DOWN, 1);
                    read_reduced(iy, in, sn, piece);
                    DPHI_RED_X_UP_GPU(r, sn);
                }

                /******************************* direction -1 *********************************/
                if (input->mask & X_DN_MASK) {
                    ix = find_neighbor(input, iy, UP, 1);
                    read_reduced(iy, in, sn, piece);
                    DPHI_RED_X_DN_GPU(r, sn);
                }

                /******************************* direction +2 *********************************/
                if (input->mask & Y_UP_MASK) {
                    ix = find_neighbor(input, iy, DOWN, 2);
                    read_reduced(iy, in, sn, piece);
                    DPHI_RED_Y_UP_GPU(r, sn);
                }

                /******************************* direction -2 *********************************/
                if (input->mask & Y_DN_MASK) {
                    ix = find_neighbor(input, iy, UP, 2);
                    read_reduced(iy, in, sn, piece);
                    DPHI_RED_Y_DN_GPU(r, sn);
                }

                /******************************* direction +3 *********************************/
                if (input->mask & Z_UP_MASK) {
                    ix = find_neighbor(input, iy, DOWN, 3);
                    read_reduced(iy, in, sn, piece);
                    DPHI_RED_Z_UP_GPU(r, sn);
                }

                /******************************* direction -3 *********************************/
                if (input->mask & Z_DN_MASK) {
                    ix = find_neighbor(input, iy, UP, 3);
                    read_reduced(iy, in, sn, piece);
                    DPHI_RED_Z_DN_GPU(r, sn);
                }

                write_assign_gpu<REAL>(0, &r, out, ix, 0, 1);
            }
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

        clover_exp(Aplus, expAplus, NN_loc);
        clover_exp(Aminus, expAminus, NN_loc);

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