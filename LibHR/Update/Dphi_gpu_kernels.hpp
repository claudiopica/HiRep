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

#define SUBBLOCKS 32

#ifndef LARGE_N
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

#else

#define prefetch_gauge(ix, u_cpx, gauge, tgt_comp, vcomp, idn_gpu)                                      \
    read_gpu<REAL>(0, &(u_cpx[0]), gauge, ix, NF * tgt_comp + vcomp, 4);                                \
    read_gpu<REAL>(0, &(u_cpx[2]), gauge, ix, NF * tgt_comp + vcomp + NF * NF, 4);                      \
    read_gpu<REAL>(0, &(u_cpx[4]), gauge, ix, NF * tgt_comp + vcomp + 2 * NF * NF, 4);                  \
    read_gpu<REAL>(0, &(u_cpx[6]), gauge, ix, NF * tgt_comp + vcomp + 3 * NF * NF, 4);                  \
    read_gpu<REAL>(0, &(u_cpx[1]), gauge, idn_gpu[4 * ix], NF * vcomp + tgt_comp, 4);                   \
    read_gpu<REAL>(0, &(u_cpx[3]), gauge, idn_gpu[4 * ix + 1], NF * vcomp + tgt_comp + NF * NF, 4);     \
    read_gpu<REAL>(0, &(u_cpx[5]), gauge, idn_gpu[4 * ix + 2], NF * vcomp + tgt_comp + 2 * NF * NF, 4); \
    read_gpu<REAL>(0, &(u_cpx[7]), gauge, idn_gpu[4 * ix + 3], NF * vcomp + tgt_comp + 3 * NF * NF, 4);

#define write_out(r0, r1, r2, r3, out, ix, tgt_comp)        \
    _complex_mul(r0, -0.5, r0);                             \
    _complex_mul(r1, -0.5, r1);                             \
    _complex_mul(r2, -0.5, r2);                             \
    _complex_mul(r3, -0.5, r3);                             \
    write_gpu<REAL>(0, &r0, out, ix, tgt_comp, 1);          \
    write_gpu<REAL>(0, &r1, out, ix, tgt_comp + NF, 1);     \
    write_gpu<REAL>(0, &r2, out, ix, tgt_comp + 2 * NF, 1); \
    write_gpu<REAL>(0, &r3, out, ix, tgt_comp + 3 * NF, 1);

#define DPHI_T_UP_GPU(ix, iy, vcomp, in, gauge, r0, r1, r2, r3, sn_cpx0, u_cpx) \
    read_gpu<REAL>(0, &(sn_cpx0), in, iy, vcomp, 1);                            \
    read_assign_gpu<REAL>(0, &(sn_cpx0), in, iy, vcomp + 2 * NF, 1);            \
    _complex_mul_assign(r0, u_cpx[0], sn_cpx0);                                 \
    _complex_mul_assign(r2, u_cpx[0], sn_cpx0);                                 \
    read_gpu<REAL>(0, &(sn_cpx0), in, iy, vcomp + NF, 1);                       \
    read_assign_gpu<REAL>(0, &(sn_cpx0), in, iy, vcomp + 3 * NF, 1);            \
    _complex_mul_assign(r1, u_cpx[0], sn_cpx0);                                 \
    _complex_mul_assign(r3, u_cpx[0], sn_cpx0);

#ifdef REPR_IS_REAL
#define DPHI_T_DN_GPU(ix, iy, vcomp, in, gauge, r0, r1, r2, r3, sn_cpx0, u_cpx) \
    read_gpu<REAL>(0, &(sn_cpx0), in, iy, vcomp, 1);                            \
    read_sub_assign_gpu<REAL>(0, &(sn_cpx0), in, iy, vcomp + 2 * NF, 1);        \
    _complex_mul_assign(r0, (u_cpx[1]), sn_cpx0);                               \
    _complex_mul_sub_assign(r2, (u_cpx[1]), sn_cpx0);                           \
    read_gpu<REAL>(0, &(sn_cpx0), in, iy, vcomp + NF, 1);                       \
    read_sub_assign_gpu<REAL>(0, &(sn_cpx0), in, iy, vcomp + 3 * NF, 1);        \
    _complex_mul_assign(r1, (u_cpx[1]), sn_cpx0);                               \
    _complex_mul_sub_assign(r3, (u_cpx[1]), sn_cpx0);
#else
#define DPHI_T_DN_GPU(ix, iy, vcomp, in, gauge, r0, r1, r2, r3, sn_cpx0, u_cpx) \
    read_gpu<REAL>(0, &(sn_cpx0), in, iy, vcomp, 1);                            \
    read_sub_assign_gpu<REAL>(0, &(sn_cpx0), in, iy, vcomp + 2 * NF, 1);        \
    _complex_mul_assign(r0, conj(u_cpx[1]), sn_cpx0);                           \
    _complex_mul_sub_assign(r2, conj(u_cpx[1]), sn_cpx0);                       \
    read_gpu<REAL>(0, &(sn_cpx0), in, iy, vcomp + NF, 1);                       \
    read_sub_assign_gpu<REAL>(0, &(sn_cpx0), in, iy, vcomp + 3 * NF, 1);        \
    _complex_mul_assign(r1, conj(u_cpx[1]), sn_cpx0);                           \
    _complex_mul_sub_assign(r3, conj(u_cpx[1]), sn_cpx0);
#endif

#define DPHI_X_UP_GPU(ix, iy, vcomp, in, gauge, r0, r1, r2, r3, sn_cpx0, u_cpx) \
    read_gpu<REAL>(0, &(sn_cpx0), in, iy, vcomp + 3 * NF, 1);                   \
    _complex_i_plus(sn_cpx0, sn_cpx0);                                          \
    read_assign_gpu<REAL>(0, &(sn_cpx0), in, iy, vcomp, 1);                     \
    _complex_mul_assign(r0, u_cpx[2], sn_cpx0);                                 \
    _complex_i_minus(sn_cpx0, sn_cpx0);                                         \
    _complex_mul_assign(r3, u_cpx[2], sn_cpx0);                                 \
    read_gpu<REAL>(0, &(sn_cpx0), in, iy, vcomp + 2 * NF, 1);                   \
    _complex_i_plus(sn_cpx0, sn_cpx0);                                          \
    read_assign_gpu<REAL>(0, &(sn_cpx0), in, iy, vcomp + NF, 1);                \
    _complex_mul_assign(r1, u_cpx[2], sn_cpx0);                                 \
    _complex_i_minus(sn_cpx0, sn_cpx0);                                         \
    _complex_mul_assign(r2, u_cpx[2], sn_cpx0);

#ifdef REPR_IS_REAL
#define DPHI_X_DN_GPU(ix, iy, vcomp, in, gauge, r0, r1, r2, r3, sn_cpx0, u_cpx) \
    read_gpu<REAL>(0, &(sn_cpx0), in, iy, vcomp + 3 * NF, 1);                   \
    _complex_i_minus(sn_cpx0, sn_cpx0);                                         \
    read_assign_gpu<REAL>(0, &(sn_cpx0), in, iy, vcomp, 1);                     \
    _complex_mul_assign(r0, (u_cpx[3]), sn_cpx0);                               \
    _complex_i_plus(sn_cpx0, sn_cpx0);                                          \
    _complex_mul_assign(r3, (u_cpx[3]), sn_cpx0);                               \
    read_gpu<REAL>(0, &(sn_cpx0), in, iy, vcomp + 2 * NF, 1);                   \
    _complex_i_minus(sn_cpx0, sn_cpx0);                                         \
    read_assign_gpu<REAL>(0, &(sn_cpx0), in, iy, vcomp + NF, 1);                \
    _complex_mul_assign(r1, (u_cpx[3]), sn_cpx0);                               \
    _complex_i_plus(sn_cpx0, sn_cpx0);                                          \
    _complex_mul_assign(r2, (u_cpx[3]), sn_cpx0);

#else
#define DPHI_X_DN_GPU(ix, iy, vcomp, in, gauge, r0, r1, r2, r3, sn_cpx0, u_cpx) \
    read_gpu<REAL>(0, &(sn_cpx0), in, iy, vcomp + 3 * NF, 1);                   \
    _complex_i_minus(sn_cpx0, sn_cpx0);                                         \
    read_assign_gpu<REAL>(0, &(sn_cpx0), in, iy, vcomp, 1);                     \
    _complex_mul_assign(r0, conj(u_cpx[3]), sn_cpx0);                           \
    _complex_i_plus(sn_cpx0, sn_cpx0);                                          \
    _complex_mul_assign(r3, conj(u_cpx[3]), sn_cpx0);                           \
    read_gpu<REAL>(0, &(sn_cpx0), in, iy, vcomp + 2 * NF, 1);                   \
    _complex_i_minus(sn_cpx0, sn_cpx0);                                         \
    read_assign_gpu<REAL>(0, &(sn_cpx0), in, iy, vcomp + NF, 1);                \
    _complex_mul_assign(r1, conj(u_cpx[3]), sn_cpx0);                           \
    _complex_i_plus(sn_cpx0, sn_cpx0);                                          \
    _complex_mul_assign(r2, conj(u_cpx[3]), sn_cpx0);

#endif

#define DPHI_Y_UP_GPU(ix, iy, vcomp, in, gauge, r0, r1, r2, r3, sn_cpx0, u_cpx) \
    read_gpu<REAL>(0, &(sn_cpx0), in, iy, vcomp, 1);                            \
    read_assign_gpu<REAL>(0, &(sn_cpx0), in, iy, vcomp + 3 * NF, 1);            \
    _complex_mul_assign(r0, u_cpx[4], sn_cpx0);                                 \
    _complex_mul_assign(r3, u_cpx[4], sn_cpx0);                                 \
    read_gpu<REAL>(0, &(sn_cpx0), in, iy, vcomp + NF, 1);                       \
    read_sub_assign_gpu<REAL>(0, &(sn_cpx0), in, iy, vcomp + 2 * NF, 1);        \
    _complex_mul_assign(r1, u_cpx[4], sn_cpx0);                                 \
    _complex_mul_sub_assign(r2, u_cpx[4], sn_cpx0);

#ifdef REPR_IS_REAL
#define DPHI_Y_DN_GPU(ix, iy, vcomp, in, gauge, r0, r1, r2, r3, sn_cpx0, u_cpx) \
    read_gpu<REAL>(0, &(sn_cpx0), in, iy, vcomp, 1);                            \
    read_sub_assign_gpu<REAL>(0, &(sn_cpx0), in, iy, vcomp + 3 * NF, 1);        \
    _complex_mul_assign(r0, (u_cpx[5]), sn_cpx0);                               \
    _complex_mul_sub_assign(r3, (u_cpx[5]), sn_cpx0);                           \
    read_gpu<REAL>(0, &(sn_cpx0), in, iy, vcomp + NF, 1);                       \
    read_assign_gpu<REAL>(0, &(sn_cpx0), in, iy, vcomp + 2 * NF, 1);            \
    _complex_mul_assign(r1, (u_cpx[5]), sn_cpx0);                               \
    _complex_mul_assign(r2, (u_cpx[5]), sn_cpx0);

#else
#define DPHI_Y_DN_GPU(ix, iy, vcomp, in, gauge, r0, r1, r2, r3, sn_cpx0, u_cpx) \
    read_gpu<REAL>(0, &(sn_cpx0), in, iy, vcomp, 1);                            \
    read_sub_assign_gpu<REAL>(0, &(sn_cpx0), in, iy, vcomp + 3 * NF, 1);        \
    _complex_mul_assign(r0, conj(u_cpx[5]), sn_cpx0);                           \
    _complex_mul_sub_assign(r3, conj(u_cpx[5]), sn_cpx0);                       \
    read_gpu<REAL>(0, &(sn_cpx0), in, iy, vcomp + NF, 1);                       \
    read_assign_gpu<REAL>(0, &(sn_cpx0), in, iy, vcomp + 2 * NF, 1);            \
    _complex_mul_assign(r1, conj(u_cpx[5]), sn_cpx0);                           \
    _complex_mul_assign(r2, conj(u_cpx[5]), sn_cpx0);

#endif

#define DPHI_Z_UP_GPU(ix, iy, vcomp, in, gauge, r0, r1, r2, r3, sn_cpx0, u_cpx) \
    read_gpu<REAL>(0, &(sn_cpx0), in, iy, vcomp + 2 * NF, 1);                   \
    _complex_i_plus(sn_cpx0, sn_cpx0);                                          \
    read_assign_gpu<REAL>(0, &(sn_cpx0), in, iy, vcomp, 1);                     \
    _complex_mul_assign(r0, u_cpx[6], sn_cpx0);                                 \
    _complex_i_minus(sn_cpx0, sn_cpx0);                                         \
    _complex_mul_assign(r2, u_cpx[6], sn_cpx0);                                 \
    read_gpu<REAL>(0, &(sn_cpx0), in, iy, vcomp + 3 * NF, 1);                   \
    _complex_i_minus(sn_cpx0, sn_cpx0);                                         \
    read_assign_gpu<REAL>(0, &(sn_cpx0), in, iy, vcomp + NF, 1);                \
    _complex_mul_assign(r1, u_cpx[6], sn_cpx0);                                 \
    _complex_i_plus(sn_cpx0, sn_cpx0);                                          \
    _complex_mul_assign(r3, u_cpx[6], sn_cpx0);

#ifdef REPR_IS_REAL
#define DPHI_Z_DN_GPU(ix, iy, vcomp, in, gauge, r0, r1, r2, r3, sn_cpx0, u_cpx) \
    read_gpu<REAL>(0, &(sn_cpx0), in, iy, vcomp + 2 * NF, 1);                   \
    _complex_i_minus(sn_cpx0, sn_cpx0);                                         \
    read_assign_gpu<REAL>(0, &(sn_cpx0), in, iy, vcomp, 1);                     \
    _complex_mul_assign(r0, (u_cpx[7]), sn_cpx0);                               \
    _complex_i_plus(sn_cpx0, sn_cpx0);                                          \
    _complex_mul_assign(r2, (u_cpx[7]), sn_cpx0);                               \
    read_gpu<REAL>(0, &(sn_cpx0), in, iy, vcomp + 3 * NF, 1);                   \
    _complex_i_plus(sn_cpx0, sn_cpx0);                                          \
    read_assign_gpu<REAL>(0, &(sn_cpx0), in, iy, vcomp + NF, 1);                \
    _complex_mul_assign(r1, (u_cpx[7]), sn_cpx0);                               \
    _complex_i_minus(sn_cpx0, sn_cpx0);                                         \
    _complex_mul_assign(r3, (u_cpx[7]), sn_cpx0);

#else
#define DPHI_Z_DN_GPU(ix, iy, vcomp, in, gauge, r0, r1, r2, r3, sn_cpx0, u_cpx) \
    read_gpu<REAL>(0, &(sn_cpx0), in, iy, vcomp + 2 * NF, 1);                   \
    _complex_i_minus(sn_cpx0, sn_cpx0);                                         \
    read_assign_gpu<REAL>(0, &(sn_cpx0), in, iy, vcomp, 1);                     \
    _complex_mul_assign(r0, conj(u_cpx[7]), sn_cpx0);                           \
    _complex_i_plus(sn_cpx0, sn_cpx0);                                          \
    _complex_mul_assign(r2, conj(u_cpx[7]), sn_cpx0);                           \
    read_gpu<REAL>(0, &(sn_cpx0), in, iy, vcomp + 3 * NF, 1);                   \
    _complex_i_plus(sn_cpx0, sn_cpx0);                                          \
    read_assign_gpu<REAL>(0, &(sn_cpx0), in, iy, vcomp + NF, 1);                \
    _complex_mul_assign(r1, conj(u_cpx[7]), sn_cpx0);                           \
    _complex_i_minus(sn_cpx0, sn_cpx0);                                         \
    _complex_mul_assign(r3, conj(u_cpx[7]), sn_cpx0);

#endif
#endif

#define read_reduced(iy, in, sn, piece, base_in)                               \
    do {                                                                       \
        const int block_offset = base_in;                                      \
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

#ifndef LARGE_N

template <typename HSPINOR_TYPE, class REAL, typename COMPLEX, typename SITE_TYPE, typename GAUGE_TYPE>
__global__ void Dphi_gpu_inner_kernel(SITE_TYPE *in, SITE_TYPE *out, const GAUGE_TYPE *gauge, int vol_in, int vol_out,
                                      int base_in, int base_out, gd_type piece, char *imask_gpu, int *iup_gpu, int *idn_gpu) {
    for (int id = blockIdx.x * blockDim.x + threadIdx.x; id < vol_out; id += gridDim.x * blockDim.x) {
        int ix = id + base_out;

        SITE_TYPE r;
        HSPINOR_TYPE sn;
        GAUGE_TYPE u;

        _spinor_zero_f(r);

        /******************************* direction +0 *********************************/
        if (imask_gpu[ix] & T_UP_MASK) {
            const int iy = iup_gpu[4 * ix];
            DPHI_T_UP_GPU(ix, iy, in, gauge, r, sn, u);
        }

        /******************************* direction -0 *********************************/
        if (imask_gpu[ix] & T_DN_MASK) {
            const int iy = idn_gpu[4 * ix];
            DPHI_T_DN_GPU(ix, iy, in, gauge, r, sn, u);
        }

        /******************************* direction +1 *********************************/
        if (imask_gpu[ix] & X_UP_MASK) {
            const int iy = iup_gpu[4 * ix + 1];
            DPHI_X_UP_GPU(ix, iy, in, gauge, r, sn, u);
        }

        /******************************* direction -1 *********************************/
        if (imask_gpu[ix] & X_DN_MASK) {
            const int iy = idn_gpu[4 * ix + 1];
            DPHI_X_DN_GPU(ix, iy, in, gauge, r, sn, u);
        }

        /******************************* direction +2 *********************************/
        if (imask_gpu[ix] & Y_UP_MASK) {
            const int iy = iup_gpu[4 * ix + 2];
            DPHI_Y_UP_GPU(ix, iy, in, gauge, r, sn, u);
        }

        /******************************* direction -2 *********************************/
        if (imask_gpu[ix] & Y_DN_MASK) {
            const int iy = idn_gpu[4 * ix + 2];
            DPHI_Y_DN_GPU(ix, iy, in, gauge, r, sn, u);
        }

        /******************************* direction +3 *********************************/
        if (imask_gpu[ix] & Z_UP_MASK) {
            const int iy = iup_gpu[4 * ix + 3];
            DPHI_Z_UP_GPU(ix, iy, in, gauge, r, sn, u);
        }

        /******************************* direction -3 *********************************/
        if (imask_gpu[ix] & Z_DN_MASK) {
            const int iy = idn_gpu[4 * ix + 3];
            DPHI_Z_DN_GPU(ix, iy, in, gauge, r, sn, u);
        }

        write_out_spinor_field<REAL>(&r, out, ix);
    }
}

#else

template <typename HSPINOR_TYPE, class REAL, typename COMPLEX, typename SITE_TYPE, typename GAUGE_TYPE>
__global__ void Dphi_gpu_inner_kernel(SITE_TYPE *in, SITE_TYPE *out, const GAUGE_TYPE *gauge, int vol_in, int vol_out,
                                      int base_in, int base_out, gd_type piece, char *imask_gpu, int *iup_gpu, int *idn_gpu) {
    for (int id = blockIdx.x * blockDim.x + threadIdx.x; id < vol_out * NF; id += gridDim.x * blockDim.x) {
        // This might degrade performance for unusually sized lattices
        int divider = SUBBLOCKS;
        while (vol_out % divider != 0) {
            divider--;
        }

        const int ix = id % divider + (id / (divider * NF)) * divider + base_out;
        const int tgt_comp = (id / divider) % NF;

        COMPLEX sn_cpx0;
#ifdef REPR_IS_REAL
        REAL u_cpx[8];
#else
        COMPLEX u_cpx[8];
#endif
        COMPLEX r0, r1, r2, r3;
        _complex_0(r0);
        _complex_0(r1);
        _complex_0(r2);
        _complex_0(r3);

        for (int vcomp = 0; vcomp < NF; vcomp++) {
            prefetch_gauge(ix, u_cpx, gauge, tgt_comp, vcomp, idn_gpu);

            /******************************* direction +0 *********************************/
            if (imask_gpu[ix] & T_UP_MASK) {
                const int iy = iup_gpu[4 * ix];
                DPHI_T_UP_GPU(ix, iy, vcomp, in, gauge, r0, r1, r2, r3, sn_cpx0, u_cpx);
            }

            /******************************* direction -0 *********************************/
            if (imask_gpu[ix] & T_DN_MASK) {
                const int iy = idn_gpu[4 * ix];
                DPHI_T_DN_GPU(ix, iy, vcomp, in, gauge, r0, r1, r2, r3, sn_cpx0, u_cpx);
            }

            /******************************* direction +1 *********************************/
            if (imask_gpu[ix] & X_UP_MASK) {
                const int iy = iup_gpu[4 * ix + 1];
                DPHI_X_UP_GPU(ix, iy, vcomp, in, gauge, r0, r1, r2, r3, sn_cpx0, u_cpx);
            }

            /******************************* direction -1 *********************************/
            if (imask_gpu[ix] & X_DN_MASK) {
                const int iy = idn_gpu[4 * ix + 1];
                DPHI_X_DN_GPU(ix, iy, vcomp, in, gauge, r0, r1, r2, r3, sn_cpx0, u_cpx);
            }

            /******************************* direction +2 *********************************/
            if (imask_gpu[ix] & Y_UP_MASK) {
                const int iy = iup_gpu[4 * ix + 2];
                DPHI_Y_UP_GPU(ix, iy, vcomp, in, gauge, r0, r1, r2, r3, sn_cpx0, u_cpx);
            }

            /******************************* direction -2 *********************************/
            if (imask_gpu[ix] & Y_DN_MASK) {
                const int iy = idn_gpu[4 * ix + 2];
                DPHI_Y_DN_GPU(ix, iy, vcomp, in, gauge, r0, r1, r2, r3, sn_cpx0, u_cpx);
            }

            /******************************* direction +3 *********************************/
            if (imask_gpu[ix] & Z_UP_MASK) {
                const int iy = iup_gpu[4 * ix + 3];
                DPHI_Z_UP_GPU(ix, iy, vcomp, in, gauge, r0, r1, r2, r3, sn_cpx0, u_cpx);
            }

            /******************************* direction -3 *********************************/
            if (imask_gpu[ix] & Z_DN_MASK) {
                const int iy = idn_gpu[4 * ix + 3];
                DPHI_Z_DN_GPU(ix, iy, vcomp, in, gauge, r0, r1, r2, r3, sn_cpx0, u_cpx);
            }
        }

        write_out(r0, r1, r2, r3, out, ix, tgt_comp);
    }
}
#endif

// Cannot run two boundary kernels at the same time -> race condition
template <typename HSPINOR_TYPE, class REAL, typename SITE_TYPE, typename GAUGE_TYPE>
__global__ void Dphi_gpu_boundary_kernel(SITE_TYPE *in, SITE_TYPE *out, const GAUGE_TYPE *gauge, int vol_in, int vol_out,
                                         int base_in, int base_out, gd_type piece, char mask, int *iup_gpu, int *idn_gpu) {
    for (int id = blockIdx.x * blockDim.x + threadIdx.x; id < vol_in; id += gridDim.x * blockDim.x) {
        int ix = 0;
        int iy = id + base_in;

        SITE_TYPE r;
        HSPINOR_TYPE sn;
        GAUGE_TYPE u;

        _spinor_zero_f(r);

        /******************************* direction +0 *********************************/
        if (mask & T_UP_MASK) {
            ix = find_neighbor(iy, DOWN, 0, iup_gpu, idn_gpu);
            read_reduced(iy, in, sn, piece, base_in);
            DPHI_RED_T_UP_GPU(r, sn);
        }

        /******************************* direction -0 *********************************/
        if (mask & T_DN_MASK) {
            ix = find_neighbor(iy, UP, 0, iup_gpu, idn_gpu);
            read_reduced(iy, in, sn, piece, base_in);
            DPHI_RED_T_DN_GPU(r, sn);
        }

        /******************************* direction +1 *********************************/
        if (mask & X_UP_MASK) {
            ix = find_neighbor(iy, DOWN, 1, iup_gpu, idn_gpu);
            read_reduced(iy, in, sn, piece, base_in);
            DPHI_RED_X_UP_GPU(r, sn);
        }

        /******************************* direction -1 *********************************/
        if (mask & X_DN_MASK) {
            ix = find_neighbor(iy, UP, 1, iup_gpu, idn_gpu);
            read_reduced(iy, in, sn, piece, base_in);
            DPHI_RED_X_DN_GPU(r, sn);
        }

        /******************************* direction +2 *********************************/
        if (mask & Y_UP_MASK) {
            ix = find_neighbor(iy, DOWN, 2, iup_gpu, idn_gpu);
            read_reduced(iy, in, sn, piece, base_in);
            DPHI_RED_Y_UP_GPU(r, sn);
        }

        /******************************* direction -2 *********************************/
        if (mask & Y_DN_MASK) {
            ix = find_neighbor(iy, UP, 2, iup_gpu, idn_gpu);
            read_reduced(iy, in, sn, piece, base_in);
            DPHI_RED_Y_DN_GPU(r, sn);
        }

        /******************************* direction +3 *********************************/
        if (mask & Z_UP_MASK) {
            ix = find_neighbor(iy, DOWN, 3, iup_gpu, idn_gpu);
            read_reduced(iy, in, sn, piece, base_in);
            DPHI_RED_Z_UP_GPU(r, sn);
        }

        /******************************* direction -3 *********************************/
        if (mask & Z_DN_MASK) {
            ix = find_neighbor(iy, UP, 3, iup_gpu, idn_gpu);
            read_reduced(iy, in, sn, piece, base_in);
            DPHI_RED_Z_DN_GPU(r, sn);
        }

        write_assign_gpu<REAL>(0, &r, out, ix, 0, 1);
    }
}

#ifdef WITH_CLOVER
template <typename VECTOR_TYPE, class REAL, typename COMPLEX, typename SITE_TYPE>
__global__ void Cphi_gpu_kernel_(SITE_TYPE *dptr, SITE_TYPE *sptr, suNfc *cl_term, double mass, int assign, int N,
                                 int block_start) {
    for (int id = blockIdx.x * blockDim.x + threadIdx.x; id < N * NF * 4; id += gridDim.x * blockDim.x) {
        const int dir = id / (N * NF);
        const int id2 = id % (N * NF);
        int divider = SUBBLOCKS;
        while (N % divider != 0) {
            divider--;
        }
        const int tgt_comp = (id2 / divider) % NF;
        const int ix = id2 % divider + (id2 / (divider * NF)) * divider;

        COMPLEX v1, v2, v3;
        COMPLEX s0, s1, out;
        COMPLEX cl_0, cl_1;

        if (dir == 0) {
            v3 = 0.0;

            for (int vcomp = 0; vcomp < NF; vcomp++) {
                read_gpu<REAL>(0, &s0, sptr, ix, vcomp, 1);
                read_gpu<REAL>(0, &s1, sptr, ix, vcomp + NF, 1);
                read_gpu<double>(0, &cl_0, cl_term, ix, NF * tgt_comp + vcomp, 4);
                read_gpu<double>(0, &cl_1, cl_term, ix, NF * tgt_comp + vcomp + NF * NF, 4);
                v1 = cl_0 * s0;
                v2 = cl_1 * s1;
                v3 += v1 + v2;
                if (vcomp == tgt_comp) { v3 += mass * s0; }
            }

            if (assign) {
                write_assign_gpu<REAL>(0, &v3, dptr, ix, tgt_comp, 1);
            } else {
                write_gpu<REAL>(0, &v3, dptr, ix, tgt_comp, 1);
            }
        }

        if (dir == 1) {
            v3 = 0.0;

            for (int vcomp = 0; vcomp < NF; vcomp++) {
                // this reread might not be ideal
                read_gpu<REAL>(0, &s0, sptr, ix, vcomp, 1);
                read_gpu<REAL>(0, &s1, sptr, ix, vcomp + NF, 1);
                read_gpu<double>(0, &cl_0, cl_term, ix, NF * tgt_comp + vcomp, 4);
                read_gpu<double>(0, &cl_1, cl_term, ix, NF * vcomp + tgt_comp + NF * NF, 4);
                v1 = conj(cl_1) * s0;
                v2 = cl_0 * s1;
                v3 += (v1 - v2);
                if (vcomp == tgt_comp) { v3 += mass * s1; }
            }

            if (assign) {
                write_assign_gpu<REAL>(0, &v3, dptr, ix, tgt_comp + NF, 1);
            } else {
                write_gpu<REAL>(0, &v3, dptr, ix, tgt_comp + NF, 1);
            }
        }

        if (dir == 2) {
            v3 = 0.0;

            for (int vcomp = 0; vcomp < NF; vcomp++) {
                read_gpu<REAL>(0, &s0, sptr, ix, vcomp + 2 * NF, 1);
                read_gpu<REAL>(0, &s1, sptr, ix, vcomp + 3 * NF, 1);
                read_gpu<double>(0, &cl_0, cl_term, ix, NF * tgt_comp + vcomp + 2 * NF * NF, 4);
                read_gpu<double>(0, &cl_1, cl_term, ix, NF * tgt_comp + vcomp + 3 * NF * NF, 4);
                v1 = cl_0 * s0;
                v2 = cl_1 * s1;
                v3 += v1 + v2;
                if (vcomp == tgt_comp) { v3 += mass * s0; }
            }

            if (assign) {
                write_assign_gpu<REAL>(0, &v3, dptr, ix, tgt_comp + 2 * NF, 1);
            } else {
                write_gpu<REAL>(0, &v3, dptr, ix, tgt_comp + 2 * NF, 1);
            }
        }

        if (dir == 3) {
            v3 = 0.0;

            for (int vcomp = 0; vcomp < NF; vcomp++) {
                read_gpu<REAL>(0, &s0, sptr, ix, vcomp + 2 * NF, 1);
                read_gpu<REAL>(0, &s1, sptr, ix, vcomp + 3 * NF, 1);
                read_gpu<double>(0, &cl_0, cl_term, ix, NF * tgt_comp + vcomp + 2 * NF * NF, 4);
                read_gpu<double>(0, &cl_1, cl_term, ix, NF * vcomp + tgt_comp + 3 * NF * NF, 4);
                v1 = conj(cl_1) * s0;
                v2 = cl_0 * s1;
                v3 += (v1 - v2);
                if (vcomp == tgt_comp) { v3 += mass * s1; }
            }

            if (assign) {
                write_assign_gpu<REAL>(0, &v3, dptr, ix, tgt_comp + 3 * NF, 1);
            } else {
                write_gpu<REAL>(0, &v3, dptr, ix, tgt_comp + 3 * NF, 1);
            }
        }
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