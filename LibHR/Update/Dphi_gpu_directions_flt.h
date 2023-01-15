/***************************************************************************\
* Copyright (c) 2008, 2022, Claudio Pica, Sofie Martins                     *   
* All rights reserved.                                                      * 
\***************************************************************************/

#ifndef DIRECTIONS_FLT_H
#define DIRECTIONS_FLT_H

#define inner_direction_flt(__macro, __mask, __ix, __iy) \
      if (imask_gpu[ix]&__mask) { \
            int local_iy = __iy - block_start_iyp; \
            suNf_hspinor_flt sn;\
            suNf_flt u;\
            __macro(vol4h, vol4h); \
      }

#define boundary_calculation_flt(__macro, __mask, __ix, __iy) \
      /*Don't invert, use buffer indices */ \
      /* Or: */ \
      if (!(imask_gpu[ix]&__mask)) { \
        int local_iy = __iy - start; \
        if (local_iy < buf_stride && local_iy >= 0) { \
            suNf_hspinor_flt sn;\
            suNf_flt u;\
            __macro(buf_stride, vol4h); \
        } \
      }

#define _T_plus_flt(__read_stride, __write_stride) \
      read_gpu_suNf_vector_flt(__read_stride, sn.c[0], in, local_iy, 0); \
      read_gpu_suNf_vector_flt(__read_stride, sn.c[1], in, local_iy, 2); \
      read_gpu_suNf_flt(__write_stride, u, gauge_ixp, local_ix, 0); \
\
      _vector_add_assign_f(sn.c[0], sn.c[1]); \
      _suNf_theta_T_multiply(sn.c[1], u, sn.c[0]); \
\
      _vector_add_assign_f(r.c[0], sn.c[1]); \
      _vector_add_assign_f(r.c[2], sn.c[1]); \
\
      read_gpu_suNf_vector_flt(__read_stride, sn.c[0], in, local_iy, 1); \
      read_gpu_suNf_vector_flt(__read_stride, sn.c[1], in, local_iy, 3); \
      _vector_add_assign_f(sn.c[0], sn.c[1]); \
\
      _suNf_theta_T_multiply(sn.c[1], u, sn.c[0]); \
\
      _vector_add_assign_f(r.c[1], sn.c[1]); \
      _vector_add_assign_f(r.c[3], sn.c[1]); \


#define _T_minus_flt(__read_stride, __write_stride) \
    read_gpu_suNf_vector_flt(__read_stride, sn.c[0], in, local_iy, 0); \
    read_gpu_suNf_vector_flt(__read_stride, sn.c[1], in, local_iy, 2); \
    read_gpu_suNf_flt(__read_stride, u, gauge_iyp, local_iy, 0); \
\
    _vector_sub_assign_f(sn.c[0], sn.c[1]); \
    _suNf_theta_T_inverse_multiply(sn.c[1], u, sn.c[0]); \
\
    _vector_add_assign_f(r.c[0], sn.c[1]); \
    _vector_sub_assign_f(r.c[2], sn.c[1]); \
\
    read_gpu_suNf_vector_flt(__read_stride, sn.c[0], in, local_iy, 1); \
    read_gpu_suNf_vector_flt(__read_stride, sn.c[1], in, local_iy, 3); \
    _vector_sub_assign_f(sn.c[0], sn.c[1]); \
\
    _suNf_theta_T_inverse_multiply(sn.c[1], u, sn.c[0]); \
\
    _vector_add_assign_f(r.c[1], sn.c[1]); \
    _vector_sub_assign_f(r.c[3], sn.c[1]); 

#define _X_plus_flt(__read_stride, __write_stride) \
      read_gpu_suNf_vector_flt(__read_stride, sn.c[0], in, local_iy, 0); \
      read_gpu_suNf_vector_flt(__read_stride, sn.c[1], in, local_iy, 3); \
      read_gpu_suNf_flt(__write_stride, u, gauge_ixp, local_ix, 1); \
\
      _vector_i_add_assign_f(sn.c[0], sn.c[1]); \
      _suNf_theta_X_multiply(sn.c[1], u, sn.c[0]); \
\
      _vector_add_assign_f(r.c[0], sn.c[1]); \
      _vector_i_sub_assign_f(r.c[3], sn.c[1]); \
\
      read_gpu_suNf_vector_flt(__read_stride, sn.c[0], in, local_iy, 1); \
      read_gpu_suNf_vector_flt(__read_stride, sn.c[1], in, local_iy, 2); \
      _vector_i_add_assign_f(sn.c[0], sn.c[1]); \
\
      _suNf_theta_X_multiply(sn.c[1], u, sn.c[0]); \
\
      _vector_add_assign_f(r.c[1], sn.c[1]); \
      _vector_i_sub_assign_f(r.c[2], sn.c[1]); 

#define _X_minus_flt(__read_stride, __write_stride) \
      read_gpu_suNf_vector_flt(__read_stride, sn.c[0], in, local_iy, 0); \
      read_gpu_suNf_vector_flt(__read_stride, sn.c[1], in, local_iy, 3); \
      read_gpu_suNf_flt(__read_stride, u, gauge_iyp, local_iy, 1); \
\
      _vector_i_sub_assign_f(sn.c[0], sn.c[1]);\
      _suNf_theta_X_inverse_multiply(sn.c[1], u, sn.c[0]);\
\
      _vector_add_assign_f(r.c[0], sn.c[1]); \
      _vector_i_add_assign_f(r.c[3], sn.c[1]); \
\
      read_gpu_suNf_vector_flt(__read_stride, sn.c[0], in, local_iy, 1); \
      read_gpu_suNf_vector_flt(__read_stride, sn.c[1], in, local_iy, 2); \
      _vector_i_sub_assign_f(sn.c[0], sn.c[1]); \
\
      _suNf_theta_X_inverse_multiply(sn.c[1], u, sn.c[0]); \
\
      _vector_add_assign_f(r.c[1], sn.c[1]); \
      _vector_i_add_assign_f(r.c[2], sn.c[1]); \

#define _Y_plus_flt(__read_stride, __write_stride) \
      read_gpu_suNf_vector_flt(__read_stride, sn.c[0], in, local_iy, 0); \
      read_gpu_suNf_vector_flt(__read_stride, sn.c[1], in, local_iy, 3); \
      _vector_add_assign_f(sn.c[0], sn.c[1]); \
\
      read_gpu_suNf_flt(__write_stride, u, gauge_ixp, local_ix, 2); \
      _suNf_theta_Y_multiply(sn.c[1], u, sn.c[0]); \
\
      _vector_add_assign_f(r.c[0], sn.c[1]); \
      _vector_add_assign_f(r.c[3], sn.c[1]); \
\
      read_gpu_suNf_vector_flt(__read_stride, sn.c[0], in, local_iy, 1); \
      read_gpu_suNf_vector_flt(__read_stride, sn.c[1], in, local_iy, 2); \
      _vector_sub_assign_f(sn.c[0], sn.c[1]); \
\
      _suNf_theta_Y_multiply(sn.c[1], u, sn.c[0]); \
\
      _vector_add_assign_f(r.c[1], sn.c[1]); \
      _vector_sub_assign_f(r.c[2], sn.c[1]);

#define _Y_minus_flt(__read_stride, __write_stride) \
      read_gpu_suNf_vector_flt(__read_stride, sn.c[0], in, local_iy, 0); \
      read_gpu_suNf_vector_flt(__read_stride, sn.c[1], in, local_iy, 3); \
      _vector_sub_assign_f(sn.c[0], sn.c[1]); \
\
      read_gpu_suNf_flt(__read_stride, u, gauge_iyp, local_iy, 2); \
      _suNf_theta_Y_inverse_multiply(sn.c[1], u, sn.c[0]); \
\
      _vector_add_assign_f(r.c[0], sn.c[1]); \
      _vector_sub_assign_f(r.c[3], sn.c[1]); \
\
      read_gpu_suNf_vector_flt(__read_stride, sn.c[0], in, local_iy, 1); \
      read_gpu_suNf_vector_flt(__read_stride, sn.c[1], in, local_iy, 2); \
      _vector_add_assign_f(sn.c[0], sn.c[1]); \
\
      _suNf_theta_Y_inverse_multiply(sn.c[1], u, sn.c[0]); \
\
      _vector_add_assign_f(r.c[1], sn.c[1]); \
      _vector_add_assign_f(r.c[2], sn.c[1]); 

#define _Z_plus_flt(__read_stride, __write_stride) \
      read_gpu_suNf_vector_flt(__read_stride, sn.c[0], in, local_iy, 0); \
      read_gpu_suNf_vector_flt(__read_stride, sn.c[1], in, local_iy, 2); \
      _vector_i_add_assign_f(sn.c[0], sn.c[1]); \
\
      read_gpu_suNf_flt(__write_stride, u, gauge_ixp, local_ix, 3); \
      _suNf_theta_Z_multiply(sn.c[1], u, sn.c[0]); \
\
      _vector_add_assign_f(r.c[0], sn.c[1]); \
      _vector_i_sub_assign_f(r.c[2], sn.c[1]); \
\
      read_gpu_suNf_vector_flt(__read_stride, sn.c[0], in, local_iy, 1); \
      read_gpu_suNf_vector_flt(__read_stride, sn.c[1], in, local_iy, 3); \
      _vector_i_sub_assign_f(sn.c[0], sn.c[1]); \
\
      _suNf_theta_Z_multiply(sn.c[1], u, sn.c[0]); \
\
      _vector_add_assign_f(r.c[1], sn.c[1]); \
      _vector_i_add_assign_f(r.c[3], sn.c[1]); \

#define _Z_minus_flt(__read_stride, __write_stride) \
      read_gpu_suNf_vector_flt(__read_stride, sn.c[0], in, local_iy, 0); \
      read_gpu_suNf_vector_flt(__read_stride, sn.c[1], in, local_iy, 2); \
      _vector_i_sub_assign_f(sn.c[0], sn.c[1]); \
\
      read_gpu_suNf_flt(__read_stride, u, gauge_iyp, local_iy, 3); \
      _suNf_theta_Z_inverse_multiply(sn.c[1], u, sn.c[0]); \
\
      _vector_add_assign_f(r.c[0], sn.c[1]); \
      _vector_i_add_assign_f(r.c[2], sn.c[1]); \
\
      read_gpu_suNf_vector_flt(__read_stride, sn.c[0], in, local_iy, 1); \
      read_gpu_suNf_vector_flt(__read_stride, sn.c[1], in, local_iy, 3); \
      _vector_i_add_assign_f(sn.c[0], sn.c[1]); \
\
      _suNf_theta_Z_inverse_multiply(sn.c[1], u, sn.c[0]); \
\
      _vector_add_assign_f(r.c[1], sn.c[1]); \
      _vector_i_sub_assign_f(r.c[3], sn.c[1]);

#endif

