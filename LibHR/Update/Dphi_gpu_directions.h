/***************************************************************************\
* Copyright (c) 2008, 2022, Claudio Pica, Sofie Martins                     *   
* All rights reserved.                                                      * 
\***************************************************************************/

#ifndef DIRECTIONS_H
#define DIRECTIONS_H

#define inner_direction(__macro, __mask, __ix_to_iy) \
      if (imask_gpu[__idx_out_global]&__mask) { \
            _find_index(__ix_to_iy); \
            __macro; \
      }

#define boundary_calculation(__macro, __mask, __ix_to_iy) \
      /*Don't invert, use buffer indices */ \
      /* Or: */ \
      if (!(imask_gpu[__idx_out_global]&__mask)) { \
        _find_index(__ix_to_iy); \
        if (__idx_in_local < __stride_in && __idx_in_local >= 0) { \
            __macro; \
        } \
      }

#define _T_plus \
      _IN_SPINOR_FIELD(sn.c[0], 0); \
      _IN_SPINOR_FIELD(sn.c[1], 2); \
      _OUT_GAUGE_FIELD(u, 0); \
      \
      _vector_add_assign_f(sn.c[0], sn.c[1]); \
      _suNf_theta_T_multiply(sn.c[1], u, sn.c[0]); \
      _vector_add_assign_f(r.c[0], sn.c[1]); \
      _vector_add_assign_f(r.c[2], sn.c[1]); \
      \
      _IN_SPINOR_FIELD(sn.c[0], 1); \
      _IN_SPINOR_FIELD(sn.c[1], 3); \
      \
      _vector_add_assign_f(sn.c[0], sn.c[1]); \
      _suNf_theta_T_multiply(sn.c[1], u, sn.c[0]); \
      _vector_add_assign_f(r.c[1], sn.c[1]); \
      _vector_add_assign_f(r.c[3], sn.c[1]); \

#define _T_minus \
      _IN_SPINOR_FIELD(sn.c[0], 0); \
      _IN_SPINOR_FIELD(sn.c[1], 2); \
      _IN_GAUGE_FIELD(u, 0); \
      \
      _vector_sub_assign_f(sn.c[0], sn.c[1]); \
      _suNf_theta_T_inverse_multiply(sn.c[1], u, sn.c[0]); \
      _vector_add_assign_f(r.c[0], sn.c[1]); \
      _vector_sub_assign_f(r.c[2], sn.c[1]); \
      \
      _IN_SPINOR_FIELD(sn.c[0], 1); \
      _IN_SPINOR_FIELD(sn.c[1], 3); \
      \
      _vector_sub_assign_f(sn.c[0], sn.c[1]); \
      _suNf_theta_T_inverse_multiply(sn.c[1], u, sn.c[0]); \
      _vector_add_assign_f(r.c[1], sn.c[1]); \
      _vector_sub_assign_f(r.c[3], sn.c[1]); 

#define _X_plus \
      _IN_SPINOR_FIELD(sn.c[0], 0); \
      _IN_SPINOR_FIELD(sn.c[1], 3); \
      _OUT_GAUGE_FIELD(u, 1); \
      \
      _vector_i_add_assign_f(sn.c[0], sn.c[1]); \
      _suNf_theta_X_multiply(sn.c[1], u, sn.c[0]); \
      _vector_add_assign_f(r.c[0], sn.c[1]); \
      _vector_i_sub_assign_f(r.c[3], sn.c[1]); \
      \
      _IN_SPINOR_FIELD(sn.c[0], 1); \
      _IN_SPINOR_FIELD(sn.c[1], 2); \
      \
      _vector_i_add_assign_f(sn.c[0], sn.c[1]); \
      _suNf_theta_X_multiply(sn.c[1], u, sn.c[0]); \
      _vector_add_assign_f(r.c[1], sn.c[1]); \
      _vector_i_sub_assign_f(r.c[2], sn.c[1]); 

#define _X_minus \
      _IN_SPINOR_FIELD(sn.c[0], 0); \
      _IN_SPINOR_FIELD(sn.c[1], 3); \
      _IN_GAUGE_FIELD(u, 1); \
      \
      _vector_i_sub_assign_f(sn.c[0], sn.c[1]);\
      _suNf_theta_X_inverse_multiply(sn.c[1], u, sn.c[0]);\
      _vector_add_assign_f(r.c[0], sn.c[1]); \
      _vector_i_add_assign_f(r.c[3], sn.c[1]); \
      \
      _IN_SPINOR_FIELD(sn.c[0], 1); \
      _IN_SPINOR_FIELD(sn.c[1], 2); \
      \
      _vector_i_sub_assign_f(sn.c[0], sn.c[1]); \
      _suNf_theta_X_inverse_multiply(sn.c[1], u, sn.c[0]); \
      _vector_add_assign_f(r.c[1], sn.c[1]); \
      _vector_i_add_assign_f(r.c[2], sn.c[1]);

#define _Y_plus \
      _IN_SPINOR_FIELD(sn.c[0], 0); \
      _IN_SPINOR_FIELD(sn.c[1], 3); \
      _OUT_GAUGE_FIELD(u, 2); \
      if (__stride_in < 80) printf("Spinor comp: %0.2e + i%0.2e\n", creal(u.c[3]), cimag(u.c[3])); \
      \
      _vector_add_assign_f(sn.c[0], sn.c[1]); \
      _suNf_theta_Y_multiply(sn.c[1], u, sn.c[0]); \
      _vector_add_assign_f(r.c[0], sn.c[1]); \
      _vector_add_assign_f(r.c[3], sn.c[1]); \
      \
      _IN_SPINOR_FIELD(sn.c[0], 1); \
      _IN_SPINOR_FIELD(sn.c[1], 2); \
      \
      _vector_sub_assign_f(sn.c[0], sn.c[1]); \
      _suNf_theta_Y_multiply(sn.c[1], u, sn.c[0]); \
      _vector_add_assign_f(r.c[1], sn.c[1]); \
      _vector_sub_assign_f(r.c[2], sn.c[1]);  \

#define _Y_minus \
      _IN_SPINOR_FIELD(sn.c[0], 0); \
      _IN_SPINOR_FIELD(sn.c[1], 3); \
      _IN_GAUGE_FIELD(u, 2); \
      \
      _vector_sub_assign_f(sn.c[0], sn.c[1]); \
      _suNf_theta_Y_inverse_multiply(sn.c[1], u, sn.c[0]); \
      _vector_add_assign_f(r.c[0], sn.c[1]); \
      _vector_sub_assign_f(r.c[3], sn.c[1]); \
      \
      _IN_SPINOR_FIELD(sn.c[0], 1); \
      _IN_SPINOR_FIELD(sn.c[1], 2); \
      \
      _vector_add_assign_f(sn.c[0], sn.c[1]); \
      _suNf_theta_Y_inverse_multiply(sn.c[1], u, sn.c[0]); \
      _vector_add_assign_f(r.c[1], sn.c[1]); \
      _vector_add_assign_f(r.c[2], sn.c[1]);

#define _Z_plus \
      _IN_SPINOR_FIELD(sn.c[0], 0); \
      _IN_SPINOR_FIELD(sn.c[1], 2); \
      _OUT_GAUGE_FIELD(u, 3); \
      \
      _vector_i_add_assign_f(sn.c[0], sn.c[1]); \
      _suNf_theta_Z_multiply(sn.c[1], u, sn.c[0]); \
      _vector_add_assign_f(r.c[0], sn.c[1]); \
      _vector_i_sub_assign_f(r.c[2], sn.c[1]); \
      \
      _IN_SPINOR_FIELD(sn.c[0], 1); \
      _IN_SPINOR_FIELD(sn.c[1], 3); \
      \
      _vector_i_sub_assign_f(sn.c[0], sn.c[1]); \
      _suNf_theta_Z_multiply(sn.c[1], u, sn.c[0]); \
      _vector_add_assign_f(r.c[1], sn.c[1]); \
      _vector_i_add_assign_f(r.c[3], sn.c[1]);

#define _Z_minus \
      _IN_SPINOR_FIELD(sn.c[0], 0); \
      _IN_SPINOR_FIELD(sn.c[1], 2); \
      _IN_GAUGE_FIELD(u, 3); \
      \
      _vector_i_sub_assign_f(sn.c[0], sn.c[1]); \
      _suNf_theta_Z_inverse_multiply(sn.c[1], u, sn.c[0]); \
      _vector_add_assign_f(r.c[0], sn.c[1]); \
      _vector_i_add_assign_f(r.c[2], sn.c[1]); \
      \
      _IN_SPINOR_FIELD(sn.c[0], 1); \
      _IN_SPINOR_FIELD(sn.c[1], 3); \
      \
      _vector_i_add_assign_f(sn.c[0], sn.c[1]); \
      _suNf_theta_Z_inverse_multiply(sn.c[1], u, sn.c[0]); \
      _vector_add_assign_f(r.c[1], sn.c[1]); \
      _vector_i_sub_assign_f(r.c[3], sn.c[1]);


#endif

