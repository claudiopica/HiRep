/***************************************************************************\
* Copyright (c) 2022, Sofie Martins                                         *   
* All rights reserved.                                                      * 
\***************************************************************************/

#ifndef DPHI_GPU_KERNELS_HPP
#define DPHI_GPU_KERNELS_HPP

#include "./Dphi_gpu_twisted_bc.h"
#include "geometry.h"
#include "libhr_core.h"

__device__ __constant__ char UP_MASK=T_UP_MASK | X_UP_MASK | Y_UP_MASK | Z_UP_MASK;
__device__ __constant__ char DN_MASK=T_DN_MASK | X_DN_MASK | Y_DN_MASK | Z_DN_MASK;
__device__ __constant__ char T_MASK=T_UP_MASK | T_DN_MASK;
__device__ __constant__ char X_MASK=X_UP_MASK | X_DN_MASK;
__device__ __constant__ char Y_MASK=Y_UP_MASK | Y_DN_MASK;
__device__ __constant__ char Z_MASK=Z_UP_MASK | Z_DN_MASK;


#define find_neighbor(input, _ix, _dir, _mu) ((_dir == UP) ? input->iup_gpu[4*(_ix) + _mu] : input->idn_gpu[4*(_ix) + _mu])

#define _FIND_BUFFER_DIRECTION(_ix, _iy, _mu, _dir, _piece, _input) \
    _iy = blockIdx.x * BLOCK_SIZE + threadIdx.x + _input->base_in[_piece-1]; \
    const char DIR_MASK = _input->imask_gpu[iy]; \
    _mu = _MU(DIR_MASK); \
    const int dir_inverted = _DIR(DIR_MASK); \
    _ix = find_neighbor(_input, _iy, dir_inverted, _mu); \
    _dir = !dir_inverted;

#define _LOOP_DIRECTIONS(_ix, _iy, _mu, _dir, _input, body) \
    for (_mu = 0; _mu < 4; ++_mu) { \
        for (_dir = UP; _dir <= DOWN; ++_dir) { \
            const char DIR_MASK = MASK(_mu, _dir); \
            if (_input->imask_gpu[_ix]&DIR_MASK) { \
                _iy = find_neighbor(_input, _ix, _dir, _mu); \
                body; \
            } \
        } \
    }

#define iup_on_gpu(_dir) int __idx_in_global = iup_d[4*(__idx_out_global) + _dir]
#define idn_on_gpu(_dir) int __idx_in_global = idn_d[4*(__idx_out_global) + _dir]
#define MASK(_mu, _dir) (1u << (2*_mu + _dir));

#define _DIR(MASK) ((MASK&UP_MASK) ? UP : DOWN)
#define _MU(MASK) ((MASK&T_MASK) ? 0 : (MASK&X_MASK) ? 1 : (MASK&Y_MASK) ? 2 : 3) 

template<typename HSPINOR_TYPE, class REAL, typename GAUGE_TYPE, typename SITE_TYPE>
__device__ void evaluate_direction(SITE_TYPE *r, SITE_TYPE *in, GAUGE_TYPE *gauge, int ix, int iy, int mu, int dir, int master_shift) {
    HSPINOR_TYPE sn;
    GAUGE_TYPE u;
    if (mu == 0) {
        if (dir == UP) {
            in_spinor_field<REAL>(&(sn.c[0]), in, iy-master_shift, 0);
            in_spinor_field<REAL>(&(sn.c[1]), in, iy-master_shift, 2);
            in_gauge_field<REAL>(&u, gauge, ix, iy, 0, dir);

            _vector_add_assign_f(sn.c[0], sn.c[1]);

            _suNf_theta_T_multiply(sn.c[1], u, sn.c[0]);

            _vector_add_assign_f((*r).c[0], sn.c[1]);
            _vector_add_assign_f((*r).c[2], sn.c[1]);

            in_spinor_field<REAL>(&(sn.c[0]), in, iy-master_shift, 1);
            in_spinor_field<REAL>(&(sn.c[1]), in, iy-master_shift, 3);
            _vector_add_assign_f(sn.c[0], sn.c[1]);

            _suNf_theta_T_multiply(sn.c[1], u, sn.c[0]);

            _vector_add_assign_f((*r).c[1], sn.c[1]);
            _vector_add_assign_f((*r).c[3], sn.c[1]); 
        } else {
            in_spinor_field<REAL>(&(sn.c[0]), in, iy-master_shift, 0);
            in_spinor_field<REAL>(&(sn.c[1]), in, iy-master_shift, 2);
            in_gauge_field<REAL>(&u, gauge, ix, iy, 0, dir);

            _vector_sub_assign_f(sn.c[0], sn.c[1]);  
            _suNf_theta_T_inverse_multiply(sn.c[1], u, sn.c[0]);  
        
            _vector_add_assign_f((*r).c[0], sn.c[1]);  
            _vector_sub_assign_f((*r).c[2], sn.c[1]);  
        
            in_spinor_field<REAL>(&(sn.c[0]), in, iy-master_shift, 1);  
            in_spinor_field<REAL>(&(sn.c[1]), in, iy-master_shift, 3);  
            _vector_sub_assign_f(sn.c[0], sn.c[1]);  
        
            _suNf_theta_T_inverse_multiply(sn.c[1], u, sn.c[0]);  
        
            _vector_add_assign_f((*r).c[1], sn.c[1]);  
            _vector_sub_assign_f((*r).c[3], sn.c[1]); 
        }
    } else if (mu == 1) {
        if (dir == UP) {
            in_spinor_field<REAL>(&(sn.c[0]), in, iy-master_shift, 0);  
            in_spinor_field<REAL>(&(sn.c[1]), in, iy-master_shift, 3);  
            in_gauge_field<REAL>(&u, gauge, ix, iy, 1, dir);  
        
            _vector_i_add_assign_f(sn.c[0], sn.c[1]);  
            _suNf_theta_X_multiply(sn.c[1], u, sn.c[0]);  
        
            _vector_add_assign_f((*r).c[0], sn.c[1]);  
            _vector_i_sub_assign_f((*r).c[3], sn.c[1]);  
        
            in_spinor_field<REAL>(&(sn.c[0]), in, iy-master_shift, 1);  
            in_spinor_field<REAL>(&(sn.c[1]), in, iy-master_shift, 2);  
            _vector_i_add_assign_f(sn.c[0], sn.c[1]);  
        
            _suNf_theta_X_multiply(sn.c[1], u, sn.c[0]);  
        
            _vector_add_assign_f((*r).c[1], sn.c[1]);  
            _vector_i_sub_assign_f((*r).c[2], sn.c[1]); 
        } else {
            in_spinor_field<REAL>(&(sn.c[0]), in, iy-master_shift, 0);  
            in_spinor_field<REAL>(&(sn.c[1]), in, iy-master_shift, 3);  
            in_gauge_field<REAL>(&u, gauge, ix, iy, 1, dir);  
        
            _vector_i_sub_assign_f(sn.c[0], sn.c[1]); 
            _suNf_theta_X_inverse_multiply(sn.c[1], u, sn.c[0]); 
        
            _vector_add_assign_f((*r).c[0], sn.c[1]);  
            _vector_i_add_assign_f((*r).c[3], sn.c[1]);  
        
            in_spinor_field<REAL>(&(sn.c[0]), in, iy-master_shift, 1);  
            in_spinor_field<REAL>(&(sn.c[1]), in, iy-master_shift, 2);  
            _vector_i_sub_assign_f(sn.c[0], sn.c[1]);  
        
            _suNf_theta_X_inverse_multiply(sn.c[1], u, sn.c[0]);  
        
            _vector_add_assign_f((*r).c[1], sn.c[1]);  
            _vector_i_add_assign_f((*r).c[2], sn.c[1]); 
        }
    } else if (mu == 2) {
        if (dir == UP) {
            in_spinor_field<REAL>(&(sn.c[0]), in, iy-master_shift, 0);  
            in_spinor_field<REAL>(&(sn.c[1]), in, iy-master_shift, 3);  
            _vector_add_assign_f(sn.c[0], sn.c[1]);  
        
            in_gauge_field<REAL>(&u, gauge, ix, iy, 2, dir);  
            _suNf_theta_Y_multiply(sn.c[1], u, sn.c[0]);  
        
            _vector_add_assign_f((*r).c[0], sn.c[1]);  
            _vector_add_assign_f((*r).c[3], sn.c[1]);  
        
            in_spinor_field<REAL>(&(sn.c[0]), in, iy-master_shift, 1);  
            in_spinor_field<REAL>(&(sn.c[1]), in, iy-master_shift, 2);  
            _vector_sub_assign_f(sn.c[0], sn.c[1]);  
        
            _suNf_theta_Y_multiply(sn.c[1], u, sn.c[0]);  
        
            _vector_add_assign_f((*r).c[1], sn.c[1]);  
            _vector_sub_assign_f((*r).c[2], sn.c[1]);
        } else {
            in_spinor_field<REAL>(&(sn.c[0]), in, iy-master_shift, 0);  
            in_spinor_field<REAL>(&(sn.c[1]), in, iy-master_shift, 3);  
            _vector_sub_assign_f(sn.c[0], sn.c[1]);  
        
            in_gauge_field<REAL>(&u, gauge, ix, iy, 2, dir);  
            _suNf_theta_Y_inverse_multiply(sn.c[1], u, sn.c[0]);  
        
            _vector_add_assign_f((*r).c[0], sn.c[1]);  
            _vector_sub_assign_f((*r).c[3], sn.c[1]);  
        
            in_spinor_field<REAL>(&(sn.c[0]), in, iy-master_shift, 1);  
            in_spinor_field<REAL>(&(sn.c[1]), in, iy-master_shift, 2);  
            _vector_add_assign_f(sn.c[0], sn.c[1]);  
        
            _suNf_theta_Y_inverse_multiply(sn.c[1], u, sn.c[0]);  
        
            _vector_add_assign_f((*r).c[1], sn.c[1]);  
            _vector_add_assign_f((*r).c[2], sn.c[1]);
        }
    } else if (mu == 3) {
        if (dir == UP) {
            in_spinor_field<REAL>(&(sn.c[0]), in, iy-master_shift, 0);  
            in_spinor_field<REAL>(&(sn.c[1]), in, iy-master_shift, 2);  
            _vector_i_add_assign_f(sn.c[0], sn.c[1]);  
        
            in_gauge_field<REAL>(&u, gauge, ix, iy, 3, dir);  
            _suNf_theta_Z_multiply(sn.c[1], u, sn.c[0]);  
        
            _vector_add_assign_f((*r).c[0], sn.c[1]);  
            _vector_i_sub_assign_f((*r).c[2], sn.c[1]);  
        
            in_spinor_field<REAL>(&(sn.c[0]), in, iy-master_shift, 1);  
            in_spinor_field<REAL>(&(sn.c[1]), in, iy-master_shift, 3);  
            _vector_i_sub_assign_f(sn.c[0], sn.c[1]);  
        
            _suNf_theta_Z_multiply(sn.c[1], u, sn.c[0]);  
        
            _vector_add_assign_f((*r).c[1], sn.c[1]);  
            _vector_i_add_assign_f((*r).c[3], sn.c[1]);  
        } else {
            in_spinor_field<REAL>(&(sn.c[0]), in, iy-master_shift, 0);  
            in_spinor_field<REAL>(&(sn.c[1]), in, iy-master_shift, 2);  
            _vector_i_sub_assign_f(sn.c[0], sn.c[1]);  
        
            in_gauge_field<REAL>(&u, gauge, ix, iy, 3, dir);  
            _suNf_theta_Z_inverse_multiply(sn.c[1], u, sn.c[0]);  
        
            _vector_add_assign_f((*r).c[0], sn.c[1]);  
            _vector_i_add_assign_f((*r).c[2], sn.c[1]);  
        
            in_spinor_field<REAL>(&(sn.c[0]), in, iy-master_shift, 1);  
            in_spinor_field<REAL>(&(sn.c[1]), in, iy-master_shift, 3);  
            _vector_i_add_assign_f(sn.c[0], sn.c[1]);  
        
            _suNf_theta_Z_inverse_multiply(sn.c[1], u, sn.c[0]);  
        
            _vector_add_assign_f((*r).c[1], sn.c[1]);  
            _vector_i_sub_assign_f((*r).c[3], sn.c[1]);
        }
    }
}    

template<typename HSPINOR_TYPE, class REAL, typename SITE_TYPE, typename GAUGE_TYPE>
__global__ void Dphi_gpu_inner_kernel(kernel_field_input* input) {
      SITE_TYPE r;
      SITE_TYPE* field_out = (SITE_TYPE*)input->field_out;
      SITE_TYPE* field_in = (SITE_TYPE*)input->field_in;
      GAUGE_TYPE* gauge = (GAUGE_TYPE*)input->gauge;

      int iy, mu, dir;
      _KERNEL_PIECE_FOR(piece) {
            _IF_IN_BOX_OUT(input, piece) {
                  _spinor_zero_f(r);
                  const int ix = blockIdx.x * BLOCK_SIZE + threadIdx.x + input->base_out[piece-1];

                  _LOOP_DIRECTIONS(ix, iy, mu, dir, input,
                        (evaluate_direction<HSPINOR_TYPE, REAL, GAUGE_TYPE, SITE_TYPE>(&r, field_in, gauge, ix, iy, mu, dir, input->master_shift_in));
                  )

                  _spinor_mul_f(r, -0.5, r);
                  int ix_spinor = ix - input->master_shift_out;
                  write_out_spinor_field<REAL>(&r, field_out, ix_spinor);
            }
      }
}

// Cannot run two boundary kernels at the same time -> race condition
template<typename HSPINOR_TYPE, class REAL, typename SITE_TYPE, typename GAUGE_TYPE>
__global__ void Dphi_gpu_boundary_kernel(kernel_field_input* input) {
      SITE_TYPE r;
      SITE_TYPE res;
      SITE_TYPE* field_out = (SITE_TYPE*)input->field_out;
      SITE_TYPE* field_in = (SITE_TYPE*)input->field_in;
      GAUGE_TYPE *gauge = (GAUGE_TYPE*)input->gauge;

      int iy, ix, mu, dir;
      _KERNEL_PIECE_FOR(piece) {
            _IF_IN_BOX_IN(input, piece) {
                _FIND_BUFFER_DIRECTION(ix, iy, mu, dir, piece, input);
                _spinor_zero_f(r);
                evaluate_direction<HSPINOR_TYPE, REAL, GAUGE_TYPE, SITE_TYPE>(&r, field_in, gauge, ix, iy, mu, dir, input->master_shift_in);
                
                const int ix_spinor = ix - input->master_shift_out;
                in_spinor_field<REAL>(&res, field_out, ix_spinor, 0);
                _spinor_mul_add_assign_f(res, -0.5, r);
                write_out_spinor_field<REAL>(&res, field_out, ix_spinor);
            }
      }
}

#endif