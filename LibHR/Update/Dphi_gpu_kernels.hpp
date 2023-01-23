/***************************************************************************\
* Copyright (c) 2022, Sofie Martins                                         *   
* All rights reserved.                                                      * 
\***************************************************************************/

#ifndef DPHI_GPU_KERNELS_HPP
#define DPHI_GPU_KERNELS_HPP

#include "./Dphi_gpu_twisted_bc.h"
#include "geometry.h"
#include "libhr_core.h"

__device__ __constant__ int UP_MASK=T_UP_MASK+X_UP_MASK+Y_UP_MASK+Z_UP_MASK;
__device__ __constant__ int DN_MASK=T_DN_MASK+X_DN_MASK+Y_DN_MASK+Z_DN_MASK;
__device__ __constant__ char T_MASK=T_UP_MASK +T_DN_MASK;
__device__ __constant__ char X_MASK=X_UP_MASK +X_DN_MASK;
__device__ __constant__ char Y_MASK=Y_UP_MASK +Y_DN_MASK;
__device__ __constant__ char Z_MASK=Z_UP_MASK +Z_DN_MASK;

#define iup_on_gpu(_dir) int __idx_in_global = iup_d[4*(__idx_out_global) + _dir]
#define idn_on_gpu(_dir) int __idx_in_global = idn_d[4*(__idx_out_global) + _dir]
#define find_neighbor(ix, _dir, _mu) ((_dir == UP) ? iup_d[4*(ix) + _mu] : idn_d[4*(ix) + _mu])
#define MASK(_mu, _dir) (1u << (2*_mu + _dir));

#define _DIR(MASK) ((MASK&UP_MASK) ? UP : DOWN)
#define _MU(MASK) ((MASK&T_MASK) ? 0 : (MASK&X_MASK) ? 1 : (MASK&Y_MASK) ? 2 : 3)

#define c_idx(i) ((i==0) ? 3 : 2)
#define is_c(i) (i == 0 || i == 3)

template<typename VECTOR_TYPE, typename COMPLEX, typename HSPINOR_TYPE, typename SITE_TYPE>
__device__ void read_hspinor_spinmatrix_dev(HSPINOR_TYPE *hspinor, SITE_TYPE *in, int iy, int comp, int mu, int dir) {
   in_spinor_field<COMPLEX>(&hspinor->c[0], in, iy, comp);
   if (is_c(mu)) {
      in_spinor_field<COMPLEX>(&hspinor->c[1], in, iy, comp+2);
   } else { 
      in_spinor_field<COMPLEX>(&hspinor->c[1], in, iy, c_idx(comp));
   }
   //if (iy== 130 && mu == 0 && dir == UP) printf("0+ GPU spinor comp: %0.2e + i%0.2e\n", creal((*hspinor).c[0].c[0]), cimag((*hspinor).c[0].c[0]));
   if (mu % 2 == 1) {
      _vector_i_plus_f((*hspinor).c[1], (*hspinor).c[1]);
   }
   if ((dir + (mu/2)*comp) == 1) {
      _vector_minus_f((*hspinor).c[1], (*hspinor).c[1]);
   }
   _vector_add_assign_f((*hspinor).c[0], (*hspinor).c[1]);
}

template<typename VECTOR_TYPE, typename COMPLEX, typename SITE_TYPE, typename HSPINOR_TYPE>
__device__ void write_hspinor_spinmatrix_dev(SITE_TYPE *spinor_out, HSPINOR_TYPE *hspinor_in, int iy, int comp, int mu, int dir) {
   _vector_add_assign_f((*spinor_out).c[comp], (*hspinor_in).c[1]);
      int cmp2 = comp;
      if ((mu + dir + (mu/2)*comp) % 2 == 1) {
         _vector_minus_f((*hspinor_in).c[1], (*hspinor_in).c[1]);
      }
      if (mu % 2 == 1) {
         _vector_i_plus_f((*hspinor_in).c[1], (*hspinor_in).c[1]);
      }
      if ((mu == 1) || (mu == 2)) {
         cmp2 = !comp;
      }
      _vector_add_assign_f((*spinor_out).c[2+cmp2], (*hspinor_in).c[1]);
}


template<typename HSPINOR_TYPE, typename GAUGE_TYPE>
__device__ void apply_gauge_dev(HSPINOR_TYPE *hspinor, GAUGE_TYPE *u, int mu, int dir) {
      if (mu == 0) {
            if (dir == UP)        { _suNf_theta_T_multiply((*hspinor).c[1], (*u), (*hspinor).c[0]); }
            else if (dir == DOWN) { _suNf_theta_T_inverse_multiply((*hspinor).c[1], (*u), (*hspinor).c[0]); }
      } else if (mu == 1) {
            if (dir == UP)        { _suNf_theta_X_multiply((*hspinor).c[1], (*u), (*hspinor).c[0]); }
            else if (dir == DOWN) { _suNf_theta_X_inverse_multiply((*hspinor).c[1], (*u), (*hspinor).c[0]); }
      } else if (mu == 2) {
            if (dir == UP)        { _suNf_theta_Y_multiply((*hspinor).c[1], (*u), (*hspinor).c[0]); }
            else if (dir == DOWN) { _suNf_theta_Y_inverse_multiply((*hspinor).c[1], (*u), (*hspinor).c[0]); }
      } else if (mu == 3) {
            if (dir == UP)        { _suNf_theta_Z_multiply((*hspinor).c[1], (*u), (*hspinor).c[0]); }
            else if (dir == DOWN) { _suNf_theta_Z_inverse_multiply((*hspinor).c[1], (*u), (*hspinor).c[0]); }
      }
} 

template<typename VECTOR_TYPE, typename COMPLEX, typename SITE_TYPE, typename HSPINOR_TYPE, typename GAUGE_TYPE>
__device__ void core_calculation_dev(int mu, DIRECTION dir, box_t* box, gd_type gd_t,
                                 SITE_TYPE *r, HSPINOR_TYPE *sn, GAUGE_TYPE *u, 
                                 const GAUGE_TYPE *gauge_in, SITE_TYPE* field_in, const int ix, const int iy) {
      //if (iy > 127 && iy < 192)printf("mu: %d, dir: %d, Buffer @ idx: %d\n", mu, dir, iy);
      const int iy_spinor = iy - _master_shift_in(box, gd_t);
      in_gauge_field<COMPLEX>(u, gauge_in, ix, iy, mu, dir);
      //if (iy == 130 && mu == 0 && dir == UP) printf("0+ GPU gauge comp: %0.2e + i%0.2e\n", creal((*u).c[0]), cimag((*u).c[0]));
      for (int hcomp = 0; hcomp <= 1; ++hcomp) {
        read_hspinor_spinmatrix_dev<VECTOR_TYPE, COMPLEX>(sn, field_in, iy_spinor, hcomp, mu, dir);
        apply_gauge_dev(sn, u, mu, dir);
        write_hspinor_spinmatrix_dev<VECTOR_TYPE, COMPLEX>(r, sn, iy_spinor, hcomp, mu, dir);       
      }
}

template<typename VECTOR_TYPE, typename COMPLEX, typename HSPINOR_TYPE, typename SITE_TYPE, typename GAUGE_TYPE>
__device__ void inner_calculation(int mu, DIRECTION dir, 
                                 const char *__restrict__ imask_gpu, 
                                 const int *__restrict__ iup_d, 
                                 const int *__restrict__ idn_d,
                                 box_t* box, gd_type gd_t, gd_type gd_part,
                                 SITE_TYPE *r, HSPINOR_TYPE *sn, GAUGE_TYPE *u, 
                                 const GAUGE_TYPE *__restrict__ gauge_in, 
                                 SITE_TYPE* field_in, int master_shift_in) {
      const char DIR_MASK = MASK(mu, dir);
      const int ix = _idx_out_global(box, gd_part);
      if (imask_gpu[ix]&DIR_MASK) {
            const int iy = find_neighbor(ix, dir, mu);
            core_calculation_dev<VECTOR_TYPE, COMPLEX>(mu, dir, box, gd_t, r, sn, u, gauge_in, field_in, ix, iy);
      }
}

template<typename VECTOR_TYPE, typename COMPLEX, typename HSPINOR_TYPE, typename SITE_TYPE, typename GAUGE_TYPE> 
__device__ void boundary_calculation_dev(int mu, DIRECTION dir, 
                                 const char *__restrict__ imask_gpu, 
                                 const int *__restrict__ iup_d, 
                                 const int *__restrict__ idn_d,
                                 box_t* box_in, box_t* box_out, gd_type gd_t, gd_type gd_part,
                                 SITE_TYPE *r, HSPINOR_TYPE *sn, GAUGE_TYPE *u, 
                                 const GAUGE_TYPE *__restrict__ gauge_in, 
                                 SITE_TYPE* field_in, int master_shift_in) {
        const char DIR_MASK = MASK(mu, dir);
        const int ix = _idx_out_global(box_out, gd_part);
        if (!(imask_gpu[ix]&DIR_MASK)) {
              const int iy = find_neighbor(ix, dir, mu);
              const int iy_loc = iy - _start_in(box_in, gd_part);
              if (iy_loc < _stride_in(box_in, gd_part) && iy_loc >= 0) {
                core_calculation_dev<VECTOR_TYPE, COMPLEX>(mu, dir, box_out, gd_t, r, sn, u, gauge_in, field_in, ix, iy);
              }
        }
}

/* Takes an even input spinor and returns an odd spinor */
/**
 * @brief 
 */
template<typename HSPINOR_TYPE, typename VECTOR_TYPE, typename COMPLEX, typename SITE_TYPE, typename GAUGE_TYPE>
__global__ void Dphi_gpu_inner_kernel(SITE_TYPE* field_in, 
                                        SITE_TYPE* field_out,
                                        const GAUGE_TYPE* __restrict__ gauge,
                                        const int* __restrict__ iup_d,
                                        const int* __restrict__ idn_d,
                                        const char* __restrict__ imask_gpu, 
                                        enum gd_type gd_t, 
                                        box_t* geometry_boxes) {
  SITE_TYPE r;
  HSPINOR_TYPE sn;
  GAUGE_TYPE u;

  box_t* box = geometry_boxes;
  for (int piece = EVEN; piece <= ODD; piece++) {
    if (gd_t & piece) {
      _spinor_zero_f(r);
      if (_idx_out_local < _box_even_volume(box)) {
        for (int mu = 0; mu < 4; ++mu) {
          inner_calculation<VECTOR_TYPE, COMPLEX>(mu, UP, imask_gpu, iup_d, idn_d, box, gd_t, (gd_type)piece, &r, &sn, &u, gauge, field_in, _master_shift_in(box, gd_t));
          inner_calculation<VECTOR_TYPE, COMPLEX>(mu, DOWN, imask_gpu, iup_d, idn_d, box, gd_t, (gd_type)piece, &r, &sn, &u, gauge, field_in, _master_shift_in(box, gd_t));
        }

        _spinor_mul_f(r, -0.5, r);
        int ix = _idx_out_global(box, (gd_type)piece);
        write_out_spinor_field<COMPLEX>(&r, field_out, ix-_master_shift_out(box, gd_t));
      }
    }
  }
}

/* Takes an even input spinor and returns an odd spinor */

//Start a kernel for each buffer piece
// use the same kernel as for the bulk calculation

// Cannot run two boundary kernels at the same time -> race condition

template<typename HSPINOR_TYPE, typename VECTOR_TYPE, typename COMPLEX, typename SITE_TYPE, typename GAUGE_TYPE>
__global__ void Dphi_gpu_boundary_kernel(SITE_TYPE* field_in, 
                                        SITE_TYPE* field_out,
                                        const GAUGE_TYPE* __restrict__ gauge,
                                        const int* __restrict__ iup_d,
                                        const int* __restrict__ idn_d,
                                        const char* __restrict__ imask_gpu, 
                                        enum gd_type gd_t, 
                                        box_t* geometry_boxes,
                                        box_t* box_in) {
    SITE_TYPE r;
    SITE_TYPE res;
    HSPINOR_TYPE sn;
    GAUGE_TYPE u;

    for (int piece = EVEN; piece <= ODD; piece++) {
      if (gd_t & piece) {
        const int iy_loc = blockIdx.x * BLOCK_SIZE + threadIdx.x;
        if (iy_loc < _box_even_volume(box_in)) {
          const int iy = iy_loc + _start_in(box_in, piece);
          const char DIR_MASK = imask_gpu[iy];
          const int mu = _MU(DIR_MASK);
          const int dir_inverted = _DIR(DIR_MASK);
          const int ix = find_neighbor(iy, dir_inverted, mu);
          const DIRECTION dir = (DIRECTION)!dir_inverted;
          _spinor_zero_f(r);
          core_calculation_dev<VECTOR_TYPE, COMPLEX>(mu, dir, geometry_boxes, gd_t, &r, &sn, &u, gauge, field_in, ix, iy);
          
          in_spinor_field<COMPLEX>(&res, field_out, ix-_master_shift_out(geometry_boxes, gd_t), 0);
          _spinor_mul_add_assign_f(res, -0.5, r);
          write_out_spinor_field<COMPLEX>(&res, field_out, ix-_master_shift_out(geometry_boxes, gd_t));
        }
      }
    }
}

#endif