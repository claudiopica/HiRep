/***************************************************************************\
* Copyright (c) 2022, Sofie Martins                                         *   
* All rights reserved.                                                      * 
\***************************************************************************/

#include "spinor_field.h"
#include "global_gpu.h"
#include "Geometry/gpu_geometry.h"
#include "suN.h"

#include "Dphi_gpu_twisted_bc.h"
#include "Dphi_gpu_directions.h"

typedef struct suNf_hspinor
{
  suNf_vector c[2];
} suNf_hspinor;

#define iup_on_gpu(_dir) int __idx_in_global = iup_d[4*(__idx_out_global) + _dir]
#define idn_on_gpu(_dir) int __idx_in_global = idn_d[4*(__idx_out_global) + _dir]

/* Takes an even input spinor and returns an odd spinor */
/**
 * @brief 
 */
__global__ void Dphi_gpu_inner_kernel(kernel_field_input* input_even, kernel_field_input* input_odd, 
                                        const suNf* __restrict__ gauge,
                                        const int* __restrict__ iup_d,
                                        const int* __restrict__ idn_d,
                                        const char* __restrict__ imask_gpu, 
                                        enum gd_type gd_t) {
  suNf_spinor r;
  _spinor_zero_f(r);
  suNf_hspinor sn;
  suNf u;

  if (gd_t & EVEN) do {
      _KERNEL_FOR_WITH_GAUGE(input_even, gauge, suNf_spinor, suNf) {
      inner_direction(_T_plus,  T_UP_MASK, iup_on_gpu(0));
      inner_direction(_T_minus, T_DN_MASK, idn_on_gpu(0));
      inner_direction(_X_plus,  X_UP_MASK, iup_on_gpu(1));
      inner_direction(_X_minus, X_DN_MASK, idn_on_gpu(1));
      inner_direction(_Y_plus,  Y_UP_MASK, iup_on_gpu(2));
      inner_direction(_Y_minus, Y_DN_MASK, idn_on_gpu(2));
      inner_direction(_Z_plus,  Z_UP_MASK, iup_on_gpu(3));
      inner_direction(_Z_minus, Z_DN_MASK, idn_on_gpu(3));

      _spinor_mul_f(r, -0.5, r);
      _WRITE_OUT_SPINOR_FIELD(r);

    }
  } while (0);

  //__syncthreads();
  _spinor_zero_f(r);

  if (gd_t & ODD) do {
    _KERNEL_FOR_WITH_GAUGE(input_odd, gauge, suNf_spinor, suNf) {
      inner_direction(_T_plus,  T_UP_MASK, iup_on_gpu(0));
      inner_direction(_T_minus, T_DN_MASK, idn_on_gpu(0));
      inner_direction(_X_plus,  X_UP_MASK, iup_on_gpu(1));
      inner_direction(_X_minus, X_DN_MASK, idn_on_gpu(1));
      inner_direction(_Y_plus,  Y_UP_MASK, iup_on_gpu(2));
      inner_direction(_Y_minus, Y_DN_MASK, idn_on_gpu(2));
      inner_direction(_Z_plus,  Z_UP_MASK, iup_on_gpu(3));
      inner_direction(_Z_minus, Z_DN_MASK, idn_on_gpu(3));

      //if (__idx_out_global==0) printf("GPU spinor %0.2e + i%0.2e\n", creal(r.c[0].c[0]), cimag(r.c[0].c[0]));
      _spinor_mul_f(r, -0.5, r);
      _WRITE_OUT_SPINOR_FIELD(r);
    } 
  } while (0);

}

/* Takes an even input spinor and returns an odd spinor */

//Start a kernel for each buffer piece
// use the same kernel as for the bulk calculation

// Cannot run two boundary kernels at the same time -> race condition
__global__ void Dphi_gpu_boundary_kernel(kernel_field_input* input_even, 
                                        kernel_field_input* input_odd,
                                        const suNf* __restrict__ gauge, 
                                        const int* __restrict__ iup_d, 
                                        const int* __restrict__ idn_d, 
                                        const char* __restrict__ imask_gpu, 
                                        enum gd_type gd_t) {
    suNf_spinor r;
    _spinor_zero_f(r);
    suNf_spinor res;
    suNf_hspinor sn;
    suNf u;
    
    if (gd_t & EVEN) do {
      _KERNEL_FOR_WITH_GAUGE(input_even, gauge, suNf_spinor, suNf) {
        _OUT_SPINOR_FIELD(res);
        boundary_calculation(_T_plus,  T_UP_MASK, iup_on_gpu(0));
        boundary_calculation(_T_minus, T_DN_MASK, idn_on_gpu(0));
        boundary_calculation(_X_plus,  X_UP_MASK, iup_on_gpu(1));
        boundary_calculation(_X_minus, X_DN_MASK, idn_on_gpu(1));
        boundary_calculation(_Y_plus,  Y_UP_MASK, iup_on_gpu(2));
        boundary_calculation(_Y_minus, Y_DN_MASK, idn_on_gpu(2));
        boundary_calculation(_Z_plus,  Z_UP_MASK, iup_on_gpu(3));
        boundary_calculation(_Z_minus, Z_DN_MASK, idn_on_gpu(3));

        _spinor_mul_f(r, -0.5, r);
        _spinor_add_assign_f(res, r);
        _WRITE_OUT_SPINOR_FIELD(res);  
      }
    } while (0);

    _spinor_zero_f(r);
    
    if (gd_t & ODD) do {
      _KERNEL_FOR_WITH_GAUGE(input_odd, gauge, suNf_spinor, suNf) {
        _OUT_SPINOR_FIELD(res);
        boundary_calculation(_T_plus,  T_UP_MASK, iup_on_gpu(0));
        boundary_calculation(_T_minus, T_DN_MASK, idn_on_gpu(0));
        boundary_calculation(_X_plus,  X_UP_MASK, iup_on_gpu(1));
        boundary_calculation(_X_minus, X_DN_MASK, idn_on_gpu(1));
        boundary_calculation(_Y_plus,  Y_UP_MASK, iup_on_gpu(2));
        boundary_calculation(_Y_minus, Y_DN_MASK, idn_on_gpu(2));
        boundary_calculation(_Z_plus,  Z_UP_MASK, iup_on_gpu(3));
        boundary_calculation(_Z_minus, Z_DN_MASK, idn_on_gpu(3));

        _spinor_mul_f(r, -0.5, r);
        _spinor_add_assign_f(res, r);
        _WRITE_OUT_SPINOR_FIELD(res);  
      }
    } while (0);
}