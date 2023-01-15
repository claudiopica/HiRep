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

/* Takes an even input spinor and returns an odd spinor */
/**
 * @brief 
 */
__global__ void Dphi_gpu_inner_kernel(suNf_spinor* __restrict__ out,
                            const suNf_spinor* __restrict__ in,
                            const suNf* __restrict__ gauge_ixp,
                            const suNf* __restrict__ gauge_iyp,
                            const int* __restrict__ iup_d,
                            const int* __restrict__ idn_d,
                            const char* __restrict__ imask_gpu,
                            const int vol4h,
                            const int block_start_ixp, 
                            const int block_start_iyp)
{
  suNf_spinor r;
  _spinor_zero_f(r);

  suNf_hspinor sn;
  suNf u;
  #ifdef FERMION_THETA
    suNf_vector vtmp;
  #endif

  int ix, iy;
  int local_ix, local_iy;
  local_ix = blockIdx.x*BLOCK_SIZE + threadIdx.x;
  if (local_ix < vol4h) {
    int ix = block_start_ixp + local_ix;

    inner_direction(_T_plus,  T_UP_MASK, ix, iup_d[4*ix]);
    inner_direction(_T_minus, T_DN_MASK, ix, idn_d[4*ix]);
    inner_direction(_X_plus,  X_UP_MASK, ix, iup_d[4*ix+1]);
    inner_direction(_X_minus, X_DN_MASK, ix, idn_d[4*ix+1]);
    inner_direction(_Y_plus,  Y_UP_MASK, ix, iup_d[4*ix+2]);
    inner_direction(_Y_minus, Y_DN_MASK, ix, idn_d[4*ix+2]);
    inner_direction(_Z_plus,  Z_UP_MASK, ix, iup_d[4*ix+3]);
    inner_direction(_Z_minus, Z_DN_MASK, ix, idn_d[4*ix+3]);

    _spinor_mul_f(r, -0.5, r);
    write_gpu_suNf_spinor(vol4h, r, out, local_ix, 0);
  }
}

/* Takes an even input spinor and returns an odd spinor */

//Start a kernel for each buffer piece
// use the same kernel as for the bulk calculation
__global__ void Dphi_gpu_boundary_kernel(suNf_spinor* __restrict__ out,
                            const suNf_spinor* __restrict__ in,
                            const suNf* __restrict__ gauge_ixp,
                            const suNf* __restrict__ gauge_iyp,
                            const int* __restrict__ iup_d,
                            const int* __restrict__ idn_d,
                            const char* __restrict__ imask_gpu,
                            const int vol4h, const int buf_stride,
                            const int start, const int start_piece)
{
  suNf_spinor r;
  suNf_spinor res;
  _spinor_zero_f(r);

  int local_ix = blockIdx.x*BLOCK_SIZE + threadIdx.x; 

  if (local_ix < vol4h) {
    int ix = local_ix + start_piece;

    read_gpu_suNf_spinor(vol4h, res, out, local_ix, 0);
    boundary_calculation(_T_plus,  T_UP_MASK, ix, iup_d[4*ix]);
    boundary_calculation(_T_minus, T_DN_MASK, ix, idn_d[4*ix]);
    boundary_calculation(_X_plus,  X_UP_MASK, ix, iup_d[4*ix+1]);
    boundary_calculation(_X_minus, X_DN_MASK, ix, idn_d[4*ix+1]);
    boundary_calculation(_Y_plus,  Y_UP_MASK, ix, iup_d[4*ix+2]);
    boundary_calculation(_Y_minus, Y_DN_MASK, ix, idn_d[4*ix+2]);
    boundary_calculation(_Z_plus,  Z_UP_MASK, ix, iup_d[4*ix+3]);
    boundary_calculation(_Z_minus, Z_DN_MASK, ix, idn_d[4*ix+3]);
    // TODO: Cannot run two of such kernels at the same time because the kernels reference neighbors in different buffers

    _spinor_mul_f(r, -0.5, r);
    _spinor_add_assign_f(res, r);
    write_gpu_suNf_spinor(vol4h, res, out, local_ix, 0);
  }
}



