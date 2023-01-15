/***************************************************************************\
* Copyright (c) 2022, Sofie Martins                                         *   
* All rights reserved.                                                      * 
\***************************************************************************/

#include "spinor_field.h"
#include "global_gpu.h"
#include "Geometry/gpu_geometry.h"
#include "suN.h"

#include "Dphi_gpu_twisted_bc.h"
#include "Dphi_gpu_directions_flt.h"

typedef struct suNf_hspinor_flt
{
  suNf_vector_flt c[2];
} suNf_hspinor_flt;

__global__ void Dphi_gpu_inner_kernel_flt(suNf_spinor_flt* __restrict__ out, 
                                        const suNf_spinor_flt* __restrict__ in, 
                                        const suNf_flt* __restrict__ gauge_ixp, 
                                        const suNf_flt* __restrict__ gauge_iyp, 
                                        const int* __restrict__ iup_d, 
                                        const int* __restrict__ idn_d, 
                                        const char* __restrict__ imask_gpu, 
                                        const int vol4h, 
                                        const int block_start_ixp, 
                                        const int block_start_iyp) 
{
  suNf_spinor_flt r;
  _spinor_zero_f(r);

  int local_ix = blockIdx.x*BLOCK_SIZE + threadIdx.x;
  if (local_ix < vol4h) {
    int ix = block_start_ixp + local_ix;

    inner_direction_flt(_T_plus_flt,  T_UP_MASK, ix, iup_d[4*ix]);
    inner_direction_flt(_T_minus_flt, T_DN_MASK, ix, idn_d[4*ix]);
    inner_direction_flt(_X_plus_flt,  X_UP_MASK, ix, iup_d[4*ix+1]);
    inner_direction_flt(_X_minus_flt, X_DN_MASK, ix, idn_d[4*ix+1]);
    inner_direction_flt(_Y_plus_flt,  Y_UP_MASK, ix, iup_d[4*ix+2]);
    inner_direction_flt(_Y_minus_flt, Y_DN_MASK, ix, idn_d[4*ix+2]);
    inner_direction_flt(_Z_plus_flt,  Z_UP_MASK, ix, iup_d[4*ix+3]);
    inner_direction_flt(_Z_minus_flt, Z_DN_MASK, ix, idn_d[4*ix+3]);

    _spinor_mul_f(r, -0.5, r);
    write_gpu_suNf_spinor_flt(vol4h, r, out, local_ix, 0);
  }
}

__global__ void Dphi_gpu_boundary_kernel_flt(suNf_spinor_flt* __restrict__ out,
                            const suNf_spinor_flt* __restrict__ in,
                            const suNf_flt* __restrict__ gauge_ixp,
                            const suNf_flt* __restrict__ gauge_iyp,
                            const int* __restrict__ iup_d,
                            const int* __restrict__ idn_d,
                            const char* __restrict__ imask_gpu,
                            const int vol4h, const int buf_stride,
                            const int start, const int start_piece)
{
  suNf_spinor_flt r;
  suNf_spinor_flt res;
  _spinor_zero_f(r);

  suNf_hspinor_flt sn;
  suNf_flt u;
  #ifdef FERMION_THETA
    suNf_vector vtmp;
  #endif

  int local_ix = blockIdx.x*BLOCK_SIZE + threadIdx.x; 
  if (local_ix < vol4h) {
    int ix = local_ix + start_piece;

    read_gpu_suNf_spinor_flt(vol4h, res, out, local_ix, 0);
    __syncthreads();
    boundary_calculation_flt(_T_plus_flt,  T_UP_MASK, ix, iup_d[4*ix]);
    boundary_calculation_flt(_T_minus_flt, T_DN_MASK, ix, idn_d[4*ix]);
    boundary_calculation_flt(_X_plus_flt,  X_UP_MASK, ix, iup_d[4*ix+1]);
    boundary_calculation_flt(_X_minus_flt, X_DN_MASK, ix, idn_d[4*ix+1]);
    boundary_calculation_flt(_Y_plus_flt,  Y_UP_MASK, ix, iup_d[4*ix+2]);
    boundary_calculation_flt(_Y_minus_flt, Y_DN_MASK, ix, idn_d[4*ix+2]);
    boundary_calculation_flt(_Z_plus_flt,  Z_UP_MASK, ix, iup_d[4*ix+3]);
    boundary_calculation_flt(_Z_minus_flt, Z_DN_MASK, ix, idn_d[4*ix+3]);
    // TODO: Cannot run two of such kernels at the same time because the kernels reference neighbors in different buffers

    _spinor_mul_f(r, -0.5, r);
    _spinor_add_assign_f(res, r);
    write_gpu_suNf_spinor_flt(vol4h, res, out, local_ix, 0);
  }
}