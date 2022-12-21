/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *
* All rights reserved.                                                      *
\***************************************************************************/

/**
 * @file Dphi_gpu.c
 * @brief GPU implementation of the Wilson-Dirac operator D and its hermitian version
 *        on a given double-precision spinor field.
 */

#ifdef WITH_GPU

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "suN.h"
#include "global.h"
#include "error.h"
#include "dirac.h"
#include "linear_algebra.h"
#include "spinor_field.h"
#include "geometry.h"
#include "communications.h"
#include "memory.h"
#include "gpu.h"
#include "hr_complex.h"

#include "Dphi_gpu_kernels.cu"
#include "Dirac_utils.c"

#ifdef __cplusplus
  extern "C" {
#endif

/**
 * @brief Applies the Wilson-Dirac operator only where the calculation
 *        does not depend on the communications. Use this calculation
 *        to mask communication latency.
 *
 * @param in                  Input spinor field, defined on GPU copy
 * @param out                 Output spinor field to save result
 */
static void Dphi_inner_gpu_(spinor_field *out, spinor_field *in) {
  int grid, iyp, write_stride, start_piece_ixp, start_piece_iyp;
  _PIECE_FOR(out->type, ixp)
  {
      iyp = (ixp+1)%2;
      write_stride = out->type->master_end[ixp] - out->type->master_start[ixp] + 1;
      grid = (write_stride-1)/BLOCK_SIZE + 1;
      start_piece_ixp = in->type->master_start[ixp];
      start_piece_iyp = in->type->master_start[iyp];

      Dphi_gpu_inner_kernel<<<grid, BLOCK_SIZE>>>(
            _GPU_FIELD_BLK(out, ixp), _GPU_FIELD_BLK(in, iyp),
            _GPU_4FIELD_BLK(u_gauge_f, ixp), _GPU_4FIELD_BLK(u_gauge_f, iyp),
            iup_gpu, idn_gpu, imask_gpu,
            write_stride, start_piece_ixp, start_piece_iyp);
      CudaCheckError();
  }
  cudaDeviceSynchronize();
}

/**
 * @brief Applies the Wilson-Dirac operator only where the calculation
 *        depends on the communications. Call this after having
 *        completed the spinor field communications.
 *
 * @param in                  Input spinor field defined on the extended
 *                            lattice on the GPU copy.
 * @param out                 Output spinor field to save result
 */
void Dphi_boundary_gpu_(spinor_field *out, spinor_field *in) {
  #if defined(WITH_MPI) && defined(WITH_NEW_GEOMETRY)
  int grid, ixp, block_stride, buffer_stride, block_start, buffer_start;
  for (int i = 0; i < in->type->nbuffers_spinor; ++i)
  {
      ixp = (i%2==0);
      block_stride = in->type->master_end[ixp] - in->type->master_start[ixp] + 1;
      buffer_stride = in->type->rbuf_len[i];
      buffer_start = in->type->rbuf_start[i];
      grid = (block_stride-1)/BLOCK_SIZE + 1;
      block_start = in->type->master_start[ixp];

      Dphi_gpu_boundary_kernel<<<grid, BLOCK_SIZE>>>(
            _GPU_FIELD_BLK(out, ixp), _BUF_GPU_FIELD_BLK(in, i),
            _GPU_4FIELD_BLK(u_gauge_f, ixp), _BUF_GPU_4FIELD_BLK(u_gauge_f, i),
            iup_gpu, idn_gpu, imask_gpu,
            block_stride, buffer_stride,
            buffer_start, block_start);
      CudaCheckError();
  }
  #endif
}

/**
 * @brief Implementation of the massless Wilson-Dirac operator on the GPU
 *
 * @param in                      Input spinor field
 * @param out                     Output spinor field to save result
 */
void Dphi_gpu_(spinor_field *out, spinor_field *in)
{
  // Input parameter validity checks
  _CHECK_GEOMETRY_EO(in, out);
  // TODO: Possibly add: Fields are not null, fields are unequal (SAM)

  init_bc_gpu();

  #ifdef WITH_MPI
      //TODO: gfield comms after copy, not necessary for every Dphi_gpu_ (SAM)
    start_sendrecv_gpu_gfield_f(u_gauge_f);
    complete_sendrecv_gpu_gfield_f(u_gauge_f);

    start_sendrecv_gpu_spinor_field_f(in);
  #endif

  // Mask communication latency with calculating the points
  // that do no rely on communications  
  Dphi_inner_gpu_(out, in);
  complete_sendrecv_gpu_spinor_field_f(in);

  // After completion of communications add missing calculations
  // that relied on communications.
  Dphi_boundary_gpu_(out, in);
}

/**
 * @brief Implementation of the Wilson-Dirac operator with mass on the GPU.
 *
 * @param in                    Input spinor field defined on the full lattice
 * @param out                   Output spinor field defined on the full lattice that is 
 *                              the result of the application of the Wilson-Dirac operator 
 *                              on the input spinor field.
 * @param m0                    Mass parameter
 */
void Dphi_gpu(double m0, spinor_field *out, spinor_field *in)
{
  //Input argument validity checks
  _CHECK_GEOMETRY_FULL(in);
  _CHECK_GEOMETRY_FULL(out);

  // Operation
  double rho = 4. + m0;
  Dphi_gpu_(out, in);
  spinor_field_mul_add_assign_f_gpu(out, rho, in);
}

/**
 * @brief Implementation of the Hermitian Wilson-Dirac operator with mass on the GPU.
 *
 * @param in                    Input spinor field defined on the full lattice
 * @param out                   Output spinor field defined on the full lattice that is 
 *                              the result of the application of the Hermitian Wilson-Dirac 
 *                              operator on the input spinor field
 * @param m0                    Mass parameter
 */
void g5Dphi_gpu(double m0, spinor_field *out, spinor_field *in)
{
  // Input argument validity checks
  _CHECK_GEOMETRY_FULL(in);
  _CHECK_GEOMETRY_FULL(out);

  // Operation
  double rho = 4. + m0;
  Dphi_gpu_(out, in);
  spinor_field_mul_add_assign_f_gpu(out, rho, in);
  spinor_field_g5_assign_f_gpu(out);
}

 /**
  * @brief Even-odd preconditioned Wilson-Dirac operator with mass on the GPU.
  *
  * @param in                 Even input spinor field
  * @param out                Even output spinor field that is the result of the application
  *                           of the even-odd preconditioned Wilson-Dirac operator on the 
  *                           input spinor field
  * @param m0                 Mass parameter
  */
void Dphi_eopre_gpu(double m0, spinor_field *out, spinor_field *in)
{
  // Input argument validity checks
  _CHECK_GEOMETRY_EVEN(in);
  _CHECK_GEOMETRY_EVEN(out);

  // alloc memory for temporary spinor field
  if (init) { init_Dirac(); }

  // Operation
  Dphi_gpu_(otmp, in);
  Dphi_gpu_(out, otmp);
  double rho = 4. + m0;
  rho*=-rho; /* this minus sign is taken into account below */
  spinor_field_mul_add_assign_f_gpu(out, rho, in);
  spinor_field_minus_f_gpu(out, out);
}

 /**
  * @brief Even-odd preconditioned Wilson-Dirac operator with mass on the GPU.
  *
  * @param in                 Odd input spinor field
  * @param out                Odd output spinor field that is the result of the application
  *                           of the even-odd preconditioned Wilson-Dirac operator on the 
  *                           input spinor field
  * @param m0                 Mass parameter
  */
void Dphi_oepre_gpu(double m0, spinor_field *out, spinor_field *in)
{
  // Input argument validity checks
  _CHECK_GEOMETRY_ODD(in);
  _CHECK_GEOMETRY_ODD(out);

  // alloc memory for temporary spinor field 
  if (init) { init_Dirac();}

  Dphi_gpu_(etmp, in);
  Dphi_gpu_(out, etmp);
  double rho = 4. + m0;
  rho*=-rho; /* this minus sign is taken into account below */
  spinor_field_mul_add_assign_f_gpu(out, rho, in);
  spinor_field_minus_f_gpu(out, out);
}

 /**
  * @brief Even-odd preconditioned Hermitian Wilson-Dirac operator with mass on the GPU.
  *
  * @param in                 Even input spinor field
  * @param out                Even output spinor field that is the result of the application
  *                           of the even-odd preconditioned Wilson-Dirac operator on the 
  *                           input spinor field
  * @param m0                 Mass parameter
  */
void g5Dphi_eopre_gpu(double m0, spinor_field *out, spinor_field *in)
{
  // Input argument validity checks
  _CHECK_GEOMETRY_EVEN(in);
  _CHECK_GEOMETRY_EVEN(out);

  #if defined(BASIC_SF) || defined(ROTATED_SF)
    SF_spinor_bcs(in);
  #endif

  // alloc memory for temporary spinor field 
  if (init) { init_Dirac();}

  // Operation
  Dphi_gpu_(otmp, in);
  Dphi_gpu_(out, otmp);
  double rho = 4. + m0;
  rho*=-rho; /* this minus sign is taken into account below */
  spinor_field_mul_add_assign_f_gpu(out, rho, in);
  spinor_field_minus_f_gpu(out, out);
  spinor_field_g5_assign_f_gpu(out);

  #if defined(BASIC_SF) || defined(ROTATED_SF)
    SF_spinor_bcs(out);
  #endif 
}

 /**
  * @brief Even-odd preconditioned squared Hermitian Wilson-Dirac operator with mass on the GPU.
  *
  * @param in                 Even input spinor field
  * @param out                Even output spinor field that is the result of the application
  *                           of the even-odd preconditioned Wilson-Dirac operator on the 
  *                           input spinor field
  * @param m0                 Mass parameter
  */
void g5Dphi_eopre_sq_gpu(double m0, spinor_field *out, spinor_field *in) {
  // alloc memory for temporary spinor field
  if (init) { init_Dirac(); }

  g5Dphi_eopre_gpu(m0, etmp, in);
  g5Dphi_eopre_gpu(m0, out, etmp);
}

/**
 * @brief Implementation of the squared Hermitian Wilson-Dirac operator with mass on the GPU.
 *
 * @param in                    Input spinor field defined on the full lattice
 * @param out                   Output spinor field defined on the full lattice that is 
 *                              the result of the application of the Hermitian Wilson-Dirac 
 *                              operator on the input spinor field
 * @param m0                    Mass parameter
 */
void g5Dphi_sq_gpu(double m0, spinor_field *out, spinor_field *in) {
  // alloc memory for temporary spinor field 
  if (init) { init_Dirac();  }

  g5Dphi_gpu(m0, gtmp, in);
  g5Dphi_gpu(m0, out, gtmp);
}

/* For WITH_GPU: Map the GPU functions to the default functions. */
void (*Dphi_) (spinor_field *out, spinor_field *in)=Dphi_gpu_;
void (*Dphi) (double m0, spinor_field *out, spinor_field *in)=Dphi_gpu;
void (*g5Dphi) (double m0, spinor_field *out, spinor_field *in)=g5Dphi_gpu;
void (*g5Dphi_sq) (double m0, spinor_field *out, spinor_field *in)=g5Dphi_sq_gpu;
unsigned long int (*getMVM) ()=getMVM_gpu;
void (*Dphi_eopre) (double m0, spinor_field *out, spinor_field *in)=Dphi_eopre_gpu;
void (*Dphi_oepre) (double m0, spinor_field *out, spinor_field *in)=Dphi_oepre_gpu;
void (*g5Dphi_eopre) (double m0, spinor_field *out, spinor_field *in)=g5Dphi_eopre_gpu;
void (*g5Dphi_eopre_sq) (double m0, spinor_field *out, spinor_field *in)=g5Dphi_eopre_sq_gpu;

#ifdef __cplusplus
  }
#endif

#endif
