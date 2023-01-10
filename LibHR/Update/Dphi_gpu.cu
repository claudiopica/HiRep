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

#include "update.h"
#include "libhr_core.h"
#include "Inverters/linear_algebra.h"
#include "error.h"
#include "io.h"
#include "memory.h"
#include "utils.h"

#include "./Dphi_gpu_kernels.hpp"

#ifdef ROTATED_SF
  extern rhmc_par _update_par; /* Update/update_rhmc.c */
#endif 

static int init=1;
static spinor_field *gtmp=NULL;
static spinor_field *etmp=NULL;
static spinor_field *otmp=NULL;

/**
 * @brief Initializes the boundary conditions for fermion twisting
 */
static void init_bc_gpu(){
  #ifdef FERMION_THETA
  static int initialized=0;
  if (!initialized){
    cudaMemcpyToSymbol(eitheta_gpu, eitheta, 4*sizeof(hr_complex), 0, cudaMemcpyHostToDevice);
    CudaCheckError();
    initialized=1;
  }
  #endif
}

/**
 * @brief the following variable is used to keep trace of
 *        matrix-vector multiplication.
 *        we count how many time the function Dphi_ is called
 */
static unsigned long int MVMcounter=0;

/**
 * @brief Getter for number of applications of the Dirac operator
 */
unsigned long int getMVM_gpu() {
	unsigned long int res=MVMcounter>>1; /* divide by two */
	return res;
}

/**
 * @brief Reset counter for number of applications of the Dirac operator
 */
void resetMVM_gpu() {
  MVMcounter = 0;
}

/**
 * @brief Free fields allocated for intermediate storage of field data.
 */
static void free_mem() {
  if (gtmp!=NULL) {
    free_spinor_field_f(gtmp);
    etmp=NULL;
  }
  if (etmp!=NULL) {
    free_spinor_field_f(etmp);
    etmp=NULL;
  }
  if (otmp!=NULL) {
    free_spinor_field_f(otmp);
    otmp=NULL;
  }
  init=1;
}

/**
 * @brief Allocate fields intended for storage of field data in intermediate
 *        steps
 */
static void init_Dirac() {
  if (init) {
    alloc_mem_t=GPU_MEM;

    gtmp=alloc_spinor_field_f(1, &glattice);
    etmp=alloc_spinor_field_f(1, &glat_even);
    otmp=alloc_spinor_field_f(1, &glat_odd);

    alloc_mem_t=std_mem_t;

    atexit(&free_mem);
    init=0;
  }
}

/**
 * @brief Applies the Wilson-Dirac operator only where the calculation
 *        does not depend on the communications. Use this calculation
 *        to mask communication latency.
 *
 * @param in                  Input spinor field, defined on GPU copy
 * @param out                 Output spinor field to save result
 */
static void Dphi_inner_gpu_(spinor_field *out, spinor_field *in) {
  int grid, iyp, gauge_ixp, gauge_iyp, write_stride, start_piece_ixp, start_piece_iyp;
  _PIECE_FOR(out->type, ixp)
  {
      if (in->type==&glattice) iyp = (ixp+1)%2;
      else iyp = 0;

      if (in->type==&glattice) {
        gauge_ixp = ixp;
        gauge_iyp = iyp;
      } else if (in->type==&glat_odd) {
        gauge_ixp = 0;
        gauge_iyp = 1;
      } else if (in->type==&glat_even) {
        gauge_ixp = 1;
        gauge_iyp = 0;
      }

      write_stride = out->type->master_end[ixp] - out->type->master_start[ixp] + 1;
      grid = (write_stride-1)/BLOCK_SIZE + 1;
      start_piece_ixp = out->type->master_start[ixp];
      start_piece_iyp = in->type->master_start[iyp];

      Dphi_gpu_inner_kernel<<<grid, BLOCK_SIZE>>>(
            _GPU_FIELD_BLK(out, ixp), _GPU_FIELD_BLK(in, iyp),
            _GPU_4FIELD_BLK(u_gauge_f, gauge_ixp), _GPU_4FIELD_BLK(u_gauge_f, gauge_iyp),
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
static void Dphi_boundary_gpu_(spinor_field *out, spinor_field *in) {
  #if defined(WITH_MPI) && defined(WITH_NEW_GEOMETRY)
  int grid, ixp, gauge_ixp, gauge_i, block_stride, buffer_stride, block_start, buffer_start;

  for (int i = 0; i < in->type->nbuffers_spinor; ++i)
  {
      if (out->type==&glat_even) {
        ixp = 0;
        gauge_ixp = 0;
        gauge_i = 2*i+1;
      } else if (out->type==&glat_odd) {
        ixp = 0;
        gauge_ixp = 1;
        gauge_i = 2*i;
      } else {
        ixp = (i%2==0);
        gauge_ixp = ixp;
        gauge_i = i;
      }

      block_stride = out->type->master_end[ixp] - out->type->master_start[ixp] + 1;
      block_start = out->type->master_start[ixp];

      buffer_stride = in->type->rbuf_len[i];// Different strides for gauge and spinor
      buffer_start = in->type->rbuf_start[i];
      grid = (block_stride-1)/BLOCK_SIZE + 1;

      Dphi_gpu_boundary_kernel<<<grid, BLOCK_SIZE>>>(
            _GPU_FIELD_BLK(out, ixp), _BUF_GPU_FIELD_BLK(in, i),
            _GPU_4FIELD_BLK(u_gauge_f, gauge_ixp), _BUF_GPU_4FIELD_BLK(u_gauge_f, gauge_i),
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

  //init_bc_gpu();

  #ifdef WITH_MPI
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
unsigned long int (*getMVM) ()=getMVM_gpu;
void (*Dphi_) (spinor_field *out, spinor_field *in)=Dphi_gpu_;
void (*Dphi) (double m0, spinor_field *out, spinor_field *in)=Dphi_gpu;
void (*g5Dphi) (double m0, spinor_field *out, spinor_field *in)=g5Dphi_gpu;
void (*g5Dphi_sq) (double m0, spinor_field *out, spinor_field *in)=g5Dphi_sq_gpu;
void (*Dphi_eopre) (double m0, spinor_field *out, spinor_field *in)=Dphi_eopre_gpu;
void (*Dphi_oepre) (double m0, spinor_field *out, spinor_field *in)=Dphi_oepre_gpu;
void (*g5Dphi_eopre) (double m0, spinor_field *out, spinor_field *in)=g5Dphi_eopre_gpu;
void (*g5Dphi_eopre_sq) (double m0, spinor_field *out, spinor_field *in)=g5Dphi_eopre_sq_gpu;

#endif
