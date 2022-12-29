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

#include "./Dphi_gpu_kernels_flt.hpp"

#ifdef ROTATED_SF
  extern rhmc_par _update_par; /* Update/update_rhmc.c */
#endif 

static int init=1;
static spinor_field_flt *gtmp=NULL;
static spinor_field_flt *etmp=NULL;
static spinor_field_flt *otmp=NULL;

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
unsigned long int getMVM_flt_gpu() {
	unsigned long int res=MVMcounter>>1; /* divide by two */
	return res;
}

/**
 * @brief Reset counter for number of applications of the Dirac operator
 */
void resetMVM_flt_gpu() {
  MVMcounter = 0;
}

/**
 * @brief Free fields allocated for intermediate storage of field data.
 */
static void free_mem() {
  if (gtmp!=NULL) {
    free_spinor_field_f_flt(gtmp);
    etmp=NULL;
  }
  if (etmp!=NULL) {
    free_spinor_field_f_flt(etmp);
    etmp=NULL;
  }
  if (otmp!=NULL) {
    free_spinor_field_f_flt(otmp);
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

    gtmp=alloc_spinor_field_f_flt(1, &glattice);
    etmp=alloc_spinor_field_f_flt(1, &glat_even);
    otmp=alloc_spinor_field_f_flt(1, &glat_odd);

    alloc_mem_t=std_mem_t;

    atexit(&free_mem);
    init=0;
  }
}

/**
 * @brief Applies the Wilson-Dirac operator only where the calculation
 *        does not depend on the communications. Use this calculation
 *        to mask communication latency. Single precision version.
 *
 * @param in                  Input spinor field, defined on GPU copy
 * @param out                 Output spinor field to save result
 */
static void Dphi_inner_gpu_flt_(spinor_field_flt *out, spinor_field_flt *in) {
  int grid, iyp, write_stride, start_piece_ixp, start_piece_iyp;
  _PIECE_FOR(out->type, ixp)
  {
      iyp = (ixp+1)%2;
      write_stride = out->type->master_end[ixp] - out->type->master_start[ixp] + 1;
      grid = (write_stride-1)/BLOCK_SIZE + 1;
      start_piece_ixp = in->type->master_start[ixp];
      start_piece_iyp = in->type->master_start[iyp];

      Dphi_gpu_inner_kernel_flt<<<grid, BLOCK_SIZE>>>(
            _GPU_FIELD_BLK(out, ixp), _GPU_FIELD_BLK(in, iyp),
            _GPU_4FIELD_BLK(u_gauge_f_flt, ixp), _GPU_4FIELD_BLK(u_gauge_f_flt, iyp),
            iup_gpu, idn_gpu, imask_gpu,
            write_stride, start_piece_ixp, start_piece_iyp);
      CudaCheckError();
  }
  cudaDeviceSynchronize();
}

/**
 * @brief Applies the Wilson-Dirac operator only where the calculation
 *        depends on the communications. Call this after having
 *        completed the spinor field communications. Single precision
 *        version.
 *
 * @param in                  Input spinor field defined on the extended
 *                            lattice on the GPU copy.
 * @param out                 Output spinor field to save result
 */
static void Dphi_boundary_gpu_flt_(spinor_field_flt *out, spinor_field_flt *in) {
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

      Dphi_gpu_boundary_kernel_flt<<<grid, BLOCK_SIZE>>>(
            _GPU_FIELD_BLK(out, ixp), _BUF_GPU_FIELD_BLK(in, i),
            _GPU_4FIELD_BLK(u_gauge_f_flt, ixp), _BUF_GPU_4FIELD_BLK(u_gauge_f_flt, i),
            iup_gpu, idn_gpu, imask_gpu,
            block_stride, buffer_stride,
            buffer_start, block_start);
      CudaCheckError();
  }
  #endif
}

/**
 * @brief Implementation of the massless Wilson-Dirac operator on the GPU
 *        Single precision version
 *
 * @param in                      Input spinor field
 * @param out                     Output spinor field to save result
 */
void Dphi_flt_gpu_(spinor_field_flt *out, spinor_field_flt *in)
{
  // Input parameter validity checks
  _CHECK_GEOMETRY_EO(in, out);
  // TODO: Possibly add: Fields are not null, fields are unequal (SAM)

  init_bc_gpu();

  #ifdef WITH_MPI
    start_sendrecv_gpu_spinor_field_f_flt(in);
  #endif

  // Mask communication latency with calculating the points
  // that do no rely on communications  
  Dphi_inner_gpu_flt_(out, in);
  complete_sendrecv_gpu_spinor_field_f_flt(in);

  // After completion of communications add missing calculations
  // that relied on communications.
  Dphi_boundary_gpu_flt_(out, in);
}

/**
 * @brief Implementation of the Wilson-Dirac operator with mass on the GPU.
 *        Single precision version.
 *
 * @param in                    Input spinor field defined on the full lattice
 * @param out                   Output spinor field defined on the full lattice that is 
 *                              the result of the application of the Wilson-Dirac operator 
 *                              on the input spinor field.
 * @param m0                    Mass parameter
 */
void Dphi_flt_gpu(double m0, spinor_field_flt *out, spinor_field_flt *in)
{
  //Input argument validity checks
  _CHECK_GEOMETRY_FULL(in);
  _CHECK_GEOMETRY_FULL(out);

  // Operation
  float rho = 4.f + (float)m0;
  Dphi_flt_gpu_(out, in);
  spinor_field_mul_add_assign_f_flt_gpu(out, rho, in);
}

/**
 * @brief Implementation of the Hermitian Wilson-Dirac operator with mass on the GPU.
 *        Single precision version.
 *
 * @param in                    Input spinor field defined on the full lattice
 * @param out                   Output spinor field defined on the full lattice that is 
 *                              the result of the application of the Hermitian Wilson-Dirac 
 *                              operator on the input spinor field
 * @param m0                    Mass parameter
 */
void g5Dphi_flt_gpu(double m0, spinor_field_flt *out, spinor_field_flt *in)
{
  // Input argument validity checks
  _CHECK_GEOMETRY_FULL(in);
  _CHECK_GEOMETRY_FULL(out);

  // Operation
  apply_BCs_on_spinor_field_flt(in);
  float rho = 4.f + (float)m0;
  Dphi_flt_gpu_(out, in);
  spinor_field_mul_add_assign_f_flt_gpu(out, rho, in);
  spinor_field_g5_assign_f_flt_gpu(out);
  apply_BCs_on_spinor_field_flt(out);
}

 /**
  * @brief Even-odd preconditioned Wilson-Dirac operator with mass on the GPU.
  *         Single precision version.
  *
  * @param in                 Even input spinor field
  * @param out                Even output spinor field that is the result of the application
  *                           of the even-odd preconditioned Wilson-Dirac operator on the 
  *                           input spinor field
  * @param m0                 Mass parameter
  */
void Dphi_eopre_flt_gpu(double m0, spinor_field_flt *out, spinor_field_flt *in)
{
  // Input argument validity checks
  _CHECK_GEOMETRY_EVEN(in);
  _CHECK_GEOMETRY_EVEN(out);

  // alloc memory for temporary spinor field
  if (init) { init_Dirac(); }

  // Operation
  apply_BCs_on_spinor_field_flt(in);
  Dphi_flt_gpu_(otmp, in);
  apply_BCs_on_spinor_field_flt(otmp);
  Dphi_flt_gpu_(out, otmp);
  float rho = 4.f + (double)m0;
  rho*=-rho; /* this minus sign is taken into account below */
  spinor_field_mul_add_assign_f_flt_gpu(out, rho, in);
  spinor_field_minus_f_flt_gpu(out, out);
  apply_BCs_on_spinor_field_flt(out);
}

/**
  * @brief Even-odd preconditioned Wilson-Dirac operator with mass on the GPU.
  *         Single precision version
  *
  * @param in                 Odd input spinor field
  * @param out                Odd output spinor field that is the result of the application
  *                           of the even-odd preconditioned Wilson-Dirac operator on the 
  *                           input spinor field
  * @param m0                 Mass parameter
  */
void Dphi_oepre_flt_gpu(double m0, spinor_field_flt *out, spinor_field_flt *in)
{
  // Input argument validity checks
  _CHECK_GEOMETRY_ODD(in);
  _CHECK_GEOMETRY_ODD(out);

  // alloc memory for temporary spinor field 
  if (init) { init_Dirac();}

  apply_BCs_on_spinor_field_flt(in);
  Dphi_flt_gpu_(etmp, in);
  apply_BCs_on_spinor_field_flt(etmp);
  Dphi_flt_gpu_(out, etmp);
  float rho = 4.f + (float)m0;
  rho*=-rho; /* this minus sign is taken into account below */
  spinor_field_mul_add_assign_f_flt_gpu(out, rho, in);
  spinor_field_minus_f_flt_gpu(out, out);
  apply_BCs_on_spinor_field_flt(out);
}

/**
  * @brief Even-odd preconditioned Hermitian Wilson-Dirac operator with mass on the GPU.
  *        Single precision version.
  *
  * @param in                 Even input spinor field
  * @param out                Even output spinor field that is the result of the application
  *                           of the even-odd preconditioned Wilson-Dirac operator on the 
  *                           input spinor field
  * @param m0                 Mass parameter
  */
void g5Dphi_eopre_flt_gpu(double m0, spinor_field_flt *out, spinor_field_flt *in)
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
  apply_BCs_on_spinor_field_flt(in);
  Dphi_flt_gpu_(otmp, in);
  apply_BCs_on_spinor_field_flt(otmp);
  Dphi_flt_gpu_(out, otmp);
  float rho = 4.f + (float)m0;
  rho*=-rho; /* this minus sign is taken into account below */
  spinor_field_mul_add_assign_f_flt_gpu(out, rho, in);
  spinor_field_minus_f_flt_gpu(out, out);
  spinor_field_g5_assign_f_flt_gpu(out);
  apply_BCs_on_spinor_field_flt(out);

  //#if defined(BASIC_SF) || defined(ROTATED_SF)
  //  SF_spinor_bcs(out);
  //#endif 
}

/**
  * @brief Even-odd preconditioned squared Hermitian Wilson-Dirac operator with mass on the GPU.
  *         Single precision version.
  *
  * @param in                 Even input spinor field
  * @param out                Even output spinor field that is the result of the application
  *                           of the even-odd preconditioned Wilson-Dirac operator on the 
  *                           input spinor field
  * @param m0                 Mass parameter
  */
void g5Dphi_eopre_sq_flt_gpu(double m0, spinor_field_flt *out, spinor_field_flt *in) {
  // alloc memory for temporary spinor field
  if (init) { init_Dirac(); }

  g5Dphi_eopre_flt_gpu(m0, etmp, in);
  g5Dphi_eopre_flt_gpu(m0, out, etmp);
}

/**
 * @brief Implementation of the squared Hermitian Wilson-Dirac operator with mass on the GPU.
 *        Single precision version.
 *
 * @param in                    Input spinor field defined on the full lattice
 * @param out                   Output spinor field defined on the full lattice that is 
 *                              the result of the application of the Hermitian Wilson-Dirac 
 *                              operator on the input spinor field
 * @param m0                    Mass parameter
 */
void g5Dphi_sq_flt_gpu(double m0, spinor_field_flt *out, spinor_field_flt *in) {
  // alloc memory for temporary spinor field 
  if (init) { init_Dirac();  }

  g5Dphi_flt_gpu(m0, gtmp, in);
  g5Dphi_flt_gpu(m0, out, gtmp);
}

#ifdef WITH_GPU
  unsigned long int (*getMVM_flt) ()=getMVM_flt_gpu;
  void (*Dphi_flt_) (spinor_field_flt *out, spinor_field_flt *in)=Dphi_flt_gpu_;
  void (*Dphi_flt) (double m0, spinor_field_flt *out, spinor_field_flt *in)=Dphi_flt_gpu;
  void (*g5Dphi_flt) (double m0, spinor_field_flt *out, spinor_field_flt *in)=g5Dphi_flt_gpu;
  void (*g5Dphi_sq_flt) (double m0, spinor_field_flt *out, spinor_field_flt *in)=g5Dphi_sq_flt_gpu;
  void (*Dphi_eopre_flt) (double m0, spinor_field_flt *out, spinor_field_flt *in)=Dphi_eopre_flt_gpu;
  void (*Dphi_oepre_flt) (double m0, spinor_field_flt *out, spinor_field_flt *in)=Dphi_oepre_flt_gpu;
  void (*g5Dphi_eopre_flt) (double m0, spinor_field_flt *out, spinor_field_flt *in)=g5Dphi_eopre_flt_gpu;
  void (*g5Dphi_eopre_sq_flt) (double m0, spinor_field_flt *out, spinor_field_flt *in)=g5Dphi_eopre_sq_flt_gpu;
#endif

#endif
