/***************************************************************************\
* Copyright (c) 2008, 2022, Claudio Pica, Sofie Martins                     *
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

static kernel_field_input* get_even_input(spinor_field *out, spinor_field *in) {
  kernel_field_input* input;
  input = (kernel_field_input*)malloc(sizeof(kernel_field_input));

  input->field_in = in->gpu_ptr;
  input->master_shift_in = in->type->master_shift;
  input->start_in = geometryBoxes->base_index;
  input->stride_in = boxEvenVolume(geometryBoxes);

  input->field_out = out->gpu_ptr;
  input->master_shift_out = out->type->master_shift;
  input->start_out = geometryBoxes->base_index_odd;
  input->stride_out = boxOddVolume(geometryBoxes);

  kernel_field_input* d_input;
  cudaMalloc((void**)&d_input, sizeof(kernel_field_input));
  cudaMemcpy(d_input, input, sizeof(kernel_field_input), cudaMemcpyHostToDevice);
  return d_input;
}

static kernel_field_input* get_even_buffer_input(spinor_field *out, spinor_field *in, int buffer_idx) {
  kernel_field_input* input;
  input = (kernel_field_input*)malloc(sizeof(kernel_field_input));

  input->field_in = in->gpu_ptr;
  input->start_in = in->type->rbuf_start[buffer_idx];
  input->stride_in = in->type->rbuf_len[buffer_idx];
  input->master_shift_in = in->type->master_shift;

  input->field_out = out->gpu_ptr;
  input->start_out = geometryBoxes->base_index_odd;
  input->stride_out = boxOddVolume(geometryBoxes);
  input->master_shift_out = out->type->master_shift;

  kernel_field_input* d_input;
  cudaMalloc((void**)&d_input, sizeof(kernel_field_input));
  cudaMemcpy(d_input, input, sizeof(kernel_field_input), cudaMemcpyHostToDevice);
  return d_input;
}

static kernel_field_input* get_odd_input(spinor_field *out, spinor_field *in) {
  kernel_field_input *input;
  input = (kernel_field_input*)malloc(sizeof(kernel_field_input));

  input->field_in = in->gpu_ptr;
  input->start_in = geometryBoxes->base_index_odd;
  input->stride_in = boxOddVolume(geometryBoxes);
  input->master_shift_in = in->type->master_shift;

  input->field_out = out->gpu_ptr;
  input->start_out = geometryBoxes->base_index;
  input->stride_out = boxEvenVolume(geometryBoxes);
  input->master_shift_out = out->type->master_shift;

  kernel_field_input* d_input;
  cudaMalloc((void**)&d_input, sizeof(kernel_field_input));
  cudaMemcpy(d_input, input, sizeof(kernel_field_input), cudaMemcpyHostToDevice);
  return d_input;
}

static kernel_field_input* get_odd_buffer_input(spinor_field *out, spinor_field *in, int buffer_idx) {
  kernel_field_input* input;
  input = (kernel_field_input*)malloc(sizeof(kernel_field_input));
  
  int i = buffer_idx;
  if (out->type==&glattice) i = buffer_idx+1;
  input->field_in = in->gpu_ptr;
  input->start_in = in->type->rbuf_start[i];
  input->stride_in = in->type->rbuf_len[i];
  input->master_shift_in = in->type->master_shift;

  input->field_out = out->gpu_ptr;
  input->start_out = geometryBoxes->base_index;
  input->stride_out = boxEvenVolume(geometryBoxes);
  input->master_shift_out = out->type->master_shift;

  kernel_field_input* d_input;
  cudaMalloc((void**)&d_input, sizeof(kernel_field_input));
  cudaMemcpy(d_input, input, sizeof(kernel_field_input), cudaMemcpyHostToDevice);
  return d_input;
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
  enum gd_type gd_t;
  if (in->type==&glattice) gd_t = GLOBAL;
  else if (in->type==&glat_odd) gd_t = ODD;
  else if (in->type==&glat_even) gd_t = EVEN;

  kernel_field_input* input_even = get_even_input(out, in);
  kernel_field_input* input_odd = get_odd_input(out, in);
  int grid = (boxEvenVolume(geometryBoxes)-1)/BLOCK_SIZE_DIRAC + 1;
  Dphi_gpu_inner_kernel<<<grid, BLOCK_SIZE>>>(input_even, input_odd, u_gauge_f->gpu_ptr, iup_gpu, idn_gpu, imask_gpu, gd_t);
  CudaCheckError();
  //cudaDeviceSynchronize();
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
  enum gd_type gd_t;
  if (in->type==&glattice) gd_t = GLOBAL;
  else if (in->type==&glat_odd) gd_t = ODD;
  else if (in->type==&glat_even) gd_t = EVEN;

  box_t *buffer_box = geometryBoxes->next;
  int i = 0;
  int buffer_index;
  int nbuffers = in->type->nbuffers_spinor;
  int increment = 1;
  if (gd_t == GLOBAL)

  if(gd_t == GLOBAL) nbuffers /= 2;
  while(buffer_box && i < nbuffers) {
    if (gd_t == GLOBAL) {
      buffer_index = 2*i;
    } else {
      buffer_index = i;
    }
    kernel_field_input* input_even = get_even_buffer_input(out, in, buffer_index);
    kernel_field_input* input_odd = get_odd_buffer_input(out, in, buffer_index);
    int grid = (boxEvenVolume(geometryBoxes)-1)/BLOCK_SIZE_DIRAC + 1;
    Dphi_gpu_boundary_kernel<<<grid, BLOCK_SIZE>>>(input_even, input_odd, u_gauge_f->gpu_ptr, iup_gpu, idn_gpu, imask_gpu, gd_t);
    buffer_box=buffer_box->next; i++;
  }
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
