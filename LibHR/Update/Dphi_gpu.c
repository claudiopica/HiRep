/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *
* All rights reserved.                                                      *
\***************************************************************************/

/*******************************************************************************
*
* File Dphi_gpu.c
*
* Action of the Wilson-Dirac operator D and hermitian g5D on a given
* double-precision spinor field
*
*******************************************************************************/


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
#include "gpu_geometry.h"
#include "communications.h"
#include "memory.h"
#include "gpu.h"
#include "hr_complex.h"

__global__ void Dphi_gpu_kernel(suNf_spinor*,
                            const suNf_spinor*,
                            const suNf*,
                            const int*,
                            const int*,
                            const int,
                            const int,
                            const int,
                            const int,
                            const int,
                            const int*,
                            const int*,
                            const int*, 
                            const int*,
                            const int*,
                            const int*,
                            const int, 
                            int);

#ifdef ROTATED_SF
#include "update.h"
extern rhmc_par _update_par; /* Update/update_rhmc.c */
#endif /* ROTATED_SF */

//__device__ __constant__ hr_complex eitheta_gpu[4];

/*
 * the following variable is used to keep trace of
 * matrix-vector multiplication.
 * we count how many time the function Dphi_ is called
 */
static unsigned long int MVMcounter=0;

unsigned long int getMVM_gpu() {
	unsigned long int res=MVMcounter>>1; /* divide by two */
	//MVMcounter=0; /* reset counter */
	return res;
}

/* r=t*u*s */
#ifdef BC_T_THETA

#define _suNf_theta_T_multiply(r, u, s)\
    _suNf_multiply(vtmp, (u), (s));\
    _vector_mulc_f((r), eitheta_gpu[0], vtmp)

#define _suNf_theta_T_inverse_multiply(r, u, s)\
    _suNf_inverse_multiply(vtmp, (u), (s));\
    _vector_mulc_star_f((r), eitheta_gpu[0], vtmp)

#else

#define _suNf_theta_T_multiply(r, u, s) _suNf_multiply((r), (u), (s))
#define _suNf_theta_T_inverse_multiply(r, u, s) _suNf_inverse_multiply((r), (u), (s))

#endif

/* r=t*u*s */
#ifdef BC_X_THETA

#define _suNf_theta_X_multiply(r, u, s)\
_suNf_multiply(vtmp, (u), (s));\
_vector_mulc_f((r), eitheta_gpu[1], vtmp)

#define _suNf_theta_X_inverse_multiply(r, u, s)\
_suNf_inverse_multiply(vtmp, (u), (s));\
_vector_mulc_star_f((r), eitheta_gpu[1], vtmp)

#else

#define _suNf_theta_X_multiply(r, u, s) _suNf_multiply((r), (u), (s))
#define _suNf_theta_X_inverse_multiply(r, u, s) _suNf_inverse_multiply((r), (u), (s))

#endif

/* r=t*u*s */
#ifdef BC_Y_THETA

#define _suNf_theta_Y_multiply(r, u, s)\
_suNf_multiply(vtmp, (u), (s));\
_vector_mulc_f((r), eitheta_gpu[2], vtmp)

#define _suNf_theta_Y_inverse_multiply(r, u, s)\
_suNf_inverse_multiply(vtmp, (u), (s));\
_vector_mulc_star_f((r), eitheta_gpu[2], vtmp)

#else

#define _suNf_theta_Y_multiply(r, u, s) _suNf_multiply((r), (u), (s))
#define _suNf_theta_Y_inverse_multiply(r, u, s) _suNf_inverse_multiply((r), (u), (s))

#endif

/* r=t*u*s */
#ifdef BC_Z_THETA

#define _suNf_theta_Z_multiply(r, u, s)\
_suNf_multiply(vtmp, (u), (s));\
_vector_mulc_f((r), eitheta_gpu[3], vtmp)

#define _suNf_theta_Z_inverse_multiply(r, u, s)\
_suNf_inverse_multiply(vtmp, (u), (s));\
_vector_mulc_star_f((r), eitheta_gpu[3], vtmp)

#else

#define _suNf_theta_Z_multiply(r, u, s) _suNf_multiply((r), (u), (s))
#define _suNf_theta_Z_inverse_multiply(r, u, s) _suNf_inverse_multiply((r), (u), (s))

#endif

typedef struct _suNf_hspinor
{
  suNf_vector c[2];
} suNf_hspinor;

#define THREADSITE 1

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

void Dphi_gpu_(spinor_field *out, spinor_field *in)
{
  unsigned int N, grid;
  const int vol4h=T*X*Y*Z/2;
  const int test_block=vol4h/2;

  init_bc_gpu();

  error((in==NULL)||(out==NULL), 1, "Dphi_gpu_ [Dphi_gpu.c]",
         "Attempt to access unallocated memory space");

  error(in==out, 1, "Dphi_gpu_ [Dphi_gpu.c]",
         "Input and output fields must be different");

  #ifndef CHECK_SPINOR_MATCHING
    error(out->type==&glat_even && in->type!=&glat_odd, 1, "Dphi_gpu_ [Dphi_gpu.c]", "Spinors don't match! (1)");
    error(out->type==&glat_odd && in->type!=&glat_even, 1, "Dphi_gpu_ [Dphi_gpu.c]", "Spinors don't match! (2)");
    error(out->type==&glattice && in->type!=&glattice, 1, "Dphi_gpu_ [Dphi_gpu.c]", "Spinors don't match! (3)");
  #endif


  #ifdef WITH_MPI
    // Sync + communications, TODO: Not all of these might be necessary
    sync_gpu_spinor_field_f(in);
    sync_gpu_spinor_field_f(out);
    sync_gpu_gfield_f(u_gauge_f);

    start_sendrecv_gpu_spinor_field_f(in);
    /*start_sendrecv_gpu_spinor_field_f(out);
    start_sendrecv_gpu_gfield_f(u_gauge_f);*/

    complete_sendrecv_gpu_spinor_field_f(in);
    /*complete_sendrecv_gpu_spinor_field_f(out);
    complete_sendrecv_gpu_gfield_f(u_gauge_f);*/

    //fill_buffers_spinor_field_f(in);
    //fill_buffers_with_zeroes_spinor_field_f(in);
  #endif

  int* master_start;
  cudaMalloc((int**)&master_start, in->type->local_master_pieces*sizeof(int));
  cudaMemcpy(master_start, in->type->master_start, in->type->local_master_pieces*sizeof(int), cudaMemcpyHostToDevice);

  int* master_end;
  cudaMalloc((int**)&master_end, in->type->local_master_pieces*sizeof(int));
  cudaMemcpy(master_end, in->type->master_end, in->type->local_master_pieces*sizeof(int), cudaMemcpyHostToDevice);

  int* rbuf_start;
  cudaMalloc((int**)&rbuf_start, in->type->nbuffers_spinor*sizeof(int));// What about gauge field buffers?
  cudaMemcpy(rbuf_start, in->type->rbuf_start, in->type->nbuffers_spinor*sizeof(int), cudaMemcpyHostToDevice);

  int* rbuf_len;
  cudaMalloc((int**)&rbuf_len, in->type->nbuffers_spinor*sizeof(int));
  cudaMemcpy(rbuf_len, in->type->rbuf_len, in->type->nbuffers_spinor*sizeof(int), cudaMemcpyHostToDevice);

  int* sbuf_start;
  cudaMalloc((int**)&sbuf_start, in->type->nbuffers_spinor*sizeof(int));
  cudaMemcpy(sbuf_start, in->type->sbuf_start, in->type->nbuffers_spinor*sizeof(int), cudaMemcpyHostToDevice);

  int* sbuf_len;
  cudaMalloc((int**)&sbuf_len, in->type->nbuffers_spinor*sizeof(int));
  cudaMemcpy(sbuf_len, in->type->sbuf_start, in->type->nbuffers_spinor*sizeof(int), cudaMemcpyHostToDevice);

  // TODO: split up into boundary and bulk calculation -> Hide communications behind the bulk calculation
  _INNER_PIECE_FOR(out->type, ixp) 
  {
      int write_stride = out->type->master_end[ixp] - out->type->master_start[ixp] + 1;
      int write_start = out->type->master_start[ixp];
      grid = (write_stride-1)/BLOCK_SIZE + 1;
      printf("write_start: %d\n", write_start);
      Dphi_gpu_kernel<<<grid, BLOCK_SIZE>>>(_GPU_FIELD_BLK(out, ixp),
                                            in->gpu_ptr,
                                            u_gauge_f->gpu_ptr,
                                            iup_gpu, 
                                            idn_gpu,
                                            write_start,
                                            write_stride,
                                            in->type->inner_master_pieces, 
                                            in->type->local_master_pieces-in->type->inner_master_pieces,
                                            in->type->nbuffers_spinor,
                                            master_start, 
                                            master_end,
                                            sbuf_start,
                                            sbuf_len,
                                            rbuf_start, 
                                            rbuf_len,
                                            ixp, PID);
      CudaCheckError();
  }

  // This is very inefficient, find a way to merge these kernels
  for (int i = 0; i < in->type->nbuffers_spinor; ++i) 
  {
    int write_stride = out->type->sbuf_len[i];
    grid = (write_stride-1)/BLOCK_SIZE + 1;
    Dphi_gpu_kernel<<<grid, BLOCK_SIZE>>>(out->gpu_ptr + out->type->sbuf_start[i], 
                                          in->gpu_ptr, 
                                          u_gauge_f->gpu_ptr, 
                                          iup_gpu, 
                                          idn_gpu, 
                                          out->type->sbuf_start[i], 
                                          write_stride, 
                                          in->type->inner_master_pieces, 
                                          in->type->local_master_pieces-in->type->inner_master_pieces, 
                                          in->type->nbuffers_spinor, 
                                          master_start, 
                                          master_end, 
                                          sbuf_start, 
                                          sbuf_len,
                                          rbuf_start, 
                                          rbuf_len, 
                                          0, PID);
  }
}


/*
 * this function takes 2 spinors defined on the whole lattice
 */
void Dphi_gpu(double m0, spinor_field *out, spinor_field *in)
{
  double rho;

  error((in==NULL)||(out==NULL), 1, "Dphi_gpu [Dphi_gpu.c]",
        "Attempt to access unallocated memory space");

  error(in==out, 1, "Dphi_gpu [Dphi_gpu.c]",
        "Input and output fields must be different");

#ifdef CHECK_SPINOR_MATCHING
  error(out->type!=&glattice || in->type!=&glattice, 1, "Dphi_gpu [Dphi_gpu.c]", "Spinors are not defined on all the lattice!");
#endif /* CHECK_SPINOR_MATCHING */

  Dphi_gpu_(out, in);

  rho = 4. + m0;
  spinor_field_mul_add_assign_f_gpu(out, rho, in);
}


void g5Dphi_gpu(double m0, spinor_field *out, spinor_field *in)
{
  double rho;

  error((in==NULL)||(out==NULL), 1, "g5Dphi_gpu [Dphi_gpu.c]",
	"Attempt to access unallocated memory space");

  error(in==out, 1, "g5Dphi_gpu [Dphi_gpu.c]",
	"Input and output fields must be different");

#ifdef CHECK_SPINOR_MATCHING
   error(out->type!=&glattice || in->type!=&glattice, 1, "g5Dphi_gpu [Dphi_gpu.c]", "Spinors are not defined on all the lattice!");
#endif /* CHECK_SPINOR_MATCHING */

  Dphi_gpu_(out, in);
  rho=4.+m0;
  spinor_field_mul_add_assign_f_gpu(out, rho, in);
  spinor_field_g5_assign_f_gpu(out);
}


static int init=1;
static spinor_field *gtmp=NULL;
static spinor_field *etmp=NULL;
static spinor_field *otmp=NULL;

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


/* Even/Odd preconditioned dirac operator
 * this function takes 2 spinors defined on the even lattice
 * Dphi in = (4+m0)^2*in - D_EO D_OE in
 *
 */
void Dphi_eopre_gpu(double m0, spinor_field *out, spinor_field *in)
{
  double rho;

  error((in==NULL)||(out==NULL), 1, "Dphi_eopre_gpu [Dphi_gpu.c]",
	"Attempt to access unallocated memory space");

  error(in==out, 1, "Dphi_eopre_gpu [Dphi_gpu.c]",
	"Input and output fields must be different");

#ifdef CHECK_SPINOR_MATCHING
  error(out->type!=&glat_even || in->type!=&glat_even, 1, "Dphi_eopre_gpu " __FILE__, "Spinors are not defined on even lattice!");
#endif /* CHECK_SPINOR_MATCHING */

  /* alloc memory for temporary spinor field */
  if (init) { init_Dirac(); }

  Dphi_gpu_(otmp, in);
  Dphi_gpu_(out, otmp);

  rho=4.0+m0;
  rho*=-rho; /* this minus sign is taken into account below */

  spinor_field_mul_add_assign_f_gpu(out, rho, in);
  spinor_field_minus_f_gpu(out, out);
}


/* Even/Odd preconditioned dirac operator
 * this function takes 2 spinors defined on the odd lattice
 * Dphi in = (4+m0)^2*in - D_OE D_EO in
 *
 */
void Dphi_oepre_gpu(double m0, spinor_field *out, spinor_field *in)
{
  double rho;

  error((in==NULL)||(out==NULL), 1, "Dphi_oepre_gpu [Dphi_gpu.c]",
	"Attempt to access unallocated memory space");

  error(in==out, 1, "Dphi_oepre_gpu [Dphi_gpu.c]",
	"Input and output fields must be different");

#ifdef CHECK_SPINOR_MATCHING
  error(out->type!=&glat_odd || in->type!=&glat_odd, 1, "Dphi_oepre_gpu " __FILE__, "Spinors are not defined on odd lattice!");
#endif /* CHECK_SPINOR_MATCHING */


  /* alloc memory for temporary spinor field */
  if (init) { init_Dirac();}

  Dphi_gpu_(etmp, in);
  Dphi_gpu_(out, etmp);

  rho=4.0+m0;
  rho*=-rho; /* this minus sign is taken into account below */

  spinor_field_mul_add_assign_f_gpu(out, rho, in);
  spinor_field_minus_f_gpu(out, out);
}


void g5Dphi_eopre_gpu(double m0, spinor_field *out, spinor_field *in)
{
  double rho;

  error((in==NULL)||(out==NULL), 1, "g5Dphi_eopre_gpu [Dphi_gp.c]",
	"Attempt to access unallocated memory space");

  error(in==out, 1, "Dphi_eopre_gpu [Dphi_gpu.c]",
	"Input and output fields must be different");

#ifdef CHECK_SPINOR_MATCHING
  error(out->type!=&glat_even || in->type!=&glat_even, 1, "g5Dphi_eopre_gpu " __FILE__, "Spinors are not defined on even lattice!");
#endif /* CHECK_SPINOR_MATCHING */

#if defined(BASIC_SF) || defined(ROTATED_SF)
  SF_spinor_bcs(in);
#endif /* defined(BASIC_SF) || defined(ROTATED_SF) */

  /* alloc memory for temporary spinor field */
  if (init) { init_Dirac();}

  Dphi_gpu_(otmp, in);
  Dphi_gpu_(out, otmp);

  rho=4.0+m0;
  rho*=-rho; /* this minus sign is taken into account below */

  spinor_field_mul_add_assign_f_gpu(out, rho, in);
  spinor_field_minus_f_gpu(out, out);
  spinor_field_g5_assign_f_gpu(out);

#if defined(BASIC_SF) || defined(ROTATED_SF)
  SF_spinor_bcs(out);
#endif /* defined(BASIC_SF) || defined(ROTATED_SF) */
}


/* g5Dphi_eopre ^2 */
void g5Dphi_eopre_sq_gpu(double m0, spinor_field *out, spinor_field *in) 
{
  /* alloc memory for temporary spinor field */
  if (init) { init_Dirac(); }

  g5Dphi_eopre_gpu(m0, etmp, in);
  g5Dphi_eopre_gpu(m0, out, etmp);
}


/* g5Dhi ^2 */
void g5Dphi_sq_gpu(double m0, spinor_field *out, spinor_field *in) 
{
  /* alloc memory for temporary spinor field */
  if (init) { init_Dirac();  }

  g5Dphi_gpu(m0, gtmp, in);
  g5Dphi_gpu(m0, out, gtmp);
}

/* ================================== KERNEL ================================== */

/* Takes an even input spinor and returns an odd spinor */
__global__ void Dphi_gpu_kernel(suNf_spinor* __restrict__ out,
                            const suNf_spinor* __restrict__ in,
                            const suNf* __restrict__ gauge,
                            const int* __restrict__ iup_d,
                            const int* __restrict__ idn_d,
                            const int write_start,//TODO: Figure these out from master start, end+ixp
                            const int write_stride,
                            const int ninnerpieces, 
                            const int nboundarypieces,
                            const int nbuffers,
                            const int* master_start, 
                            const int* master_end,
                            const int* sbuf_start, 
                            const int* sbuf_len,
                            const int* rbuf_start, 
                            const int* rbuf_len,
                            const int ixp, int PID)
{

  suNf_spinor r;
  suNf_hspinor sn;
  suNf u;
  #ifdef FERMION_THETA
    suNf_vector vtmp;
  #endif

  int ix, iy, local_iy, iyp;
  int read_start, read_stride, tmp;
  int local_ix = blockIdx.x * BLOCK_SIZE + threadIdx.x;

  if (local_ix < write_stride) {
    ix = local_ix + write_start;
    /******************************* direction +0 *********************************/
    iy=iup_d[4*ix]; 
    read_start = -1;
    read_stride = -1;

    int i;
    for (i = 0; i < ninnerpieces+1; ++i) {
      if (i == ninnerpieces) break;
      if (iy >= master_start[i] && iy <= master_end[i])
      {
        read_start = master_start[i];
        read_stride = master_end[i] - master_start[i] + 1;
      }
    }
  
    if (i == ninnerpieces) /* index on boundary */
    {
      for (i = 0; i < nboundarypieces+1; ++i) 
      {
        if (i == ninnerpieces) break;
        if (iy >= sbuf_start[i] && iy < (sbuf_start[i] + sbuf_len[i]))
        {
          read_start = sbuf_start[i];
          read_stride = sbuf_len[i];
        }
      }
    } 

    if (i == nboundarypieces) /* index in receive buffers */
    { 
      for (i = 0; i < nbuffers; ++i) 
      { 
        if (iy >= rbuf_start[i] && iy < (rbuf_start[i] + rbuf_len[i])) 
        { 
          read_start = rbuf_start[i]; 
          read_stride = rbuf_len[i]; 
        } 
      } 
    } 

    printf("iy: %d, read_start: %d, read_stride: %d\n", iy, read_start, read_stride);

    const suNf_spinor *in_blk = in + read_start;
    local_iy = iy - read_start;

    read_gpu_suNf_vector(read_stride, sn.c[0], in_blk, local_iy, 0);
    read_gpu_suNf_vector(read_stride, sn.c[1], in_blk, local_iy, 2);
    read_gpu_suNf(write_stride, u, gauge + 4*write_start, local_ix, 0);

    if (ix == 5) {
      printf("Read start: %d (should be 288)\n", read_start);
      printf("sbuf PID: %d, GPU buf: %0.2e + i%0.2e\n", 
            PID, creal(sn.c[0].c[0]), cimag(sn.c[0].c[0]));
      double *send_buffer = (double*)(in + 288);
      printf("PID: %d, sbuf: %0.2e\n", PID, send_buffer[0]);
    }

    _vector_add_assign_f(sn.c[0], sn.c[1]);
    _suNf_theta_T_multiply(r.c[0], u, sn.c[0]);

    r.c[2]=r.c[0];

    read_gpu_suNf_vector(read_stride, sn.c[0], in_blk, local_iy, 1);
    read_gpu_suNf_vector(read_stride, sn.c[1], in_blk, local_iy, 3);

    _vector_add_assign_f(sn.c[0], sn.c[1]);
    _suNf_theta_T_multiply(r.c[1], u, sn.c[0]);

    r.c[3]=r.c[1];

    //if (ix == 69) printf("GPU res(0): %0.2e + i%0.2e\n", 
    //  creal(r.c[0].c[0]), cimag(r.c[0].c[0]));

    __syncthreads();
    /******************************* direction -0 *********************************/
    iy=idn_d[4*ix];

    for (i = 0; i < ninnerpieces+1; ++i) {
      if (i == ninnerpieces) break;
      if (iy >= master_start[i] && iy <= master_end[i])
      {
        read_start = master_start[i];
        read_stride = master_end[i] - master_start[i] + 1;
      }
    }
  
    if (i == ninnerpieces) /* index in buffer */ 
    { 
      for (i = 0; i < nbuffers; ++i) 
      { 
        if (iy >= rbuf_start[i] && iy < (rbuf_start[i] + rbuf_len[i])) 
        { 
          read_start = rbuf_start[i]; 
          read_stride = rbuf_len[i]; 
        } 
      } 
    } 

    in_blk = in + read_start;
    local_iy = iy - read_start;

    read_gpu_suNf_vector(read_stride, sn.c[0], in_blk, local_iy, 0);
    read_gpu_suNf_vector(read_stride, sn.c[1], in_blk, local_iy, 2);
    read_gpu_suNf(read_stride, u, gauge + 4*read_start, iy-read_start, 0);

    if (iy == 5) printf("rbuf PID: %d, GPU buf: %0.2e + i%0.2e\n", 
                PID, creal(sn.c[0].c[0]), cimag(sn.c[0].c[0]));

    //if (ix == 69) printf("GPU gauge: %0.2e + i%0.2e\n", 
    //            creal(u.c[0]), cimag(u.c[0]));

    _vector_sub_assign_f(sn.c[0], sn.c[1]);
    _suNf_theta_T_inverse_multiply(sn.c[1], u, sn.c[0]);

    _vector_add_assign_f(r.c[0], sn.c[1]);
    _vector_sub_assign_f(r.c[2], sn.c[1]);

    read_gpu_suNf_vector(read_stride, sn.c[0], in_blk, local_iy, 1);
    read_gpu_suNf_vector(read_stride, sn.c[1], in_blk, local_iy, 3);
    _vector_sub_assign_f(sn.c[0], sn.c[1]);

    _suNf_theta_T_inverse_multiply(sn.c[1], u, sn.c[0]);

    _vector_add_assign_f(r.c[1], sn.c[1]);
    _vector_sub_assign_f(r.c[3], sn.c[1]);

    //if (ix == 69) printf("GPU res(0): %0.2e + i%0.2e\n", 
    //  creal(r.c[0].c[0]), cimag(r.c[0].c[0]));

    __syncthreads();
    /******************************* direction +1 *********************************/
    iy=iup_d[4*ix+1];

    for (i = 0; i < ninnerpieces+1; ++i) {
      if (i == ninnerpieces) break;
      if (iy >= master_start[i] && iy <= master_end[i])
      {
        read_start = master_start[i];
        read_stride = master_end[i] - master_start[i] + 1;
      }
    }
  
    if (i == ninnerpieces) /* index in buffer */ 
    { 
      for (i = 0; i < nbuffers; ++i) 
      { 
        if (iy >= rbuf_start[i] && iy < (rbuf_start[i] + rbuf_len[i])) 
        { 
          read_start = rbuf_start[i]; 
          read_stride = rbuf_len[i]; 
        } 
      } 
    } 

    in_blk = in + read_start;
    local_iy = iy - read_start;

    read_gpu_suNf_vector(read_stride, sn.c[0], in_blk, local_iy, 0);
    read_gpu_suNf_vector(read_stride, sn.c[1], in_blk, local_iy, 3);
    read_gpu_suNf(write_stride, u, gauge + 4*write_start, local_ix, 1);

    _vector_i_add_assign_f(sn.c[0], sn.c[1]);
    _suNf_theta_X_multiply(sn.c[1], u, sn.c[0]);

    _vector_add_assign_f(r.c[0], sn.c[1]);
    _vector_i_sub_assign_f(r.c[3], sn.c[1]);

    read_gpu_suNf_vector(read_stride, sn.c[0], in_blk, local_iy, 1);
    read_gpu_suNf_vector(read_stride, sn.c[1], in_blk, local_iy, 2);
    _vector_i_add_assign_f(sn.c[0], sn.c[1]);

    _suNf_theta_X_multiply(sn.c[1], u, sn.c[0]);

    _vector_add_assign_f(r.c[1], sn.c[1]);
    _vector_i_sub_assign_f(r.c[2], sn.c[1]);

    __syncthreads();
    /******************************* direction -1 *********************************/
    iy=idn_d[4*ix+1];

    for (i = 0; i < ninnerpieces+1; ++i) {
      if (i == ninnerpieces) break;
      if (iy >= master_start[i] && iy <= master_end[i])
      {
        read_start = master_start[i];
        read_stride = master_end[i] - master_start[i] + 1;
      }
    }
  
    if (i == ninnerpieces) /* index in buffer */ 
    { 
      for (i = 0; i < nbuffers; ++i) 
      { 
        if (iy >= rbuf_start[i] && iy < (rbuf_start[i] + rbuf_len[i])) 
        { 
          read_start = rbuf_start[i]; 
          read_stride = rbuf_len[i]; 
        } 
      } 
    } 

    in_blk = in + read_start;
    local_iy = iy - read_start;

    read_gpu_suNf_vector(read_stride, sn.c[0], in_blk, local_iy, 0);
    read_gpu_suNf_vector(read_stride, sn.c[1], in_blk, local_iy, 3);
    read_gpu_suNf(read_stride, u, gauge + 4*read_start, iy-read_start, 1);

    _vector_i_sub_assign_f(sn.c[0], sn.c[1]);
    _suNf_theta_X_inverse_multiply(sn.c[1], u, sn.c[0]);

    _vector_add_assign_f(r.c[0], sn.c[1]);
    _vector_i_add_assign_f(r.c[3], sn.c[1]);

    read_gpu_suNf_vector(read_stride, sn.c[0], in_blk, local_iy, 1);
    read_gpu_suNf_vector(read_stride, sn.c[1], in_blk, local_iy, 2);
    _vector_i_sub_assign_f(sn.c[0], sn.c[1]);

    _suNf_theta_X_inverse_multiply(sn.c[1], u, sn.c[0]);

    _vector_add_assign_f(r.c[1], sn.c[1]);
    _vector_i_add_assign_f(r.c[2], sn.c[1]);

    __syncthreads();
    /******************************* direction +2 *********************************/
    iy=iup_d[4*ix+2];

    for (i = 0; i < ninnerpieces+1; ++i) {
      if (i == ninnerpieces) break;
      if (iy >= master_start[i] && iy <= master_end[i])
      {
        read_start = master_start[i];
        read_stride = master_end[i] - master_start[i] + 1;
      }
    }
  
    if (i == ninnerpieces) /* index in buffer */ 
    { 
      for (i = 0; i < nbuffers; ++i) 
      { 
        if (iy >= rbuf_start[i] && iy < (rbuf_start[i] + rbuf_len[i])) 
        { 
          read_start = rbuf_start[i]; 
          read_stride = rbuf_len[i]; 
        } 
      } 
    } 

    in_blk = in + read_start;
    local_iy = iy - read_start;

    read_gpu_suNf_vector(read_stride, sn.c[0], in_blk, local_iy, 0);
    read_gpu_suNf_vector(read_stride, sn.c[1], in_blk, local_iy, 3);
    //if (ix == 69) printf("GPU res spinor comp: %0.2e + i%0.2e\n", 
    //                      creal(sn.c[0].c[0]), cimag(sn.c[0].c[0]));

    read_gpu_suNf(write_stride, u, gauge + 4*write_start, local_ix, 2);

    //if (ix == 69) printf("GPU gauge: %0.2e + i%0.2e\n", 
    //                    creal(u.c[0]), cimag(u.c[0]));
    _vector_add_assign_f(sn.c[0], sn.c[1]);
    _suNf_theta_Y_multiply(sn.c[1], u, sn.c[0]);

    _vector_add_assign_f(r.c[0], sn.c[1]);
    _vector_add_assign_f(r.c[3], sn.c[1]);

    read_gpu_suNf_vector(read_stride, sn.c[0], in_blk, local_iy, 1);
    read_gpu_suNf_vector(read_stride, sn.c[1], in_blk, local_iy, 2);
    _vector_sub_assign_f(sn.c[0], sn.c[1]);

    _suNf_theta_Y_multiply(sn.c[1], u, sn.c[0]);

    _vector_add_assign_f(r.c[1], sn.c[1]);
    _vector_sub_assign_f(r.c[2], sn.c[1]);

    __syncthreads();
    /******************************* direction -2 *********************************/
    iy=idn_d[4*ix+2];

    for (i = 0; i < ninnerpieces+1; ++i) {
      if (i == ninnerpieces) break;
      if (iy >= master_start[i] && iy <= master_end[i])
      {
        read_start = master_start[i];
        read_stride = master_end[i] - master_start[i] + 1;
      }
    }
  
    if (i == ninnerpieces) /* index in buffer */ 
    { 
      for (i = 0; i < nbuffers; ++i) 
      { 
        if (iy >= rbuf_start[i] && iy < (rbuf_start[i] + rbuf_len[i])) 
        { 
          read_start = rbuf_start[i]; 
          read_stride = rbuf_len[i]; 
        } 
      } 
    } 

    in_blk = in + read_start;
    local_iy = iy - read_start;

    read_gpu_suNf_vector(read_stride, sn.c[0], in_blk, local_iy, 0);
    read_gpu_suNf_vector(read_stride, sn.c[1], in_blk, local_iy, 3);
    _vector_sub_assign_f(sn.c[0], sn.c[1]);

    read_gpu_suNf(read_stride, u, gauge + 4*read_start, iy-read_start, 2);
    _suNf_theta_Y_inverse_multiply(sn.c[1], u, sn.c[0]);

    _vector_add_assign_f(r.c[0], sn.c[1]);
    _vector_sub_assign_f(r.c[3], sn.c[1]);

    read_gpu_suNf_vector(read_stride, sn.c[0], in_blk, local_iy, 1);
    read_gpu_suNf_vector(read_stride, sn.c[1], in_blk, local_iy, 2);
    _vector_add_assign_f(sn.c[0], sn.c[1]);

    _suNf_theta_Y_inverse_multiply(sn.c[1], u, sn.c[0]);

    _vector_add_assign_f(r.c[1], sn.c[1]);
    _vector_add_assign_f(r.c[2], sn.c[1]);

    __syncthreads();
    /******************************* direction +3 *********************************/
    iy=iup_d[4*ix+3];

    for (i = 0; i < ninnerpieces+1; ++i) {
      if (i == ninnerpieces) break;
      if (iy >= master_start[i] && iy <= master_end[i])
      {
        read_start = master_start[i];
        read_stride = master_end[i] - master_start[i] + 1;
      }
    }
  
    if (i == ninnerpieces) /* index in buffer */ 
    { 
      for (i = 0; i < nbuffers; ++i) 
      { 
        if (iy >= rbuf_start[i] && iy < (rbuf_start[i] + rbuf_len[i])) 
        { 
          read_start = rbuf_start[i]; 
          read_stride = rbuf_len[i]; 
        } 
      } 
    } 

    in_blk = in + read_start;
    local_iy = iy - read_start;

    read_gpu_suNf_vector(read_stride, sn.c[0], in_blk, local_iy, 0);
    read_gpu_suNf_vector(read_stride, sn.c[1], in_blk, local_iy, 2);
    _vector_i_add_assign_f(sn.c[0], sn.c[1]);

    read_gpu_suNf(write_stride, u, gauge + 4*write_stride, local_ix, 3);
    _suNf_theta_Z_multiply(sn.c[1], u, sn.c[0]);

    _vector_add_assign_f(r.c[0], sn.c[1]);
    _vector_i_sub_assign_f(r.c[2], sn.c[1]);

    read_gpu_suNf_vector(read_stride, sn.c[0], in_blk, local_iy, 1);
    read_gpu_suNf_vector(read_stride, sn.c[1], in_blk, local_iy, 3);
    _vector_i_sub_assign_f(sn.c[0], sn.c[1]);

    _suNf_theta_Z_multiply(sn.c[1], u, sn.c[0]);

    _vector_add_assign_f(r.c[1], sn.c[1]);
    _vector_i_add_assign_f(r.c[3], sn.c[1]);

    __syncthreads();
    /******************************* direction -3 *********************************/
    iy=idn_d[4*ix+3];

    for (i = 0; i < ninnerpieces+1; ++i) {
      if (i == ninnerpieces) break;
      if (iy >= master_start[i] && iy <= master_end[i])
      {
        read_start = master_start[i];
        read_stride = master_end[i] - master_start[i] + 1;
      }
    }
  
    if (i == ninnerpieces) /* index in buffer */ 
    { 
      for (i = 0; i < nbuffers; ++i) 
      { 
        if (iy >= rbuf_start[i] && iy < (rbuf_start[i] + rbuf_len[i])) 
        { 
          read_start = rbuf_start[i]; 
          read_stride = rbuf_len[i]; 
        } 
      } 
    } 

    in_blk = in + read_start;
    local_iy = iy - read_start;

    read_gpu_suNf_vector(read_stride, sn.c[0], in_blk, local_iy, 0);
    read_gpu_suNf_vector(read_stride, sn.c[1], in_blk, local_iy, 2);
    _vector_i_sub_assign_f(sn.c[0], sn.c[1]);

    read_gpu_suNf(read_stride, u, gauge + 4*read_start, iy-read_start, 3);
    _suNf_theta_Z_inverse_multiply(sn.c[1], u, sn.c[0]);

    _vector_add_assign_f(r.c[0], sn.c[1]);
    _vector_i_add_assign_f(r.c[2], sn.c[1]);

    read_gpu_suNf_vector(read_stride, sn.c[0], in_blk, local_iy, 1);
    read_gpu_suNf_vector(read_stride, sn.c[1], in_blk, local_iy, 3);
    _vector_i_add_assign_f(sn.c[0], sn.c[1]);

    _suNf_theta_Z_inverse_multiply(sn.c[1], u, sn.c[0]);

    _vector_add_assign_f(r.c[1], sn.c[1]);
    _vector_i_sub_assign_f(r.c[3], sn.c[1]);

    __syncthreads();
    /******************************** end of directions *********************************/
    _spinor_mul_f(r, -0.5, r);

    write_gpu_suNf_vector(write_stride, r.c[0], out, local_ix, 0);
    write_gpu_suNf_vector(write_stride, r.c[1], out, local_ix, 1);
    write_gpu_suNf_vector(write_stride, r.c[2], out, local_ix, 2);
    write_gpu_suNf_vector(write_stride, r.c[3], out, local_ix, 3);
  }
}

void (*Dphi_) (spinor_field *out, spinor_field *in)=Dphi_gpu_;
void (*Dphi) (double m0, spinor_field *out, spinor_field *in)=Dphi_gpu;
void (*g5Dphi) (double m0, spinor_field *out, spinor_field *in)=g5Dphi_gpu;
void (*g5Dphi_sq) (double m0, spinor_field *out, spinor_field *in)=g5Dphi_sq_gpu;
unsigned long int (*getMVM) ()=getMVM_gpu;
void (*Dphi_eopre) (double m0, spinor_field *out, spinor_field *in)=Dphi_eopre_gpu;
void (*Dphi_oepre) (double m0, spinor_field *out, spinor_field *in)=Dphi_oepre_gpu;
void (*g5Dphi_eopre) (double m0, spinor_field *out, spinor_field *in)=g5Dphi_eopre_gpu;
void (*g5Dphi_eopre_sq) (double m0, spinor_field *out, spinor_field *in)=g5Dphi_eopre_sq_gpu;

#endif
