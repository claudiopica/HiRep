/***************************************************************************\
* Copyright (c) 20012, Ulrik Ishøj Søndergaard							*   
* All rights reserved.											       * 
\***************************************************************************/

/*******************************************************************************
*
* File single_double_utils_gpu.c
*
* Functions for conversion from single to double precision and viceversa with gpu
*
*******************************************************************************/
#ifdef WITH_GPU
#include <stdlib.h>
#include "utils.h"
#include "suN.h"
#include "error.h"
#include "global.h"
#include "spinor_field.h"
#include "memory.h"

/*
void assign_s2sd(spinor_field *out, spinor_field_flt *in) {
  spinor_field_copy_from_gpu_f_flt(in);
  assign_s2sd_cpu(out,in);
  spinor_field_copy_to_gpu_f(out);
}

void assign_sd2s(spinor_field_flt *out, spinor_field *in) {
  spinor_field_copy_from_gpu_f(in);
  assign_sd2s_cpu(out,in);
  spinor_field_copy_to_gpu_f_flt(out);
}
*/


#include "gpu.h"
 
 
 __global__ void assign_s2sd_kernel(suNf_spinor* out, suNf_spinor_flt* in, const int N)
 {
 int ix = blockIdx.x*BLOCK_SIZE + threadIdx.x;
 ix = min(ix,N-1);						
 ((double*)(out))[ix]=(double)((float*)(in))[ix];
 }
 
 __global__ void assign_sd2s_kernel(suNf_spinor_flt* out, suNf_spinor* in, const int N)
 {
 int ix = blockIdx.x*BLOCK_SIZE + threadIdx.x;
 ix = min(ix,N-1);						
 ((float*)(out))[ix]=(float)((double*)(in))[ix];
 }


__global__ void assign_u2ud_quaternions_kernel(suNg* gauge, suNg_flt* gauge_flt, int N){ //Only for quaternions
  suNg u;
  suNg_flt u_flt;
  int ix = blockIdx.x*BLOCK_SIZE+ threadIdx.x;
  int i;
  if (ix>=N/2) {
    gauge+=2*N;
    ix -= N/2;
  }
  for (i=0;i<4;++i){
    _suNg_flt_read_gpu(N/2,u_flt,gauge_flt,ix,i);
    u.c[0]=(double) u_flt.c[0];
    u.c[1]=(double) u_flt.c[1];
    u.c[2]=(double) u_flt.c[2];
    u.c[3]=(double) u_flt.c[3];
    _suNg_write_gpu(N/2,u,gauge,ix,i);
  }
}

__global__ void assign_ud2u_quaternions_kernel(suNg* gauge, suNg_flt* gauge_flt, int N){ //Only for quaternions
  suNg u;
  suNg_flt u_flt;
  int ix = blockIdx.x*BLOCK_SIZE+ threadIdx.x;
  int i;
  if (ix>=N/2) {
    gauge+=2*N;
    gauge_flt+=2*N;
    ix -= N/2;
  }
  for (i=0;i<4;++i){
    _suNg_read_gpu(N/2,u,gauge,ix,i);
    u_flt.c[0]=(float) u.c[0];
    u_flt.c[1]=(float) u.c[1];
    u_flt.c[2]=(float) u.c[2];
    u_flt.c[3]=(float) u.c[3];
    _suNg_flt_write_gpu(N/2,u_flt,gauge_flt,ix,i);
  }
}


void assign_s2sd(spinor_field *out, spinor_field_flt *in) {
  int size, grid;
  size = in->type->master_end[0] -  in->type->master_start[0] + 1;
  size *= sizeof(suNf_spinor)/sizeof(double);				
  grid = size/BLOCK_SIZE + ((size % BLOCK_SIZE == 0) ? 0 : 1); 	// NF is defined in suN_types.h
  assign_s2sd_kernel<<<grid,BLOCK_SIZE>>>(START_SP_ADDRESS_GPU(out),START_SP_ADDRESS_GPU(in), size);
  CudaCheckError();	
}

void assign_sd2s(spinor_field_flt *out, spinor_field *in) {
  int size, grid;
  size = in->type->master_end[0] -  in->type->master_start[0] + 1;
  size *= sizeof(suNf_spinor)/sizeof(double);				
  grid = size/BLOCK_SIZE + ((size % BLOCK_SIZE == 0) ? 0 : 1); 	// NF is defined in suN_types.h
  assign_sd2s_kernel<<<grid,BLOCK_SIZE>>>(START_SP_ADDRESS_GPU(out),START_SP_ADDRESS_GPU(in), size);
  CudaCheckError();
}

void assign_u2ud(){
  int N = u_gauge->type->master_end[0] -  u_gauge->type->master_start[0] + 1;
  int grid = N/BLOCK_SIZE + ((N % BLOCK_SIZE == 0) ? 0 : 1);
#ifdef WITH_QUATERNIONS
  assign_u2ud_quaternions_kernel<<<grid,BLOCK_SIZE>>>(u_gauge->gpu_ptr,u_gauge_flt->gpu_ptr,N);
#else
#error : single <-> double : GPU only with quaternions
#endif
}

void assign_ud2u(){
  int N = u_gauge->type->master_end[0] -  u_gauge->type->master_start[0] + 1;
  int grid = N/BLOCK_SIZE + ((N % BLOCK_SIZE == 0) ? 0 : 1);
#ifdef WITH_QUATERNIONS
  assign_ud2u_quaternions_kernel<<<grid,BLOCK_SIZE>>>(u_gauge->gpu_ptr,u_gauge_flt->gpu_ptr,N);
#else
#error : single <-> double : GPU only with quaternions
#endif
}

void assign_u2ud_f(){
#ifdef WITH_QUATERNIONS
  assign_u2ud();
#else
#error : single <-> double : GPU only with quaternions
#endif 
}

void assign_ud2u_f(){
#ifdef WITH_QUATERNIONS
  assign_ud2u();
#else
#error : single <-> double : GPU only with quaternions
#endif 
}


#endif //WITH_GPU
