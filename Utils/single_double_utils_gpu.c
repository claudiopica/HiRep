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
 ix = min(ix,N-1);						// NF is defined in suN_types.h
 ((double*)(out))[ix]=(double)((float*)(in))[ix];
 }
 
 __global__ void assign_sd2s_kernel(suNf_spinor_flt* out, suNf_spinor* in, const int N)
 {
 int ix = blockIdx.x*BLOCK_SIZE + threadIdx.x;
 ix = min(ix,N-1);						// NF is defined in suN_types.h
 ((float*)(out))[ix]=(float)((double*)(in))[ix];
 }

void assign_s2sd(spinor_field *out, spinor_field_flt *in) {

	 int size, grid;
	 size=4*2*NF*T*X*Y*Z;				// 4 spinor indices, 2 reals in complex, NF gauge indices
	 grid = size/BLOCK_SIZE + ((size % BLOCK_SIZE == 0) ? 0 : 1); 	// NF is defined in suN_types.h
	 assign_s2sd_kernel<<<grid,BLOCK_SIZE>>>(START_SP_ADDRESS_GPU(out),START_SP_ADDRESS_GPU(in), size);
	 CudaCheckError();	
}

void assign_sd2s(spinor_field_flt *out, spinor_field *in) {

	 int size, grid;
	 size=4*2*NF*T*X*Y*Z;
	 grid = size/BLOCK_SIZE + ((size % BLOCK_SIZE == 0) ? 0 : 1); 	// NF is defined in suN_types.h
	 assign_sd2s_kernel<<<grid,BLOCK_SIZE>>>(START_SP_ADDRESS_GPU(out),START_SP_ADDRESS_GPU(in), size);
	 CudaCheckError();
	 
}
