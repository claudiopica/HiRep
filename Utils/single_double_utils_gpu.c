/***************************************************************************\
* Copyright (c) 2012, Ulrik Ishøj Søndergaard, Claudio Pica                *   
* All rights reserved.                                                     * 
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


#include "gpu.h"
 
 
 __global__ void assign_f2d_kernel(double* out, float* in, const int N)
 {
     int ix = blockIdx.x*BLOCK_SIZE + threadIdx.x;
     ix = min(ix,N-1);						
     out[ix]=(double)in[ix];
 }
 
 __global__ void assign_d2f_kernel(float* out, double* in, const int N)
 {
     int ix = blockIdx.x*BLOCK_SIZE + threadIdx.x;
     ix = min(ix,N-1);						
     out[ix]=(float)in[ix];
 }


void assign_s2sd(spinor_field *out, spinor_field_flt *in) {
    int size, grid;
    double *r;
    float *rf;
    
    _TWO_SPINORS_MATCHING(in,out);
    
    r=(double*)(_GPU_FIELD_BLK(out,0)); //first block only 
    rf=(float*)(_GPU_FIELD_BLK(in,0));
    size = in->type->master_end[0] -  in->type->master_start[0] + 1;
    size *= sizeof(suNf_spinor)/sizeof(double); ////lenght of the block in real numbers	
    
    grid = size/BLOCK_SIZE + ((size % BLOCK_SIZE == 0) ? 0 : 1); 
    assign_f2d_kernel<<<grid,BLOCK_SIZE>>>(r, rf, size);
    CudaCheckError();	
}

void assign_sd2s(spinor_field_flt *out, spinor_field *in) {
    int size, grid;
    double *r;
    float *rf;
    
    _TWO_SPINORS_MATCHING(in,out);
    
    r=(double*)(_GPU_FIELD_BLK(in,0)); //first block only
    rf=(float*)(_GPU_FIELD_BLK(out,0));
    size = in->type->master_end[0] -  in->type->master_start[0] + 1;
    size *= sizeof(suNf_spinor)/sizeof(double); ////lenght of the block in real numbers	
    
    grid = size/BLOCK_SIZE + ((size % BLOCK_SIZE == 0) ? 0 : 1); 
    assign_d2f_kernel<<<grid,BLOCK_SIZE>>>(rf, r, size);
    CudaCheckError();	
}

void assign_u2ud(){
    int size, grid;
    double *r;
    float *rf;
    
    r=(double*)(_GPU_FIELD_BLK(u_gauge,0)); //first block only
    rf=(float*)(_GPU_FIELD_BLK(u_gauge_flt,0));
    size=glattice.master_end[0] - glattice.master_start[0] + 1;
    size*=4*sizeof(suNg)/sizeof(double); //lenght of the block in real numbers
    
    grid = size/BLOCK_SIZE + ((size % BLOCK_SIZE == 0) ? 0 : 1); 
    assign_f2d_kernel<<<grid,BLOCK_SIZE>>>(r, rf, size);
    CudaCheckError();	
}

void assign_ud2u(){
    int size, grid;
    double *r;
    float *rf;
    
    r=(double*)(_GPU_FIELD_BLK(u_gauge,0)); //first block only
    rf=(float*)(_GPU_FIELD_BLK(u_gauge_flt,0));
    size=glattice.master_end[0] - glattice.master_start[0] + 1;
    size*=4*sizeof(suNg)/sizeof(double); //lenght of the block in real numbers
    
    grid = size/BLOCK_SIZE + ((size % BLOCK_SIZE == 0) ? 0 : 1); 
    assign_d2f_kernel<<<grid,BLOCK_SIZE>>>(rf, r, size);
    CudaCheckError();	
}

void assign_u2ud_f(){
    int size, grid;
    double *r;
    float *rf;
    
    r=(double*)(_GPU_FIELD_BLK(u_gauge_f,0)); //first block only
    rf=(float*)(_GPU_FIELD_BLK(u_gauge_f_flt,0));
    size=glattice.master_end[0] - glattice.master_start[0] + 1;
    size*=4*sizeof(suNf)/sizeof(double); //lenght of the block in real numbers
    
    grid = size/BLOCK_SIZE + ((size % BLOCK_SIZE == 0) ? 0 : 1); 
    assign_f2d_kernel<<<grid,BLOCK_SIZE>>>(r, rf, size);
    CudaCheckError();	
}

void assign_ud2u_f(){
    int size, grid;
    double *r;
    float *rf;
    
    r=(double*)(_GPU_FIELD_BLK(u_gauge_f,0)); //first block only
    rf=(float*)(_GPU_FIELD_BLK(u_gauge_f_flt,0));
    size=glattice.master_end[0] - glattice.master_start[0] + 1;
    size*=4*sizeof(suNf)/sizeof(double); //lenght of the block in real numbers
    
    grid = size/BLOCK_SIZE + ((size % BLOCK_SIZE == 0) ? 0 : 1); 
    assign_d2f_kernel<<<grid,BLOCK_SIZE>>>(rf, r, size);
    CudaCheckError();	
}


#endif //WITH_GPU
