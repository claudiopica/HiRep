/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *   
* All rights reserved.                                                      * 
\***************************************************************************/

#include "global.h"
#include "observables.h"
#include "update.h"
#include "suN.h"
#include "linear_algebra.h"
#include "error.h"
#include "representation.h"
#include "memory.h"
#include "gpu.h"
#ifdef WITH_GPU

extern rhmc_par _update_par;

__global__ void zero_scalar_field_gpu(double* loc_action, int N){
  int ix = blockIdx.x*BLOCK_SIZE + threadIdx.x;
  ix = min(ix,N-1);
  loc_action[ix]=0.;
}

__global__ void minus_scalar_field_gpu(double* loc_action, int N){
  int ix = blockIdx.x*BLOCK_SIZE + threadIdx.x;
  ix = min(ix,N-1);
  loc_action[ix] = -loc_action[ix];
}

__global__ void local_momenta_gpu(double* loc_action, suNg_algebra_vector* momenta, int N){
  int i;
  int ix = blockIdx.x*BLOCK_SIZE + threadIdx.x;
  double tmp,a=0.;
  ix = min(ix,N-1);
  for (i=0;i<4;++i){
    suNg_algebra_vector av1 = momenta[coord_to_index(ix,i)];
    _algebra_vector_sqnorm_g(tmp,av1);
    a+=tmp;
  }
  a*=0.5*_FUND_NORM2;
  loc_action[ix]+=a;
}


#define _suNf_read_spinor_gpu(stride,v,in,iy,x)\
iz=(iy)+((x)*3)*(stride);\
(v).c[0]=((complex*)(in))[iz]; iz+=(stride); \
(v).c[1]=((complex*)(in))[iz]; iz+=(stride);\
(v).c[2]=((complex*)(in))[iz]

__global__ void  local_pseudo_fermions_gpu(double* loc_action, suNf_spinor* phi1, suNf_spinor* phi2, int stride, int N){
  int ix = blockIdx.x*BLOCK_SIZE + threadIdx.x;
  int iz;
  suNf_spinor s1,s2;
  double tmp;
  ix = min(ix,N-1);
  _suNf_read_spinor_gpu(stride,s1.c[0],phi1,ix,0);
  _suNf_read_spinor_gpu(stride,s1.c[1],phi1,ix,1);
  _suNf_read_spinor_gpu(stride,s1.c[2],phi1,ix,2);
  _suNf_read_spinor_gpu(stride,s1.c[3],phi1,ix,3);
  _suNf_read_spinor_gpu(stride,s2.c[0],phi2,ix,0);
  _suNf_read_spinor_gpu(stride,s2.c[1],phi2,ix,1);
  _suNf_read_spinor_gpu(stride,s2.c[2],phi2,ix,2);
  _suNf_read_spinor_gpu(stride,s2.c[3],phi2,ix,3);
  _spinor_prod_re_f(tmp,s1,s2);
  loc_action[ix]+=tmp;
}


#define _suNg_read_gpu(stride,v,in,iy,x)\
iw=(iy)+((x)*4)*(stride);\
(v).c[0]=((double*)(in))[iw]; iw+=(stride); \
(v).c[1]=((double*)(in))[iw]; iw+=(stride);\
(v).c[2]=((double*)(in))[iw]; iw+=(stride);\
(v).c[3]=((double*)(in))[iw]

__device__ double plaq_gpu(const suNg* gauge, int ix,int mu,int nu, const int *iup, int vol4h)
{
  int iy,iz,iw;
  double p=0;
  suNg v1,v2,v3,v4, w1, w2, w3 ;
  
  
  iy=iup(ix,mu);
  iz=iup(ix,nu);
  
  if (ix<vol4h) {
    iy-=vol4h;
    iz-=vol4h;
    _suNg_read_gpu(vol4h,v1,gauge,ix,mu);
    _suNg_read_gpu(vol4h,v2,gauge+4*vol4h,iy,nu);
    _suNg_read_gpu(vol4h,v3,gauge+4*vol4h,iz,mu);
    _suNg_read_gpu(vol4h,v4,gauge,ix,nu);
  }
  else {
    ix-=vol4h;
    _suNg_read_gpu(vol4h,v1,gauge+4*vol4h,ix,mu);
    _suNg_read_gpu(vol4h,v2,gauge,iy,nu);
    _suNg_read_gpu(vol4h,v3,gauge,iz,mu);
    _suNg_read_gpu(vol4h,v4,gauge+4*vol4h,ix,nu);
  }
  
  _suNg_times_suNg(w1,v1,v2);
  _suNg_times_suNg(w2,v4,v3);
  _suNg_times_suNg_dagger(w3,w1,w2);      
  
  _suNg_trace_re(p,w3);
  
  return p;
}

__global__ void local_gauge_action_gpu(const suNg* gauge,double *loc_action, const int *iup, const int vol4h, double beta)
{
	double pa; 
	int ix = blockIdx.x*blockDim.x+ threadIdx.x;
	ix = min(ix,vol4h*2-1);
	
	pa=plaq_gpu(gauge, ix,1,0,iup,vol4h);
	pa+=plaq_gpu(gauge, ix,2,0,iup,vol4h);
	pa+=plaq_gpu(gauge, ix,2,1,iup,vol4h);
	pa+=plaq_gpu(gauge, ix,3,0,iup,vol4h);
	pa+=plaq_gpu(gauge, ix,3,1,iup,vol4h);
	pa+=plaq_gpu(gauge, ix,3,2,iup,vol4h);
	
	loc_action[ix]+=-pa*beta/(double)NG;	
}

void local_gauge_action_gpu_caller(scalar_field *loc_action, double beta){
	const int vol4=T*X*Y*Z;
	int grid = vol4/BLOCK_SIZE + ((vol4 % BLOCK_SIZE == 0) ? 0 : 1);
	
	local_gauge_action_gpu<<<grid,BLOCK_SIZE>>>(u_gauge->gpu_ptr,loc_action->gpu_ptr, iup_gpu, vol4/2, beta);
	CudaCheckError();
}

double scalar_field_sum(scalar_field* sf){
  int N = sf->type->master_end[0] -  sf->type->master_start[0] + 1;  
  return global_sum_gpu(START_SP_ADDRESS_GPU(sf),N);
}

/*
 * compute the local action at every site for the HMC (Done at GPU)
 * H = | momenta |^2 + S_g + < phi1, phi2>
 */

void local_hmc_action(local_action_type type,
                      scalar_field *loc_action,
                      suNg_av_field *momenta,
                      spinor_field *phi1,
                      spinor_field *phi2) { 

  unsigned int N, grid, N_sp, grid_sp,i;
	
  gfield_copy_to_gpu(u_gauge);
  suNg_av_field_copy_to_gpu(momenta);
	
  /* check input types */
#ifndef CHECK_SPINOR_MATCHING
  _TWO_SPINORS_MATCHING(u_gauge,loc_action); /* check that action is defined on the global lattice */
  _TWO_SPINORS_MATCHING(loc_action,momenta);
  _TWO_SPINORS_MATCHING(loc_action,phi1);
  _TWO_SPINORS_MATCHING(loc_action,phi2);
#endif

  N = loc_action->type->master_end[0] -  loc_action->type->master_start[0] + 1;
  lprintf("LOCAL_ACTION_GPU",0,"Volume 1 %d, volume 2 %d\n",N,T*X*Y*Z);  
  grid = N/BLOCK_SIZE + ((N % BLOCK_SIZE == 0) ? 0 : 1);
  switch(type) {
  case NEW:
    {
      zero_scalar_field_gpu<<<grid,BLOCK_SIZE>>>(loc_action->gpu_ptr,N);
    }
    break;
  case DELTA:
    {
      minus_scalar_field_gpu<<<grid,BLOCK_SIZE>>>(loc_action->gpu_ptr,N);
    }
    break;
  default:
    error(1,1,"local_hmc_action","Invalid type");
  }
  //  lprintf("LOCAL_ACTION_GPU",0,"After init: scalar_field_sum: %1.10g\n",scalar_field_sum(loc_action));

  //Momenta
  local_momenta_gpu<<<grid,BLOCK_SIZE>>>(loc_action->gpu_ptr, momenta->gpu_ptr , N);
  
#ifndef ROTATED_SF
  local_gauge_action_gpu<<<grid,BLOCK_SIZE>>>(u_gauge->gpu_ptr,loc_action->gpu_ptr, iup_gpu, N/2, _update_par.beta);
  CudaCheckError();
  
#else /* ROTATED_SF */
  "ROTATED_SF not yet implemented for GPU"
#endif /* ROTATED_SF */

  
  N_sp = phi1[0].type->master_end[0] -  phi1[0].type->master_start[0] + 1; 
  //N_sp *= sizeof(suNf_spinor)/sizeof(complex);
  grid_sp = N_sp/BLOCK_SIZE + ((N_sp % BLOCK_SIZE == 0) ? 0 : 1);  
  /* pseudofermion fields can be defined only on even sites is the preconditioning is used */
  
  /* Fermions */
  for (i=0;i<_update_par.n_pf;++i){
    local_pseudo_fermions_gpu<<<grid_sp,BLOCK_SIZE>>>( loc_action->gpu_ptr, 
    						       START_SP_ADDRESS_GPU(&phi1[i]), 
						       START_SP_ADDRESS_GPU(&phi2[i]),N_sp,N_sp);
  }
  //  lprintf("LOCAL_ACTION_GPU",0,"At end: scalar_field_sum: %1.10g\n",scalar_field_sum(loc_action));
   
}

#endif
