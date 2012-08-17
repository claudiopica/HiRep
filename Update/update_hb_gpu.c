/***************************************************************************\
* Copyright (c) 2008, Ari Hietanen                                         *   
\***************************************************************************/

/*******************************************************************************
*
* File update.c
*
* Update programs
*
*******************************************************************************/
#ifdef WITH_GPU
#define PROJECT_INTERVAL 10

#include "suN.h"
#include "utils.h"
#include "global.h"
#include "update.h"
#include "communications.h"
#include "gpu.h"

__device__ void normalize_gpu(suNg_vector *v){
   double fact;
   _vector_prod_re_g(fact,*v,*v);
   fact=1.0f/sqrt(fact);
   _vector_mul_g(*v, fact, *v);
}


__global__ void project_gauge_field_gpu(suNg* gauge, int N){ 
#ifdef WITH_QUATERNIONS
  int ix = blockIdx.x*BLOCK_SIZE+ threadIdx.x;
  int i;
  double norm;
  suNg u;
  ix = min(ix,N);
  if (ix>=N/2) {
    gauge+=2*N;
    ix -= N/2;
  }
  for (i=0;i<4;++i){
    _suNg_read_gpu(N/2,u,gauge,ix,i);
    _suNg_sqnorm(norm,u);
    norm = 1./sqrt(0.5*norm);
    _suNg_mul(u,norm,u);
    _suNg_write_gpu(N/2,u,gauge,ix,i);
  }
#else //WITH_QUATERNIONS
  int ix = blockIdx.x*BLOCK_SIZE+ threadIdx.x;
  int d,i,j;
  suNg_vector *v1,*v2;
  complex z;
  suNg u;
  ix = min(ix,N);
  if (ix>=N/2) {
    gauge+=2*N;
    ix -= N/2;
  }
  for (d=0;d<4;++d){
    _suNg_read_gpu(N/2,u,gauge,ix,d);
    v1 = (suNg_vector*)(&u);
    v2= v1+1;
    normalize_gpu(v1);
    for (i=1;i<NG;++i){
      for (j=i;j>0;--j){
	_vector_prod_re_g(z.re,*v1, *v2);
	_vector_prod_im_g(z.im,*v1, *v2);
	_vector_project_g(*v2, z, *v1);
	++v1;
      }
      normalize_gpu(v2);
      ++v2;
      v1=(suNg_vector*)(&u);
    }
    _suNg_write_gpu(N/2,u,gauge,ix,d);
  }

#endif //WITH_QUATERNIONS
}

void project_gauge_field(void){
  int N = u_gauge->type->master_end[0] -  u_gauge->type->master_start[0] + 1;
  int grid = (N-1)/BLOCK_SIZE + 1;
  project_gauge_field_gpu<<<grid,BLOCK_SIZE>>>(u_gauge->gpu_ptr,N);
  //  start_gf_sendrecv(u_gauge);
}



#endif
