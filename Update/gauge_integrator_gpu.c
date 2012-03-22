#include "global.h"
#include "suN.h"
#include "utils.h"
#include "update.h"
#include "representation.h"
#include "logger.h"
#include "communications.h"
#include "gpu.h"
#include <assert.h>


#define _suNg_read_gpu(stride,v,in,iy,x)\
iw=(iy)+((x)*4)*(stride);\
(v).c[0]=((double*)(in))[iw]; iw+=(stride); \
(v).c[1]=((double*)(in))[iw]; iw+=(stride);\
(v).c[2]=((double*)(in))[iw]; iw+=(stride);\
(v).c[3]=((double*)(in))[iw]

#define _suNg_write_gpu(stride,v,in,iy,x)\
iw=(iy)+((x)*4)*(stride);\
((double*)(in))[iw]=(v).c[0]; iw+=(stride);\
((double*)(in))[iw]=(v).c[1]; iw+=(stride);\
((double*)(in))[iw]=(v).c[2]; iw+=(stride);\
((double*)(in))[iw]=(v).c[3]

#define _suNg_av_read_gpu(stride,v,in,iy,x)\
iw=(iy)+((x)*3)*(stride);\
(v).c[0]=((double*)(in))[iw]; iw+=(stride); \
(v).c[1]=((double*)(in))[iw]; iw+=(stride);\
(v).c[2]=((double*)(in))[iw]

__global__ void gauge_integrator_gpu_kernel(suNg* gauge, suNg_algebra_vector* momenta, double dt, int N){ //Only for quaternions
  int ix = blockIdx.x*BLOCK_SIZE+ threadIdx.x;
  int iw,i;
  int vol4h = N/2;
  suNg u1,u2,u3;
  suNg_algebra_vector h;
  ix = min(ix,N-1);
  if (ix>=N/2){
    ix -= vol4h;
    gauge += 4*vol4h;
    momenta += 4*vol4h;
  }
  for (i=0;i<4;++i){
    _suNg_av_read_gpu(vol4h,h,momenta,ix,i);
    _suNg_read_gpu(vol4h,u1,gauge,ix,i);
    _suNg_exp(dt,h,u2);
    _suNg_times_suNg(u3,u2,u1);
    _suNg_write_gpu(vol4h,u3,gauge,ix,i);
  }
}


void gauge_integrator_gpu(suNg_av_field *momenta, double dt){
  int N = T*X*Y*Z;//u_gauge->type->master_end[0] -  u_gauge->type->master_start[0] + 1;
  int grid = N/BLOCK_SIZE + ((N % BLOCK_SIZE == 0) ? 0 : 1);
  gauge_integrator_gpu_kernel<<<grid,BLOCK_SIZE>>>(u_gauge->gpu_ptr,momenta->gpu_ptr,dt,N);
}

#undef _suNg_read_gpu
#undef _suNg_write_gpu
#undef _suNg_av_read_gpu
