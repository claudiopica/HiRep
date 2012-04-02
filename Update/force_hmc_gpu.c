#ifdef WITH_GPU

#include "global.h"
#include "update.h"
#include "suN.h"
#include "utils.h"
#include "representation.h"
#include "logger.h"
#include "linear_algebra.h"
#include "memory.h"
#include "communications.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>


#define _F_DIR0(u,chi1,chi2)				      \
  _vector_add_f(ptmp,(chi2).c[0],(chi2).c[2]);	      \
  _suNf_multiply(p.c[0],pu,ptmp);	      \
  _vector_add_f(ptmp,(chi2).c[1],(chi2).c[3]);	      \
  _suNf_multiply(p.c[1],pu,ptmp);	      \
  _vector_sub_f(p.c[2],(chi1).c[0],(chi1).c[2]);	      \
  _vector_sub_f(p.c[3],(chi1).c[1],(chi1).c[3]);	      \
  _suNf_FMAT((u),p)

#define _F_DIR1(u,chi1,chi2)				      \
  _vector_i_add_f(ptmp,(chi2).c[0],(chi2).c[3]);	      \
  _suNf_multiply(p.c[0],pu,ptmp);	      \
  _vector_i_add_f(ptmp,(chi2).c[1],(chi2).c[2]);	      \
  _suNf_multiply(p.c[1],pu,ptmp);	      \
  _vector_i_sub_f(p.c[2],(chi1).c[0],(chi1).c[3]);	      \
  _vector_i_sub_f(p.c[3],(chi1).c[1],(chi1).c[2]);	      \
  _suNf_FMAT((u),p)

#define _F_DIR2(u,chi1,chi2)				      \
  _vector_add_f(ptmp,(chi2).c[0],(chi2).c[3]);	      \
  _suNf_multiply(p.c[0],pu,ptmp);	      \
  _vector_sub_f(ptmp,(chi2).c[1],(chi2).c[2]);	      \
  _suNf_multiply(p.c[1],pu,ptmp);	      \
  _vector_sub_f(p.c[2],(chi1).c[0],(chi1).c[3]);	      \
  _vector_add_f(p.c[3],(chi1).c[1],(chi1).c[2]);	      \
  _suNf_FMAT((u),p)

#define _F_DIR3(u,chi1,chi2)				      \
  _vector_i_add_f(ptmp,(chi2).c[0],(chi2).c[2]);	      \
  _suNf_multiply(p.c[0],pu,ptmp);	      \
  _vector_i_sub_f(ptmp,(chi2).c[1],(chi2).c[3]);	      \
  _suNf_multiply(p.c[1],pu,ptmp);	      \
  _vector_i_sub_f(p.c[2],(chi1).c[0],(chi1).c[2]);	      \
  _vector_i_add_f(p.c[3],(chi1).c[1],(chi1).c[3]);	      \
  _suNf_FMAT((u),p)

#define _suNf_read_spinor_gpu(stride,v,in,iy,x)\
iw=(iy)+((x)*3)*(stride);\
(v).c[0]=((complex*)(in))[iw]; iw+=(stride); \
(v).c[1]=((complex*)(in))[iw]; iw+=(stride);\
(v).c[2]=((complex*)(in))[iw]

#define _suNf_read_full_spinor_gpu(stride,sp,in,iy)\
  _suNf_read_spinor_gpu(stride,(sp).c[0],(in),iy,0);	\
  _suNf_read_spinor_gpu(stride,(sp).c[1],(in),iy,1);	\
  _suNf_read_spinor_gpu(stride,(sp).c[2],(in),iy,2);	\
  _suNf_read_spinor_gpu(stride,(sp).c[3],(in),iy,3);

#define _suNg_av_read_gpu(stride,v,in,iy,x)\
iw=(iy)+((x)*3)*(stride);\
(v).c[0]=((double*)(in))[iw]; iw+=(stride); \
(v).c[1]=((double*)(in))[iw]; iw+=(stride);\
(v).c[2]=((double*)(in))[iw]

#define _suNg_av_write_gpu(stride,v,out,iy,x)\
iw=(iy)+((x)*3)*(stride);\
((double*)(out))[iw]=(v).c[0]; iw+=(stride); \
((double*)(out))[iw]=(v).c[1]; iw+=(stride);\
((double*)(out))[iw]=(v).c[2]

#define _suNf_read_gpu(stride,v,in,iy,x)\
iw=(iy)+((x)*4)*(stride);\
(v).c[0]=((double*)(in))[iw]; iw+=(stride);\
(v).c[1]=((double*)(in))[iw]; iw+=(stride);\
(v).c[2]=((double*)(in))[iw]; iw+=(stride);\
(v).c[3]=((double*)(in))[iw]

#define _algebra_vector_mul_add_assign_gpu_g(stride,v,iy,x,r,in)\
iw=(iy)+((x)*3)*(stride);\
((double*)(v))[iw]+=(in).c[0]*(r); iw+=(stride); \
((double*)(v))[iw]+=(in).c[1]*(r); iw+=(stride);\
((double*)(v))[iw]+=(in).c[2]*(r)


__global__ void force_hmc_gpu_kernel(suNg_algebra_vector* force, suNf_spinor *Xs, suNf_spinor *Ys, suNf *gauge, double fs, const int *iup, int N){
  suNf_spinor chi1,chi2,p;
  suNf_vector ptmp;
  suNg_algebra_vector f;
  suNf_spinor *Xss, *Yss;
  int iw,iy,ig;
  int ix = blockIdx.x*BLOCK_SIZE+ threadIdx.x;
  int vol4h = N/2;
  int shift, shift2;
  suNf_FMAT s1;
  suNf pu;
  ix = min(ix,N-1);
  
  shift  = (ix<vol4h)? 0 : vol4h; //Even ? 0 : vol4h
  ig = ix-shift; //always in [0 : vol4h -1]. Normalized ix
  shift2 = vol4h-shift; //shift factor for iy, with opposite parity w.r.t. ix
  force += 4*shift; //normalized force
  gauge += 4*shift; //and gauge

  Xss=Xs+shift2;
  Xs+=shift;
  Yss=Ys+shift2;
  Ys+=shift;

  //Dir 0
  iy = iup(ix,0) - shift2;
  _suNf_read_full_spinor_gpu(vol4h,chi1,Xs,ig);
  _suNf_read_full_spinor_gpu(vol4h,chi2,Yss,iy);
  _suNf_read_gpu(vol4h,pu,gauge,ig,0);
  _suNf_FMAT_zero(s1);
  _F_DIR0(s1,chi1,chi2);
  _suNf_read_full_spinor_gpu(vol4h,chi1,Ys,ig);
  _suNf_read_full_spinor_gpu(vol4h,chi2,Xss,iy);
  _F_DIR0(s1,chi1,chi2);

  _algebra_project(f,s1);
  _algebra_vector_mul_add_assign_gpu_g(vol4h,force,ig,0,fs,f);

  //Dir 1
  iy = iup(ix,1) - shift2;
  _suNf_read_full_spinor_gpu(vol4h,chi1,Xs,ix);
  _suNf_read_full_spinor_gpu(vol4h,chi2,Yss,iy);
  _suNf_read_gpu(vol4h,pu,gauge,ig,1);
  _suNf_FMAT_zero(s1);
  _F_DIR1(s1,chi1,chi2);
  _suNf_read_full_spinor_gpu(vol4h,chi1,Ys,ix);
  _suNf_read_full_spinor_gpu(vol4h,chi2,Xss,iy);
  _F_DIR1(s1,chi1,chi2);

  _algebra_project(f,s1);
  _algebra_vector_mul_add_assign_gpu_g(vol4h,force,ig,1,fs,f);

  //Dir 2
  iy = iup(ix,2) - shift2;
  _suNf_read_full_spinor_gpu(vol4h,chi1,Xs,ix);
  _suNf_read_full_spinor_gpu(vol4h,chi2,Yss,iy);
  _suNf_read_gpu(vol4h,pu,gauge,ig,2);
  _suNf_FMAT_zero(s1);
  _F_DIR2(s1,chi1,chi2);
  _suNf_read_full_spinor_gpu(vol4h,chi1,Ys,ix);
  _suNf_read_full_spinor_gpu(vol4h,chi2,Xss,iy);
  _F_DIR2(s1,chi1,chi2);

  _algebra_project(f,s1);
  _algebra_vector_mul_add_assign_gpu_g(vol4h,force,ig,2,fs,f);

  //Dir 3
  iy = iup(ix,3) - shift2;
  _suNf_read_full_spinor_gpu(vol4h,chi1,Xs,ix);
  _suNf_read_full_spinor_gpu(vol4h,chi2,Yss,iy);
  _suNf_read_gpu(vol4h,pu,gauge,ig,3);
  _suNf_FMAT_zero(s1);
  _F_DIR3(s1,chi1,chi2);
  _suNf_read_full_spinor_gpu(vol4h,chi1,Ys,ix);
  _suNf_read_full_spinor_gpu(vol4h,chi2,Xss,iy);
  _F_DIR3(s1,chi1,chi2);

  _algebra_project(f,s1);
  _algebra_vector_mul_add_assign_gpu_g(vol4h,force,ig,3,fs,f);

}


void force_hmc_gpu(suNg_av_field* force, spinor_field *Xs, spinor_field *Ys, double dt, force_hmc_par *par){
  int N = T*X*Y*Z;//u_gauge->type->master_end[0] -  u_gauge->type->master_start[0] + 1;
  int grid = N/BLOCK_SIZE + ((N % BLOCK_SIZE == 0) ? 0 : 1);
  double fs;
  suNg_av_field_copy_to_gpu(force);
#ifdef UPDATE_EO
  if(par->hasenbusch != 1) fs = -dt*(_REPR_NORM2/_FUND_NORM2);
  else fs = -par->bD*dt*(_REPR_NORM2/_FUND_NORM2);
#else
  if(par->hasenbusch != 1) fs = dt*(_REPR_NORM2/_FUND_NORM2);
  else fs = par->bD*dt*(_REPR_NORM2/_FUND_NORM2);
#endif
//  force_hmc_gpu_kernel<<<grid,BLOCK_SIZE>>>(force->gpu_ptr,START_SP_ADDRESS_GPU(Xs),START_SP_ADDRESS_GPU(Ys),u_gauge_f->gpu_ptr,fs,iup_gpu,N);

  force_hmc_gpu_kernel<<<grid,BLOCK_SIZE>>>(force->gpu_ptr,Xs->gpu_ptr,Ys->gpu_ptr,u_gauge_f->gpu_ptr,fs,iup_gpu,N);
  CudaCheckError();
  suNg_av_field_copy_from_gpu(force);
}

#undef _suNf_read_spinor_gpu
#undef _suNf_read_full_spinor_gpu
#undef _suNg_av_read_gpu
#undef _suNg_av_write_gpu
#undef _suNf_read_gpu


#undef _F_DIR0
#undef _F_DIR1
#undef _F_DIR2
#undef _F_DIR3

#endif //WITH_GPU
