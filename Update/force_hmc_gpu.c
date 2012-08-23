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

__constant__ complex eitheta_f_gpu[4];

#ifdef BC_T_THETA
#define _T_theta_mulc(r) _vector_mulc_f(ptmp,eitheta_f_gpu[0],(r)); (r)=ptmp
#else
#define _T_theta_mulc(r)
#endif
#ifdef BC_X_THETA
#define _X_theta_mulc(r) _vector_mulc_f(ptmp,eitheta_f_gpu[1],(r)); (r)=ptmp
#else
#define _X_theta_mulc(r)
#endif
#ifdef BC_Y_THETA
#define _Y_theta_mulc(r) _vector_mulc_f(ptmp,eitheta_f_gpu[2],(r)); (r)=ptmp
#else
#define _Y_theta_mulc(r)
#endif
#ifdef BC_Z_THETA
#define _Z_theta_mulc(r) _vector_mulc_f(ptmp,eitheta_f_gpu[3],(r)); (r)=ptmp
#else
#define _Z_theta_mulc(r)
#endif

#define _F_DIR0(u,chi1,chi2)				      \
  _vector_add_f(ptmp,(chi2).c[0],(chi2).c[2]);		      \
  _suNf_multiply(p.c[0],pu,ptmp);			      \
  _T_theta_mulc(p.c[0]);				      \
  _vector_add_f(ptmp,(chi2).c[1],(chi2).c[3]);		      \
  _suNf_multiply(p.c[1],pu,ptmp);			      \
  _T_theta_mulc(p.c[1]);                                      \
  _vector_sub_f(p.c[2],(chi1).c[0],(chi1).c[2]);	      \
  _vector_sub_f(p.c[3],(chi1).c[1],(chi1).c[3]);	      \
  _suNf_FMAT((u),p)

#define _F_DIR1(u,chi1,chi2)				      \
  _vector_i_add_f(ptmp,(chi2).c[0],(chi2).c[3]);	      \
  _suNf_multiply(p.c[0],pu,ptmp);			      \
  _X_theta_mulc(p.c[0]);                                      \
  _vector_i_add_f(ptmp,(chi2).c[1],(chi2).c[2]);	      \
  _suNf_multiply(p.c[1],pu,ptmp);			      \
  _X_theta_mulc(p.c[1]);                                      \
  _vector_i_sub_f(p.c[2],(chi1).c[0],(chi1).c[3]);	      \
  _vector_i_sub_f(p.c[3],(chi1).c[1],(chi1).c[2]);	      \
  _suNf_FMAT((u),p)

#define _F_DIR2(u,chi1,chi2)				      \
  _vector_add_f(ptmp,(chi2).c[0],(chi2).c[3]);		      \
  _suNf_multiply(p.c[0],pu,ptmp);			      \
  _Y_theta_mulc(p.c[0]);                                      \
  _vector_sub_f(ptmp,(chi2).c[1],(chi2).c[2]);		      \
  _suNf_multiply(p.c[1],pu,ptmp);			      \
  _Y_theta_mulc(p.c[1]);                                      \
  _vector_sub_f(p.c[2],(chi1).c[0],(chi1).c[3]);	      \
  _vector_add_f(p.c[3],(chi1).c[1],(chi1).c[2]);	      \
  _suNf_FMAT((u),p)

#define _F_DIR3(u,chi1,chi2)				      \
  _vector_i_add_f(ptmp,(chi2).c[0],(chi2).c[2]);	      \
  _suNf_multiply(p.c[0],pu,ptmp);			      \
  _Z_theta_mulc(p.c[0]);                                      \
  _vector_i_sub_f(ptmp,(chi2).c[1],(chi2).c[3]);	      \
  _suNf_multiply(p.c[1],pu,ptmp);			      \
  _Z_theta_mulc(p.c[1]);                                      \
  _vector_i_sub_f(p.c[2],(chi1).c[0],(chi1).c[2]);	      \
  _vector_i_add_f(p.c[3],(chi1).c[1],(chi1).c[3]);	      \
  _suNf_FMAT((u),p)


#define _suNf_read_full_spinor_gpu(stride,sp,in,iy)\
  _suNf_read_spinor_gpu(stride,(sp).c[0],(in),iy,0);	\
  _suNf_read_spinor_gpu(stride,(sp).c[1],(in),iy,1);	\
  _suNf_read_spinor_gpu(stride,(sp).c[2],(in),iy,2);	\
  _suNf_read_spinor_gpu(stride,(sp).c[3],(in),iy,3);


__global__ void force_hmc_gpu_kernel(suNg_algebra_vector* force, suNf_spinor *Xs, suNf_spinor *Ys, suNf *gauge, double fs, const int *iup, int N){
  suNf_spinor chi1,chi2,p;
  suNf_vector ptmp;
  suNg_algebra_vector f;
  suNf_spinor *Xss, *Yss;
  int iy,ig;
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
  _suNf_read_full_spinor_gpu(vol4h,chi1,Xs,ig);
  _suNf_read_full_spinor_gpu(vol4h,chi2,Yss,iy);
  _suNf_read_gpu(vol4h,pu,gauge,ig,1);
  _suNf_FMAT_zero(s1);
  _F_DIR1(s1,chi1,chi2);
  _suNf_read_full_spinor_gpu(vol4h,chi1,Ys,ig);
  _suNf_read_full_spinor_gpu(vol4h,chi2,Xss,iy);
  _F_DIR1(s1,chi1,chi2);

  _algebra_project(f,s1);
  _algebra_vector_mul_add_assign_gpu_g(vol4h,force,ig,1,fs,f);

  //Dir 2
  iy = iup(ix,2) - shift2;
  _suNf_read_full_spinor_gpu(vol4h,chi1,Xs,ig);
  _suNf_read_full_spinor_gpu(vol4h,chi2,Yss,iy);
  _suNf_read_gpu(vol4h,pu,gauge,ig,2);
  _suNf_FMAT_zero(s1);
  _F_DIR2(s1,chi1,chi2);
  _suNf_read_full_spinor_gpu(vol4h,chi1,Ys,ig);
  _suNf_read_full_spinor_gpu(vol4h,chi2,Xss,iy);
  _F_DIR2(s1,chi1,chi2);

  _algebra_project(f,s1);
  _algebra_vector_mul_add_assign_gpu_g(vol4h,force,ig,2,fs,f);

  //Dir 3
  iy = iup(ix,3) - shift2;
  _suNf_read_full_spinor_gpu(vol4h,chi1,Xs,ig);
  _suNf_read_full_spinor_gpu(vol4h,chi2,Yss,iy);
  _suNf_read_gpu(vol4h,pu,gauge,ig,3);
  _suNf_FMAT_zero(s1);
  _F_DIR3(s1,chi1,chi2);
  _suNf_read_full_spinor_gpu(vol4h,chi1,Ys,ig);
  _suNf_read_full_spinor_gpu(vol4h,chi2,Xss,iy);
  _F_DIR3(s1,chi1,chi2);

  _algebra_project(f,s1);
  _algebra_vector_mul_add_assign_gpu_g(vol4h,force,ig,3,fs,f);

}

static void init_bc_gpu(){
#ifdef FERMION_THETA
  static int initialized=0;
  if (!initialized){
    cudaMemcpyToSymbol("eitheta_f_gpu", eitheta, 4*sizeof(complex), 0, cudaMemcpyHostToDevice);
    initialized=1;
  }
#endif
}

void force_hmc_gpu(suNg_av_field* force, spinor_field *Xs, spinor_field *Ys, double dfs){
  int N = T*X*Y*Z;//u_gauge->type->master_end[0] -  u_gauge->type->master_start[0] + 1;
  int grid = N/BLOCK_SIZE + ((N % BLOCK_SIZE == 0) ? 0 : 1);
  init_bc_gpu();
  force_hmc_gpu_kernel<<<grid,BLOCK_SIZE>>>(force->gpu_ptr,Xs->gpu_ptr,Ys->gpu_ptr,u_gauge_f->gpu_ptr,dfs,iup_gpu,N);
  CudaCheckError();
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
