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
#include "communications.h"
#include "memory.h"
#include "gpu.h"
#include "hr_complex.h"
#include <iostream>

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

unsigned long int getMVM() {
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

/* Takes an even input spinor and returns an odd spinor */
__global__ void Dphi_gpu_oe(suNf_spinor* __restrict__ out, const suNf_spinor* __restrict__ in,
                            const suNf* __restrict__ gauge, const int* __restrict__ iup_d, const int* __restrict__ idn_d,
                            const int vol4h)
{
  suNf_spinor r;
  suNf_hspinor sn;
  suNf u;
#ifdef FERMION_THETA
  suNf_vector vtmp;
#endif

  int iy;
  int ix = blockIdx.x*BLOCK_SIZE + threadIdx.x;

  if (ix < vol4h) {
    /******************************* direction +0 *********************************/
    iy=iup_d[4*(ix+vol4h)];

    _suNf_read_spinor_gpu(vol4h, sn.c[0], in, iy, 0);
    _suNf_read_spinor_gpu(vol4h, sn.c[1], in, iy, 2);
    _suNf_read_gpu(vol4h, u, gauge, ix+16*vol4h, 0);


    _vector_add_assign_f(sn.c[0], sn.c[1]);
    _suNf_theta_T_multiply(r.c[0], u, sn.c[0]);

    r.c[2]=r.c[0];

    _suNf_read_spinor_gpu(vol4h, sn.c[0], in, iy, 1);
    _suNf_read_spinor_gpu(vol4h, sn.c[1], in, iy, 3);
    _vector_add_assign_f(sn.c[0], sn.c[1]);

    _suNf_theta_T_multiply(r.c[1], u, sn.c[0]);

    r.c[3]=r.c[1];

    __syncthreads();
    /******************************* direction -0 *********************************/
    iy=idn_d[4*(ix+vol4h)];

    _suNf_read_spinor_gpu(vol4h, sn.c[0], in, iy, 0);
    _suNf_read_spinor_gpu(vol4h, sn.c[1], in, iy, 2);
    _suNf_read_gpu(vol4h, u, gauge, iy, 0);

    _vector_sub_assign_f(sn.c[0], sn.c[1]);
    _suNf_theta_T_inverse_multiply(sn.c[1], u, sn.c[0]);

    _vector_add_assign_f(r.c[0], sn.c[1]);
    _vector_sub_assign_f(r.c[2], sn.c[1]);

    _suNf_read_spinor_gpu(vol4h, sn.c[0], in, iy, 1);
    _suNf_read_spinor_gpu(vol4h, sn.c[1], in, iy, 3);
    _vector_sub_assign_f(sn.c[0], sn.c[1]);

    _suNf_theta_T_inverse_multiply(sn.c[1], u, sn.c[0]);

    _vector_add_assign_f(r.c[1], sn.c[1]);
    _vector_sub_assign_f(r.c[3], sn.c[1]);

    __syncthreads();
    /******************************* direction +1 *********************************/
    iy=iup_d[4*(ix+vol4h)+1];

    _suNf_read_spinor_gpu(vol4h, sn.c[0], in, iy, 0);
    _suNf_read_spinor_gpu(vol4h, sn.c[1], in, iy, 3);
    _suNf_read_gpu(vol4h, u, gauge, ix+16*vol4h, 1);

    _vector_i_add_assign_f(sn.c[0], sn.c[1]);
    _suNf_theta_X_multiply(sn.c[1], u, sn.c[0]);

    _vector_add_assign_f(r.c[0], sn.c[1]);
    _vector_i_sub_assign_f(r.c[3], sn.c[1]);

    _suNf_read_spinor_gpu(vol4h, sn.c[0], in, iy, 1);
    _suNf_read_spinor_gpu(vol4h, sn.c[1], in, iy, 2);
    _vector_i_add_assign_f(sn.c[0], sn.c[1]);

    _suNf_theta_X_multiply(sn.c[1], u, sn.c[0]);

    _vector_add_assign_f(r.c[1], sn.c[1]);
    _vector_i_sub_assign_f(r.c[2], sn.c[1]);

    __syncthreads();
    /******************************* direction -1 *********************************/
    iy=idn_d[4*(ix+vol4h)+1];

    _suNf_read_spinor_gpu(vol4h, sn.c[0], in, iy, 0);
    _suNf_read_spinor_gpu(vol4h, sn.c[1], in, iy, 3);
    _suNf_read_gpu(vol4h, u, gauge, iy, 1);

    _vector_i_sub_assign_f(sn.c[0], sn.c[1]);
    _suNf_theta_X_inverse_multiply(sn.c[1], u, sn.c[0]);

    _vector_add_assign_f(r.c[0], sn.c[1]);
    _vector_i_add_assign_f(r.c[3], sn.c[1]);

    _suNf_read_spinor_gpu(vol4h, sn.c[0], in, iy, 1);
    _suNf_read_spinor_gpu(vol4h, sn.c[1], in, iy, 2);
    _vector_i_sub_assign_f(sn.c[0], sn.c[1]);

    _suNf_theta_X_inverse_multiply(sn.c[1], u, sn.c[0]);

    _vector_add_assign_f(r.c[1], sn.c[1]);
    _vector_i_add_assign_f(r.c[2], sn.c[1]);

    __syncthreads();
    /******************************* direction +2 *********************************/
    iy=iup_d[4*(ix+vol4h)+2];

    _suNf_read_spinor_gpu(vol4h, sn.c[0], in, iy, 0);
    _suNf_read_spinor_gpu(vol4h, sn.c[1], in, iy, 3);
    _vector_add_assign_f(sn.c[0], sn.c[1]);

    _suNf_read_gpu(vol4h, u, gauge, ix+16*vol4h, 2);
    _suNf_theta_Y_multiply(sn.c[1], u, sn.c[0]);

    _vector_add_assign_f(r.c[0], sn.c[1]);
    _vector_add_assign_f(r.c[3], sn.c[1]);

    _suNf_read_spinor_gpu(vol4h, sn.c[0], in, iy, 1);
    _suNf_read_spinor_gpu(vol4h, sn.c[1], in, iy, 2);
    _vector_sub_assign_f(sn.c[0], sn.c[1]);

    _suNf_theta_Y_multiply(sn.c[1], u, sn.c[0]);

    _vector_add_assign_f(r.c[1], sn.c[1]);
    _vector_sub_assign_f(r.c[2], sn.c[1]);

    __syncthreads();
    /******************************* direction -2 *********************************/
    iy=idn_d[4*(ix+vol4h)+2];

    _suNf_read_spinor_gpu(vol4h, sn.c[0], in, iy, 0);
    _suNf_read_spinor_gpu(vol4h, sn.c[1], in, iy, 3);
    _vector_sub_assign_f(sn.c[0], sn.c[1]);

    _suNf_read_gpu(vol4h, u, gauge, iy, 2);
    _suNf_theta_Y_inverse_multiply(sn.c[1], u, sn.c[0]);

    _vector_add_assign_f(r.c[0], sn.c[1]);
    _vector_sub_assign_f(r.c[3], sn.c[1]);

    _suNf_read_spinor_gpu(vol4h, sn.c[0], in, iy, 1);
    _suNf_read_spinor_gpu(vol4h, sn.c[1], in, iy, 2);
    _vector_add_assign_f(sn.c[0], sn.c[1]);

    _suNf_theta_Y_inverse_multiply(sn.c[1], u, sn.c[0]);

    _vector_add_assign_f(r.c[1], sn.c[1]);
    _vector_add_assign_f(r.c[2], sn.c[1]);

    __syncthreads();
    /******************************* direction +3 *********************************/
    iy=iup_d[4*(ix+vol4h)+3];

    _suNf_read_spinor_gpu(vol4h, sn.c[0], in, iy, 0);
    _suNf_read_spinor_gpu(vol4h, sn.c[1], in, iy, 2);
    _vector_i_add_assign_f(sn.c[0], sn.c[1]);

    _suNf_read_gpu(vol4h, u, gauge, ix+16*vol4h, 3);
    _suNf_theta_Z_multiply(sn.c[1], u, sn.c[0]);

    _vector_add_assign_f(r.c[0], sn.c[1]);
    _vector_i_sub_assign_f(r.c[2], sn.c[1]);

    _suNf_read_spinor_gpu(vol4h, sn.c[0], in, iy, 1);
    _suNf_read_spinor_gpu(vol4h, sn.c[1], in, iy, 3);
    _vector_i_sub_assign_f(sn.c[0], sn.c[1]);

    _suNf_theta_Z_multiply(sn.c[1], u, sn.c[0]);

    _vector_add_assign_f(r.c[1], sn.c[1]);
    _vector_i_add_assign_f(r.c[3], sn.c[1]);

    __syncthreads();
    /******************************* direction -3 *********************************/
    iy=idn_d[4*(ix+vol4h)+3];

    _suNf_read_spinor_gpu(vol4h, sn.c[0], in, iy, 0);
    _suNf_read_spinor_gpu(vol4h, sn.c[1], in, iy, 2);
    _vector_i_sub_assign_f(sn.c[0], sn.c[1]);

    _suNf_read_gpu(vol4h, u, gauge, iy, 3);
    _suNf_theta_Z_inverse_multiply(sn.c[1], u, sn.c[0]);

    _vector_add_assign_f(r.c[0], sn.c[1]);
    _vector_i_add_assign_f(r.c[2], sn.c[1]);

    _suNf_read_spinor_gpu(vol4h, sn.c[0], in, iy, 1);
    _suNf_read_spinor_gpu(vol4h, sn.c[1], in, iy, 3);
    _vector_i_add_assign_f(sn.c[0], sn.c[1]);

    _suNf_theta_Z_inverse_multiply(sn.c[1], u, sn.c[0]);

    _vector_add_assign_f(r.c[1], sn.c[1]);
    _vector_i_sub_assign_f(r.c[3], sn.c[1]);

    __syncthreads();
    /******************************** end of directions *********************************/
    _spinor_mul_f(r, -0.5, r);

    _suNf_write_spinor_gpu(vol4h, r.c[0], out, ix, 0);
    _suNf_write_spinor_gpu(vol4h, r.c[1], out, ix, 1);
    _suNf_write_spinor_gpu(vol4h, r.c[2], out, ix, 2);
    _suNf_write_spinor_gpu(vol4h, r.c[3], out, ix, 3);
  }
}



__global__ void Dphi_gpu_eo(suNf_spinor* __restrict__ out, const suNf_spinor* __restrict__ in,
                            const suNf* __restrict__ gauge, const int* __restrict__ iup_d, 
                            const int* __restrict__ idn_d, const int vol4h)
{
  suNf_spinor r;
  suNf_hspinor sn;
  suNf u;
#ifdef FERMION_THETA
  suNf_vector vtmp;
#endif

  int iy;
  int ix = blockIdx.x*BLOCK_SIZE + threadIdx.x;

  if (ix < vol4h) { 

    /******************************* direction +0 *********************************/
    iy=iup_d[4*ix];

    _suNf_read_spinor_gpu(vol4h, sn.c[0], in, iy-vol4h, 0);
    _suNf_read_spinor_gpu(vol4h, sn.c[1], in, iy-vol4h, 2);
    _suNf_read_gpu(vol4h, u, gauge, ix, 0);

    _vector_add_assign_f(sn.c[0], sn.c[1]);
    _suNf_theta_T_multiply(r.c[0], u, sn.c[0]);

    r.c[2]=r.c[0];

    _suNf_read_spinor_gpu(vol4h, sn.c[0], in, iy-vol4h, 1);
    _suNf_read_spinor_gpu(vol4h, sn.c[1], in, iy-vol4h, 3);
    _vector_add_assign_f(sn.c[0], sn.c[1]);

    _suNf_theta_T_multiply(r.c[1], u, sn.c[0]);

    r.c[3]=r.c[1];

    __syncthreads();
    /******************************* direction -0 *********************************/
    iy=idn_d[4*ix];

    _suNf_read_spinor_gpu(vol4h, sn.c[0], in, iy-vol4h, 0);
    _suNf_read_spinor_gpu(vol4h, sn.c[1], in, iy-vol4h, 2);
    _suNf_read_gpu(vol4h, u, gauge, iy+15*vol4h, 0);

    _vector_sub_assign_f(sn.c[0], sn.c[1]);

    _suNf_theta_T_inverse_multiply(sn.c[1], u, sn.c[0]);

    _vector_add_assign_f(r.c[0], sn.c[1]);
    _vector_sub_assign_f(r.c[2], sn.c[1]);
    
    _suNf_read_spinor_gpu(vol4h, sn.c[0], in, iy-vol4h, 1);
    _suNf_read_spinor_gpu(vol4h, sn.c[1], in, iy-vol4h, 3);
    _vector_sub_assign_f(sn.c[0], sn.c[1]);

    _suNf_theta_T_inverse_multiply(sn.c[1], u, sn.c[0]);
    
    _vector_add_assign_f(r.c[1], sn.c[1]);
    _vector_sub_assign_f(r.c[3], sn.c[1]);

    __syncthreads();
    /******************************* direction +1 *********************************/
    iy=iup_d[4*ix+1];

    _suNf_read_spinor_gpu(vol4h, sn.c[0], in, iy-vol4h, 0);
    _suNf_read_spinor_gpu(vol4h, sn.c[1], in, iy-vol4h, 3);
    _suNf_read_gpu(vol4h, u, gauge, ix, 1);

    _vector_i_add_assign_f(sn.c[0], sn.c[1]);
    _suNf_theta_X_multiply(sn.c[1], u, sn.c[0]);

    _vector_add_assign_f(r.c[0], sn.c[1]);
    _vector_i_sub_assign_f(r.c[3], sn.c[1]);

    _suNf_read_spinor_gpu(vol4h, sn.c[0], in, iy-vol4h, 1);
    _suNf_read_spinor_gpu(vol4h, sn.c[1], in, iy-vol4h, 2);
    _vector_i_add_assign_f(sn.c[0], sn.c[1]);

    _suNf_theta_X_multiply(sn.c[1], u, sn.c[0]);

    _vector_add_assign_f(r.c[1], sn.c[1]);
    _vector_i_sub_assign_f(r.c[2], sn.c[1]);


    __syncthreads();
    /******************************* direction -1 *********************************/
    iy=idn_d[4*ix+1];

    _suNf_read_spinor_gpu(vol4h, sn.c[0], in, iy-vol4h, 0);
    _suNf_read_spinor_gpu(vol4h, sn.c[1], in, iy-vol4h, 3);
    _suNf_read_gpu(vol4h, u, gauge, iy+15*vol4h, 1);

    _vector_i_sub_assign_f(sn.c[0], sn.c[1]);
    _suNf_theta_X_inverse_multiply(sn.c[1], u, sn.c[0]);

    _vector_add_assign_f(r.c[0], sn.c[1]);
    _vector_i_add_assign_f(r.c[3], sn.c[1]);

    _suNf_read_spinor_gpu(vol4h, sn.c[0], in, iy-vol4h, 1);
    _suNf_read_spinor_gpu(vol4h, sn.c[1], in, iy-vol4h, 2);
    _vector_i_sub_assign_f(sn.c[0], sn.c[1]);

    _suNf_theta_X_inverse_multiply(sn.c[1], u, sn.c[0]);

    _vector_add_assign_f(r.c[1], sn.c[1]);
    _vector_i_add_assign_f(r.c[2], sn.c[1]);


    __syncthreads();
    /******************************* direction +2 *********************************/
    iy=iup_d[4*ix+2];

    _suNf_read_spinor_gpu(vol4h, sn.c[0], in, iy-vol4h, 0);
    _suNf_read_spinor_gpu(vol4h, sn.c[1], in, iy-vol4h, 3);
    _vector_add_assign_f(sn.c[0], sn.c[1]);

    _suNf_read_gpu(vol4h, u, gauge, ix, 2);
    _suNf_theta_Y_multiply(sn.c[1], u, sn.c[0]);

    _vector_add_assign_f(r.c[0], sn.c[1]);
    _vector_add_assign_f(r.c[3], sn.c[1]);

    _suNf_read_spinor_gpu(vol4h, sn.c[0], in, iy-vol4h, 1);
    _suNf_read_spinor_gpu(vol4h, sn.c[1], in, iy-vol4h, 2);
    _vector_sub_assign_f(sn.c[0], sn.c[1]);

    _suNf_theta_Y_multiply(sn.c[1], u, sn.c[0]);

    _vector_add_assign_f(r.c[1], sn.c[1]);
    _vector_sub_assign_f(r.c[2], sn.c[1]);

    __syncthreads();
    /******************************* direction -2 *********************************/
    iy=idn_d[4*ix+2];

    _suNf_read_spinor_gpu(vol4h, sn.c[0], in, iy-vol4h, 0);
    _suNf_read_spinor_gpu(vol4h, sn.c[1], in, iy-vol4h, 3);
    _vector_sub_assign_f(sn.c[0], sn.c[1]);

    _suNf_read_gpu(vol4h, u, gauge, iy+15*vol4h, 2);
    _suNf_theta_Y_inverse_multiply(sn.c[1], u, sn.c[0]);

    _vector_add_assign_f(r.c[0], sn.c[1]);
    _vector_sub_assign_f(r.c[3], sn.c[1]);

    _suNf_read_spinor_gpu(vol4h, sn.c[0], in, iy-vol4h, 1);
    _suNf_read_spinor_gpu(vol4h, sn.c[1], in, iy-vol4h, 2);
    _vector_add_assign_f(sn.c[0], sn.c[1]);

    _suNf_theta_Y_inverse_multiply(sn.c[1], u, sn.c[0]);

    _vector_add_assign_f(r.c[1], sn.c[1]);
    _vector_add_assign_f(r.c[2], sn.c[1]);

    __syncthreads();
    /******************************* direction +3 *********************************/
    iy=iup_d[4*ix+3];

    _suNf_read_spinor_gpu(vol4h, sn.c[0], in, iy-vol4h, 0);
    _suNf_read_spinor_gpu(vol4h, sn.c[1], in, iy-vol4h, 2);
    _vector_i_add_assign_f(sn.c[0], sn.c[1]);

    _suNf_read_gpu(vol4h, u, gauge, ix, 3);
    _suNf_theta_Z_multiply(sn.c[1], u, sn.c[0]);

    _vector_add_assign_f(r.c[0], sn.c[1]);
    _vector_i_sub_assign_f(r.c[2], sn.c[1]);

    _suNf_read_spinor_gpu(vol4h, sn.c[0], in, iy-vol4h, 1);
    _suNf_read_spinor_gpu(vol4h, sn.c[1], in, iy-vol4h, 3);
    _vector_i_sub_assign_f(sn.c[0], sn.c[1]);

    _suNf_theta_Z_multiply(sn.c[1], u, sn.c[0]);

    _vector_add_assign_f(r.c[1], sn.c[1]);
    _vector_i_add_assign_f(r.c[3], sn.c[1]);

    __syncthreads();
    /******************************* direction -3 *********************************/
    iy=idn_d[4*ix+3];

    _suNf_read_spinor_gpu(vol4h, sn.c[0], in, iy-vol4h, 0);
    _suNf_read_spinor_gpu(vol4h, sn.c[1], in, iy-vol4h, 2);
    _vector_i_sub_assign_f(sn.c[0], sn.c[1]);

    _suNf_read_gpu(vol4h, u, gauge, iy+15*vol4h, 3);
    _suNf_theta_Z_inverse_multiply(sn.c[1], u, sn.c[0]);

    _vector_add_assign_f(r.c[0], sn.c[1]);
    _vector_i_add_assign_f(r.c[2], sn.c[1]);

    _suNf_read_spinor_gpu(vol4h, sn.c[0], in, iy-vol4h, 1);
    _suNf_read_spinor_gpu(vol4h, sn.c[1], in, iy-vol4h, 3);
    _vector_i_add_assign_f(sn.c[0], sn.c[1]);

    _suNf_theta_Z_inverse_multiply(sn.c[1], u, sn.c[0]);

    _vector_add_assign_f(r.c[1], sn.c[1]);
    _vector_i_sub_assign_f(r.c[3], sn.c[1]);

    __syncthreads();
    /******************************** end of directions *********************************/
    _spinor_mul_f(r, -0.5, r);

    _suNf_write_spinor_gpu(vol4h, r.c[0], out, ix, 0);
    _suNf_write_spinor_gpu(vol4h, r.c[1], out, ix, 1);
    _suNf_write_spinor_gpu(vol4h, r.c[2], out, ix, 2);
    _suNf_write_spinor_gpu(vol4h, r.c[3], out, ix, 3);
  }
}

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


void Dphi_(spinor_field *out, spinor_field *in)
{
  int N, grid;
  const int vol4h=T*X*Y*Z/2;

  init_bc_gpu();

  error((in==NULL)||(out==NULL), 1, "Dphi_ [Dphi_gpu.c]",
         "Attempt to access unallocated memory space");

  error(in==out, 1, "Dphi_ [Dphi_gpu.c]",
         "Input and output fields must be different");

#ifndef CHECK_SPINOR_MATCHING
   error(out->type==&glat_even && in->type!=&glat_odd, 1, "Dphi_ [Dphi_gpu.c]", "Spinors don't match! (1)");
   error(out->type==&glat_odd && in->type!=&glat_even, 1, "Dphi_ [Dphi_gpu.c]", "Spinors don't match! (2)");
   error(out->type==&glattice && in->type!=&glattice, 1, "Dphi_ [Dphi_gpu.c]", "Spinors don't match! (3)");
#endif

  N = vol4h;
  grid = (N-1)/BLOCK_SIZE + 1;
  grid *= THREADSITE; //1 for KEPLER; 4 for FERMI

  if(in->type==&glat_odd) {
    std::cout << "Dphi_ glat_odd" << std::endl;
    ++MVMcounter;
    std::cout << "Before Dphi_gpu_eo:" << sqrt(spinor_field_sqnorm_f(in)) << std::endl;
    std::cout << "Before Dphi_gpu_eo sanity:" << sqrt(spinor_field_sqnorm_f(in)) << std::endl;
    Dphi_gpu_eo<<<grid, BLOCK_SIZE>>>(START_SP_ADDRESS_GPU(out), START_SP_ADDRESS_GPU(in), u_gauge_f->gpu_ptr, iup_gpu, idn_gpu, vol4h);
    CudaCheckError();
    std::cout << "After Dphi_gpu_eo:" << sqrt(spinor_field_sqnorm_f(in)) << std::endl;
  } else if (in->type==&glat_even) {
    std::cout << "Dphi_ glat_even" << std::endl;
    std::cout << "Before Dphi_gpu_oe:" << sqrt(spinor_field_sqnorm_f(in)) << std::endl;
    std::cout << "Before Dphi_gpu_oe sanity:" << sqrt(spinor_field_sqnorm_f(in)) << std::endl;
    ++MVMcounter;
    Dphi_gpu_oe<<<grid, BLOCK_SIZE>>>(START_SP_ADDRESS_GPU(out), START_SP_ADDRESS_GPU(in), u_gauge_f->gpu_ptr, iup_gpu, idn_gpu, vol4h);
    CudaCheckError();
    std::cout << "After Dphi_gpu_oe:" << sqrt(spinor_field_sqnorm_f(in)) << std::endl;
  } else if (in->type==&glattice) {
    std::cout << "Dphi_ glattice" << std::endl;
    in->type=&glat_even;
    out->type=&glat_odd;
    ++MVMcounter;
    Dphi_gpu_oe<<<grid, BLOCK_SIZE>>>(START_SP_ADDRESS_GPU(out), START_SP_ADDRESS_GPU(in), u_gauge_f->gpu_ptr, iup_gpu, idn_gpu, vol4h);
    CudaCheckError();
    in->type=&glat_odd;
    out->type=&glat_even;
    ++MVMcounter;
    Dphi_gpu_eo<<<grid, BLOCK_SIZE>>>(START_SP_ADDRESS_GPU(out), START_SP_ADDRESS_GPU(in), u_gauge_f->gpu_ptr, iup_gpu, idn_gpu, vol4h);
    CudaCheckError();

    in->type=&glattice;
    out->type=&glattice;
  } else {
    error(1, 1, "Dphi_ [Dphi_gpu.c]", "Wrong input spinor geometry!");
  }
   CudaCheckError();
}


/*
 * this function takes 2 spinors defined on the whole lattice
 */
void Dphi(double m0, spinor_field *out, spinor_field *in)
{
  double rho;

  error((in==NULL)||(out==NULL), 1, "Dphi [Dphi_gpu.c]",
        "Attempt to access unallocated memory space");

  error(in==out, 1, "Dphi [Dphi_gpu.c]",
        "Input and output fields must be different");

#ifdef CHECK_SPINOR_MATCHING
  error(out->type!=&glattice || in->type!=&glattice, 1, "Dphi [Dphi_gpu.c]", "Spinors are not defined on all the lattice!");
#endif /* CHECK_SPINOR_MATCHING */

  Dphi_(out, in);

  rho = 4. + m0;
  spinor_field_mul_add_assign_f(out, rho, in);
}


void g5Dphi(double m0, spinor_field *out, spinor_field *in)
{
  double rho;

  error((in==NULL)||(out==NULL), 1, "g5Dphi [Dphi_gpu.c]",
	"Attempt to access unallocated memory space");

  error(in==out, 1, "g5Dphi [Dphi_gpu.c]",
	"Input and output fields must be different");

#ifdef CHECK_SPINOR_MATCHING
   error(out->type!=&glattice || in->type!=&glattice, 1, "g5Dphi [Dphi_gpu.c]", "Spinors are not defined on all the lattice!");
#endif /* CHECK_SPINOR_MATCHING */

  Dphi_(out, in);
  rho=4.+m0;
  spinor_field_mul_add_assign_f(out, rho, in);
  spinor_field_g5_assign_f(out);
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
void Dphi_eopre(double m0, spinor_field *out, spinor_field *in)
{
  double rho;

  error((in==NULL)||(out==NULL), 1, "Dphi_eopre [Dphi_gpu.c]",
	"Attempt to access unallocated memory space");

  error(in==out, 1, "Dphi_eopre [Dphi_gpu.c]",
	"Input and output fields must be different");

#ifdef CHECK_SPINOR_MATCHING
  error(out->type!=&glat_even || in->type!=&glat_even, 1, "Dphi_eopre " __FILE__, "Spinors are not defined on even lattice!");
#endif /* CHECK_SPINOR_MATCHING */

  /* alloc memory for temporary spinor field */
  if (init) { init_Dirac(); }

  Dphi_(otmp, in);
  Dphi_(out, otmp);

  rho=4.0+m0;
  rho*=-rho; /* this minus sign is taken into account below */

  spinor_field_mul_add_assign_f(out, rho, in);
  spinor_field_minus_f(out, out);
}


/* Even/Odd preconditioned dirac operator
 * this function takes 2 spinors defined on the odd lattice
 * Dphi in = (4+m0)^2*in - D_OE D_EO in
 *
 */
void Dphi_oepre(double m0, spinor_field *out, spinor_field *in)
{
  double rho;

  error((in==NULL)||(out==NULL), 1, "Dphi_oepre [Dphi_gpu.c]",
	"Attempt to access unallocated memory space");

  error(in==out, 1, "Dphi_oepre [Dphi_gpu.c]",
	"Input and output fields must be different");

#ifdef CHECK_SPINOR_MATCHING
  error(out->type!=&glat_odd || in->type!=&glat_odd, 1, "Dphi_oepre " __FILE__, "Spinors are not defined on odd lattice!");
#endif /* CHECK_SPINOR_MATCHING */


  /* alloc memory for temporary spinor field */
  if (init) { init_Dirac();}

  Dphi_(etmp, in);
  Dphi_(out, etmp);

  rho=4.0+m0;
  rho*=-rho; /* this minus sign is taken into account below */

  spinor_field_mul_add_assign_f(out, rho, in);
  spinor_field_minus_f(out, out);
}


void g5Dphi_eopre(double m0, spinor_field *out, spinor_field *in)
{
  double rho;
  std::cout << "In g5Dphi_eopre()" << std::endl;

  error((in==NULL)||(out==NULL), 1, "g5Dphi_eopre [Dphi_gp.c]",
	"Attempt to access unallocated memory space");

  error(in==out, 1, "Dphi_eopre [Dphi_gpu.c]",
	"Input and output fields must be different");

#ifdef CHECK_SPINOR_MATCHING
  error(out->type!=&glat_even || in->type!=&glat_even, 1, "g5Dphi_eopre " __FILE__, "Spinors are not defined on even lattice!");
#endif /* CHECK_SPINOR_MATCHING */

#if defined(BASIC_SF) || defined(ROTATED_SF)
  SF_spinor_bcs(in);
#endif /* defined(BASIC_SF) || defined(ROTATED_SF) */

  /* alloc memory for temporary spinor field */
  if (init) { init_Dirac();}

  std::cout << "Before Dphi_:" << sqrt(spinor_field_sqnorm_f(in)) << std::endl;
  std::cout << "   Before in ptr:" << (in)->gpu_ptr << std::endl;
  std::cout << "   Before out ptr:" << (out)->gpu_ptr << std::endl;
  std::cout << "   Before otmp ptr:" << (otmp)->gpu_ptr << std::endl;
  Dphi_(otmp, in);
  std::cout << "Between Dphi_:" << sqrt(spinor_field_sqnorm_f(in)) << std::endl;
  std::cout << "   Between in ptr:" << (in)->gpu_ptr << std::endl;
  std::cout << "   Between out ptr:" << (out)->gpu_ptr << std::endl;
  std::cout << "   Between otmp ptr:" << (otmp)->gpu_ptr << std::endl;
  Dphi_(out, otmp);
  std::cout << "After Dphi_:" << sqrt(spinor_field_sqnorm_f(in)) << std::endl;
  std::cout << "   After in ptr:" << (in)->gpu_ptr << std::endl;
  std::cout << "   After out ptr:" << (out)->gpu_ptr << std::endl;
  std::cout << "   After otmp ptr:" << (otmp)->gpu_ptr << std::endl;

  rho=4.0+m0;
  rho*=-rho; /* this minus sign is taken into account below */

  spinor_field_mul_add_assign_f(out, rho, in);
  spinor_field_minus_f(out, out);
  spinor_field_g5_assign_f(out);

#if defined(BASIC_SF) || defined(ROTATED_SF)
  SF_spinor_bcs(out);
#endif /* defined(BASIC_SF) || defined(ROTATED_SF) */
}


/* g5Dphi_eopre ^2 */
void g5Dphi_eopre_sq(double m0, spinor_field *out, spinor_field *in) {
  /* alloc memory for temporary spinor field */
  std::cout << "In g5Dphi_eopre_sq" << std::endl;
  if (init) { init_Dirac(); }

  g5Dphi_eopre(m0, etmp, in);
  g5Dphi_eopre(m0, out, etmp);
}


/* g5Dhi ^2 */
void g5Dphi_sq(double m0, spinor_field *out, spinor_field *in) {
  /* alloc memory for temporary spinor field */
  if (init) { init_Dirac();  }

  g5Dphi(m0, gtmp, in);
  g5Dphi(m0, out, gtmp);
}

#endif
