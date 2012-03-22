/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *   
* All rights reserved.                                                      * 
\***************************************************************************/
#ifdef WITH_GPU

#include "global.h"
#include "update.h"
#include "suN.h"
#include "utils.h"
#include "representation.h"
#include "logger.h"
#include "communications.h"

#include <stdio.h>
#include <math.h>

extern rhmc_par _update_par;

#define _suNg_read_gpu(stride,v,in,iy,x)\
iw=(iy)+((x)*4)*(stride);\
(v).c[0]=((double*)(in))[iw]; iw+=(stride); \
(v).c[1]=((double*)(in))[iw]; iw+=(stride);\
(v).c[2]=((double*)(in))[iw]; iw+=(stride);\
(v).c[3]=((double*)(in))[iw]

#define _suNg_av_read_gpu(stride,v,in,iy,x)\
iw=(iy)+((x)*3)*(stride);\
(v).c[0]=((double*)(in))[iw]; iw+=(stride); \
(v).c[1]=((double*)(in))[iw]; iw+=(stride);\
(v).c[2]=((double*)(in))[iw]

#define _suNg_av_write_gpu(stride,v,in,iy,x)\
iw=(iy)+((x)*3)*(stride);\
(v).c[0]=((double*)(in))[iw]; iw+=(stride); \
(v).c[1]=((double*)(in))[iw]; iw+=(stride);\
(v).c[2]=((double*)(in))[iw]


void force0(double dt, suNg_av_field *force, void *vpar){
  static suNg s1,s2;
  static suNg_algebra_vector f;
  double forcestat[2]={0.,0.}; /* used for computation of avr and max force */
  double nsq;
  int mu,x;

  /* check input types */
  _TWO_SPINORS_MATCHING(u_gauge,force);

  int N = T*X*Y*Z;//u_gauge->type->master_end[0] -  u_gauge->type->master_start[0] + 1;
  int grid = N/BLOCK_SIZE + ((N % BLOCK_SIZE == 0) ? 0 : 1);

  

  _MASTER_FOR(&glattice,i) {
    for (mu=0; mu<4; ++mu) {
      staples(i,mu,&s1);
      _suNg_times_suNg_dagger(s2,*_4FIELD_AT(u_gauge,i,mu),s1);
    
      /* the projection itself takes the TA: proj(M) = proj(TA(M)) */
      _fund_algebra_project(f,s2);
    
      _algebra_vector_mul_add_assign_g(*_4FIELD_AT(force,i,mu), dt*(-_update_par.beta/((double)(NG))), f);

      _algebra_vector_sqnorm_g(nsq,f);
      forcestat[0]+=sqrt(nsq);
      for(x=0;x<NG*NG-1;++x){
	if(forcestat[1]<fabs(*(((double*)&f)+x))) forcestat[1]=fabs(*(((double*)&f)+x));
      }
    }
  }
	
  global_sum(forcestat,2);
  forcestat[0]*=dt*_update_par.beta/((double)(NG*4*GLB_T*GLB_X*GLB_Y*GLB_Z));
  forcestat[1]*=dt*_update_par.beta/((double)NG);
  lprintf("FORCE0",50,"avr |force| = %1.8e maxforce = %1.8e\n",forcestat[0],forcestat[1]);
  
#if defined(BASIC_SF) || defined(ROTATED_SF)
	SF_force_bcs(force);
#endif /* BASIC_SF || ROTATED_SF */

  }

#undef _suNg_read_gpu
#undef _suNg_av_read_gpu
#undef _suNg_av_write_gpu

#endif

/*
void Force(double dt, suNg_av_field *force, spinor_field *pf){
  Force0(dt, force);
  Force_rhmc_f(dt, force, pf);
}
*/

