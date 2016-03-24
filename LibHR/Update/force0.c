/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *   
* All rights reserved.                                                      * 
\***************************************************************************/

#include "global.h"
#include "update.h"
#include "suN.h"
#include "utils.h"
#include "representation.h"
#include "logger.h"
#include "communications.h"

#include <stdio.h>
#include <math.h>


#define _print_avect(a) printf("(%3.5f,%3.5f,%3.5f,%3.5f,%3.5f,%3.5f,%3.5f,%3.5f)\n",(a).c1,(a).c2,(a).c3,(a).c4,(a).c5,(a).c6,(a).c7,(a).c8)

#define _print_mat(a) printf("(%3.5f,%3.5f,%3.5f)\n(%3.5f,%3.5f,%3.5f)\n(%3.5f,%3.5f,%3.5f)\n",(a).c1_1.re,(a).c1_2.re,(a).c1_3.re,(a).c2_1.re,(a).c2_2.re,(a).c2_3.re,(a).c3_1.re,(a).c3_2.re,(a).c3_3.re);printf("(%3.5f,%3.5f,%3.5f)\n(%3.5f,%3.5f,%3.5f)\n(%3.5f,%3.5f,%3.5f)\n",(a).c1_1.im,(a).c1_2.im,(a).c1_3.im,(a).c2_1.im,(a).c2_2.im,(a).c2_3.im,(a).c3_1.im,(a).c3_2.im,(a).c3_3.im)


void force0(double dt, suNg_av_field *force, void *vpar){
  #ifdef TIMING
  struct timeval start, end;
  struct timeval etime;
  
  #ifdef TIMING_WITH_BARRIERS
  MPI_Barrier(GLB_COMM);
  #endif
  gettimeofday(&start,0);
  #endif

  /* check input types */
  _TWO_SPINORS_MATCHING(u_gauge,force);
  
  double beta = *((double*)vpar);
#ifdef MEASURE_FORCE0
  double forcestat[2]={0.,0.}; /* used for computation of avr and max force */
  _MASTER_FOR_SUM(&glattice,i,forcestat[0],forcestat[1]) {
#else
  _MASTER_FOR(&glattice,i) {
#endif
    suNg s1,s2;
    suNg_algebra_vector f;
    for (int mu=0; mu<4; ++mu) {
#ifdef TLSYM
      
      double c0,c1;
      c1 = -1./12.;
      c0 = 1. - 8.*c1;
      static suNg s3,s4;
      /* pure Wilson part */
      staples(i,mu,&s1);
      _suNg_times_suNg_dagger(s2,*_4FIELD_AT(u_gauge,i,mu),s1);
      
      /* the projection itself takes the TA: proj(M) = proj(TA(M)) */
      _fund_algebra_project(f,s2);
      _algebra_vector_mul_add_assign_g(*_4FIELD_AT(force,i,mu), c0*dt*(-beta/((double)(NG))), f);
      
      /* rectangular part */
      rect_staples_1x2(i,mu,&s3);
      _suNg_times_suNg_dagger(s4,*_4FIELD_AT(u_gauge,i,mu),s3);
      
      /* the projection itself takes the TA: proj(M) = proj(TA(M)) */
      _fund_algebra_project(f,s4);
    
      _algebra_vector_mul_add_assign_g(*_4FIELD_AT(force,i,mu), c1*dt*(-beta/((double)(NG))), f);
#else
      staples(i,mu,&s1);
      _suNg_times_suNg_dagger(s2,*_4FIELD_AT(u_gauge,i,mu),s1);

      /* the projection itself takes the TA: proj(M) = proj(TA(M)) */
      _fund_algebra_project(f,s2);
      _algebra_vector_mul_add_assign_g(*_4FIELD_AT(force,i,mu), dt*(-beta/((double)(NG))), f);
#endif 

#ifdef MEASURE_FORCE0
      double nsq;
      _algebra_vector_sqnorm_g(nsq,f);
      forcestat[0]+=sqrt(nsq);
      for(int x=0;x<sizeof(suNg_algebra_vector)/sizeof(double);++x){
        if(forcestat[1]<fabs(f.c[x])) forcestat[1]=fabs(f.c[x]);
      }
#endif
    }
  }

#ifdef MEASURE_FORCE0
  //  if(logger_getlevel("FORCE-STAT")>=10){
    global_sum(forcestat,1);
    global_max(forcestat+1,1);
    
    forcestat[0]*=beta/((4.*NG)*GLB_VOLUME);
    forcestat[1]*=beta/((double)(NG));
    //    lprintf("FORCE-STAT",10," force0 : dt= %1.8e avr |force|= %1.8e maxforce= %1.8e \n",dt,forcestat[0],forcestat[1]);
    lprintf("FORCE_STAT",20,"GF: avr dt |force| = %1.8e dt maxforce = %1.8e, dt = %1.8e \n",forcestat[0]*dt,forcestat[1]*dt,dt);
    force_ave[0]+=dt*forcestat[0];
    force_max[0]+=dt*forcestat[1];    
    //  }
#endif
  
  apply_BCs_on_momentum_field(force);

  #ifdef TIMING
  #ifdef TIMING_WITH_BARRIERS
  MPI_Barrier(GLB_COMM);
  #endif
  gettimeofday(&end,0);
  timeval_subtract(&etime,&end,&start);
  lprintf("TIMING",0,"force0 %.6f s\n",1.*etime.tv_sec+1.e-6*etime.tv_usec);
  #endif
}

