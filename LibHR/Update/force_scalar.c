/***************************************************************************\
* Copyright (c) 2017
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


static void force_scalar_s(double dt, void *vpar){
  #ifdef TIMING
  struct timeval start, end;
  struct timeval etime;
  
  #ifdef TIMING_WITH_BARRIERS
  MPI_Barrier(GLB_COMM);
  #endif
  gettimeofday(&start,0);
  #endif

  force_scalar_par *par = (force_scalar_par*)vpar;
  /* check input types */
  _TWO_SPINORS_MATCHING(u_scalar,par->momenta);
  
#ifdef MEASURE_FORCE0
  double forcestat[2]={0.,0.}; /* used for computation of avr and max force */
  _MASTER_FOR_SUM(&glattice,i,forcestat[0],forcestat[1]) {
#else
  _MASTER_FOR(&glattice,i) {
#endif
    suNg_vector f;
    _vector_zero_g(f);
    for (int mu=0; mu<4; ++mu) {
	suNg U_xmu = *pu_gauge(i, mu); // U_mu(x)
	suNg U_down = *pu_gauge(idn(i,mu),mu); // U_mu(x-mu)
	suNg_vector Splus = *pu_scalar(iup(i,mu)); // S(x+mu)
	suNg_vector Sminus = *pu_scalar(idn(i,mu)); // S(x-mu)
	suNg_vector SSum, Down, Up;
	_suNg_inverse_multiply(Down,U_down,Sminus);//U_mu(x-mu)^+S(x-mu) 
	_suNg_multiply(Up,U_xmu,Splus);//U_mu(x)S(x+mu)
	_vector_add_g(SSum, Up, Down);
	suNg_vector USStar;
	vector_star(&USStar, &SSum);
	_vector_sub_assign_g(f, USStar);
    }
    double Msq = par->mass;
    Msq = Msq*Msq + 8.0;
    double lambda = par->lambda;
	
    suNg_vector S = *pu_scalar(i);
    suNg_vector SStar;
    vector_star(&SStar, &S);
    double S2;
    _vector_prod_re_g(S2,S,S);
    _vector_lc_add_assign_g(f,Msq,SStar,2.0*lambda*S2,SStar);
    
    _vector_mul_add_assign_g(*_FIELD_AT(*par->momenta,i),-dt,f);
  }

  #ifdef TIMING
  #ifdef TIMING_WITH_BARRIERS
  MPI_Barrier(GLB_COMM);
  #endif
  gettimeofday(&end,0);
  timeval_subtract(&etime,&end,&start);
  lprintf("TIMING",0,"force0 %.6f s\n",1.*etime.tv_sec+1.e-6*etime.tv_usec);
  #endif
}

static void outer_product(suNg* u, suNg_vector* v1, suNg_vector* v2){
	for(int i=0;i<NG*NG;i++){
		int row,col;
		row = i/NG;
		col = i%NG;
		_complex_mul_star(u->c[i],v1->c[row],v2->c[col]);
	}
} 

static void force_scalar_g(double dt, void *vpar){
  #ifdef TIMING
  struct timeval start, end;
  struct timeval etime;
  
  #ifdef TIMING_WITH_BARRIERS
  MPI_Barrier(GLB_COMM);
  #endif
  gettimeofday(&start,0);
  #endif

  force_scalar_par *par = (force_scalar_par*)vpar;
  /* check input types */
  _TWO_SPINORS_MATCHING(u_gauge,par->g_momenta);
  
#ifdef MEASURE_FORCE0
  double forcestat[2]={0.,0.}; /* used for computation of avr and max force */
  _MASTER_FOR_SUM(&glattice,i,forcestat[0],forcestat[1]) {
#else
  _MASTER_FOR(&glattice,i) {
#endif
    suNg s1,s2;
    suNg_algebra_vector f;
    for (int mu=0; mu<4; ++mu) {
//	    outer_product(&s1,pu_scalar(i),pu_scalar(iup(i,mu)) );
	    outer_product(&s1,pu_scalar(iup(i,mu)),pu_scalar(i));
	    _suNg_times_suNg(s2,*_4FIELD_AT(u_gauge,i,mu),s1);

	    /* the projection itself takes the TA: proj(M) = proj(TA(M)) */
	    _fund_algebra_project(f,s2);
	    _algebra_vector_mul_add_assign_g(*_4FIELD_AT(*par->g_momenta,i,mu), -dt*2.0, f);

    }
  }

  apply_BCs_on_momentum_field(*par->g_momenta);

  #ifdef TIMING
  #ifdef TIMING_WITH_BARRIERS
  MPI_Barrier(GLB_COMM);
  #endif
  gettimeofday(&end,0);
  timeval_subtract(&etime,&end,&start);
  lprintf("TIMING",0,"force0 %.6f s\n",1.*etime.tv_sec+1.e-6*etime.tv_usec);
  #endif
}

void force_scalar(double dt, void* par){
	force_scalar_s(dt, par);
	force_scalar_g(dt, par);
}
