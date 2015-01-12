/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *   
* All rights reserved.                                                      * 
\***************************************************************************/

#include "global.h"
#include "update.h"
#include "suN.h"
#include "utils.h"
#include "dirac.h"
#include "inverters.h"
#include "rational_functions.h"
#include "representation.h"
#include "logger.h"
#include "linear_algebra.h"
#include "memory.h"
#include "communications.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>


#define _print_avect(a) printf("(%3.5e,%3.5e,%3.5e,%3.5e,%3.5e,%3.5e,%3.5e,%3.5e)\n",(a).c1,(a).c2,(a).c3,(a).c4,(a).c5,(a).c6,(a).c7,(a).c8)

#define _print_mat(a) printf("(%3.5f,%3.5f,%3.5f)\n(%3.5f,%3.5f,%3.5f)\n(%3.5f,%3.5f,%3.5f)\n",(a).c1_1.re,(a).c1_2.re,(a).c1_3.re,(a).c2_1.re,(a).c2_2.re,(a).c2_3.re,(a).c3_1.re,(a).c3_2.re,(a).c3_3.re);printf("(%3.5f,%3.5f,%3.5f)\n(%3.5f,%3.5f,%3.5f)\n(%3.5f,%3.5f,%3.5f)\n",(a).c1_1.im,(a).c1_2.im,(a).c1_3.im,(a).c2_1.im,(a).c2_2.im,(a).c2_3.im,(a).c3_1.im,(a).c3_2.im,(a).c3_3.im)

/* we need to compute  Tr  U(x,mu) g_5*(1-g_mu) chi2 # chi1^+
 * where # indicates the tensor product and Tr is the trace on Lorentz space.
 * the strategy is the following:
 * given the form of g_5(1-g_mu) one can compute only the first two lorentz
 * components of the spinor; so we first apply g_5(1-g_mu) to chi2 to find the first
 * two components; then we multiply these two vectors by U(x,mu) and
 * store the result in p.c[0], p.c[1]; when computing the trace we can factorize p.c[0] and p.c[1]
 * as they both multiply two components of chi1^+; we store these factors in p.c[2] and p.c[3].
 * the tensor product is performed by the macro 
 * _suNf_FMAT(u,p): u = p.c[0] # p.c[2]^+ + p.c[1] # p.c[3]^+
 */

/* these macros use the variables ptmp, p */

#ifdef BC_T_THETA
#define _T_theta_mulc(r) _vector_mulc_f(ptmp,eitheta[0],(r)); (r)=ptmp
#else
#define _T_theta_mulc(r)
#endif
#ifdef BC_X_THETA
#define _X_theta_mulc(r) _vector_mulc_f(ptmp,eitheta[1],(r)); (r)=ptmp
#else
#define _X_theta_mulc(r)
#endif
#ifdef BC_Y_THETA
#define _Y_theta_mulc(r) _vector_mulc_f(ptmp,eitheta[2],(r)); (r)=ptmp
#else
#define _Y_theta_mulc(r)
#endif
#ifdef BC_Z_THETA
#define _Z_theta_mulc(r) _vector_mulc_f(ptmp,eitheta[3],(r)); (r)=ptmp
#else
#define _Z_theta_mulc(r)
#endif

#define _F_DIR0(u,chi1,chi2)				      \
  _vector_add_f(ptmp,(chi2)->c[0],(chi2)->c[2]);		      \
  _suNf_multiply(p.c[0],*(pu_gauge_f(x,0)),ptmp);		      \
  _T_theta_mulc(p.c[0]);                                      \
  _vector_add_f(ptmp,(chi2)->c[1],(chi2)->c[3]);		      \
  _suNf_multiply(p.c[1],*(pu_gauge_f(x,0)),ptmp);		      \
  _T_theta_mulc(p.c[1]);                                      \
  _vector_sub_f(p.c[2],(chi1)->c[0],(chi1)->c[2]);	      \
  _vector_sub_f(p.c[3],(chi1)->c[1],(chi1)->c[3]);	      \
  _suNf_FMAT((u),p)

#define _F_DIR1(u,chi1,chi2)				      \
  _vector_i_add_f(ptmp,(chi2)->c[0],(chi2)->c[3]);		      \
  _suNf_multiply(p.c[0],*(pu_gauge_f(x,1)),ptmp);		      \
  _X_theta_mulc(p.c[0]);                                      \
  _vector_i_add_f(ptmp,(chi2)->c[1],(chi2)->c[2]);		      \
  _suNf_multiply(p.c[1],*(pu_gauge_f(x,1)),ptmp);		      \
  _X_theta_mulc(p.c[1]);                                      \
  _vector_i_sub_f(p.c[2],(chi1)->c[0],(chi1)->c[3]);	      \
  _vector_i_sub_f(p.c[3],(chi1)->c[1],(chi1)->c[2]);	      \
  _suNf_FMAT((u),p)

#define _F_DIR2(u,chi1,chi2)				      \
  _vector_add_f(ptmp,(chi2)->c[0],(chi2)->c[3]);		      \
  _suNf_multiply(p.c[0],*(pu_gauge_f(x,2)),ptmp);		      \
  _Y_theta_mulc(p.c[0]);                                      \
  _vector_sub_f(ptmp,(chi2)->c[1],(chi2)->c[2]);		      \
  _suNf_multiply(p.c[1],*(pu_gauge_f(x,2)),ptmp);		      \
  _Y_theta_mulc(p.c[1]);                                      \
  _vector_sub_f(p.c[2],(chi1)->c[0],(chi1)->c[3]);	      \
  _vector_add_f(p.c[3],(chi1)->c[1],(chi1)->c[2]);	      \
  _suNf_FMAT((u),p)

#define _F_DIR3(u,chi1,chi2)				      \
  _vector_i_add_f(ptmp,(chi2)->c[0],(chi2)->c[2]);		      \
  _suNf_multiply(p.c[0],*(pu_gauge_f(x,3)),ptmp);		      \
  _Z_theta_mulc(p.c[0]);                                      \
  _vector_i_sub_f(ptmp,(chi2)->c[1],(chi2)->c[3]);		      \
  _suNf_multiply(p.c[1],*(pu_gauge_f(x,3)),ptmp);		      \
  _Z_theta_mulc(p.c[1]);                                      \
  _vector_i_sub_f(p.c[2],(chi1)->c[0],(chi1)->c[2]);	      \
  _vector_i_add_f(p.c[3],(chi1)->c[1],(chi1)->c[3]);	      \
  _suNf_FMAT((u),p)


  
static unsigned int n_ws=0;
static spinor_field *chi=NULL, *Hchi=NULL;

void init_force_rhmc() {
}

void free_force_rhmc() {
  free_spinor_field_f(chi);
}


  

void force_rhmc(double dt, suNg_av_field *force, void *vpar){
    
  #ifdef TIMING
  struct timeval start, end;
  struct timeval start1, end1;
  struct timeval etime;
  
  #ifdef TIMING_WITH_BARRIERS
  MPI_Barrier(GLB_COMM);
  #endif
  gettimeofday(&start,0);
  #endif

  static mshift_par inv_par;
  #ifdef UPDATE_EO
  spinor_field delta, sigma;
  #endif
  
  force_rhmc_par *par = (force_rhmc_par*)vpar;
  spinor_field *pf = par->pf;
  rational_app *ratio = par->ratio;
  
  if(n_ws<ratio->order+1) {
    if(chi!=NULL) free_spinor_field_f(chi);
    n_ws = ratio->order+1;
    chi = alloc_spinor_field_f(n_ws,&glattice);
    Hchi = chi+n_ws-1;
  }
  
  /* check input types */
  #ifndef CHECK_SPINOR_MATCHING
  _TWO_SPINORS_MATCHING(u_gauge,force);
  #endif
  
  /* Compute (H^2-b[n])^-1 * pf */
  /* set up cg parameters */
  inv_par.n = ratio->order;
  inv_par.shift = ratio->b;
  inv_par.err2= par->inv_err2; /* this should be high for reversibility */
  inv_par.max_iter=0; /* no limit */
  
  #ifdef UPDATE_EO
  /* change the type of chi[n] */
  for (int n=0; n<ratio->order; ++n) {
    chi[n].type=&glat_even;
  }
  Hchi->type=&glat_even; 
  #endif
  
  
  for (int k=0; k<par->n_pf; ++k) {
    /* compute inverse vectors chi[i] = (H^2 - b[i])^1 * pf */
    if (inv_par.n==1) { spinor_field_zero_f(chi); }

    #ifdef TIMING
    #ifdef TIMING_WITH_BARRIERS
    MPI_Barrier(GLB_COMM);
    #endif
    gettimeofday(&start1,0);
    #endif

    set_dirac_mass(par->mass);
    cg_mshift(&inv_par, &H2, &pf[k], chi);
    
    #ifdef TIMING
    #ifdef TIMING_WITH_BARRIERS
    MPI_Barrier(GLB_COMM);
    #endif
    gettimeofday(&end1,0);
    timeval_subtract(&etime,&end1,&start1);
    lprintf("TIMING",0,"cg_mshift in force_rhmc %.6f s\n",1.*etime.tv_sec+1.e-6*etime.tv_usec);
    #endif

    for (int n=0; n<ratio->order; ++n) {
      #ifdef UPDATE_EO
      /* change temporarely the type of chi[n] and Hchi */
      g5Dphi_eopre(par->mass, Hchi, &chi[n]);
      /* start_sf_sendrecv(Hchi); this is not needed since it is done by the following Dphi_ call*/
      /* copy the spinor field structures of chi[n] and Hchi */
      /* in this way we can reuse the odd part of the fields with new names */
      delta=*Hchi; 
      delta.ptr=Hchi->ptr+glat_odd.master_shift;
      delta.type=&glat_odd;
      sigma=chi[n]; 
      sigma.type=&glat_odd;
      sigma.ptr=chi[n].ptr+glat_odd.master_shift;
      Dphi_(&delta,Hchi);
      Dphi_(&sigma,&chi[n]);
      
      start_sf_sendrecv(&delta);
      start_sf_sendrecv(&sigma);
      #else
      g5Dphi(par->mass, Hchi, &chi[n]);
      start_sf_sendrecv(Hchi);
      #endif
      
#ifdef MEASURE_FORCERHMC
      /* reset force stat counters */
      double forcestat[2]={0.,0}; /* used for computation of avr and max force */
#endif
      
      _PIECE_FOR(&glattice,xp) {
        
        suNg_algebra_vector f;
        suNf_vector ptmp;
        suNf_spinor p;
        suNf_FMAT s1;

        if(xp==glattice.inner_master_pieces) {
          _OMP_PRAGMA( master )
          /* wait for spinor to be transfered */
          {
#ifdef UPDATE_EO
            complete_sf_sendrecv(&delta);
            complete_sf_sendrecv(&sigma);
#else
            complete_sf_sendrecv(Hchi);
#endif
          }
          _OMP_PRAGMA( barrier )
        }
        
#ifdef MEASURE_FORCEHMC
        _SITE_FOR_SUM(&glattice,xp,x,forcestat[0],forcestat[1]) {
#else
        _SITE_FOR(&glattice,xp,x) {
#endif
          
          for (int mu=0; mu<4; ++mu) {
            int y;
            suNf_spinor *chi1, *chi2;
            _suNf_FMAT_zero(s1);
            switch (mu) {
            case 0:
              y=iup(x,0);
              chi1=_FIELD_AT(Hchi,x);
              chi2=_FIELD_AT(&chi[n],y);
              _F_DIR0(s1,chi1,chi2);
              chi1=_FIELD_AT(&chi[n],x);
              chi2=_FIELD_AT(Hchi,y);
              _F_DIR0(s1,chi1,chi2);
              break;
            case 1:
              y=iup(x,1);
              chi1=_FIELD_AT(Hchi,x);
              chi2=_FIELD_AT(&chi[n],y);
              _F_DIR1(s1,chi1,chi2);
              chi1=_FIELD_AT(&chi[n],x);
              chi2=_FIELD_AT(Hchi,y);
              _F_DIR1(s1,chi1,chi2);
              break;
            case 2:
              y=iup(x,2);
              chi1=_FIELD_AT(Hchi,x);
              chi2=_FIELD_AT(&chi[n],y);
              _F_DIR2(s1,chi1,chi2);
              chi1=_FIELD_AT(&chi[n],x);
              chi2=_FIELD_AT(Hchi,y);
              _F_DIR2(s1,chi1,chi2);
              break;
            default: /* DIR 3 */
              y=iup(x,3);
              chi1=_FIELD_AT(Hchi,x);
              chi2=_FIELD_AT(&chi[n],y);
              _F_DIR3(s1,chi1,chi2);
              chi1=_FIELD_AT(&chi[n],x);
              chi2=_FIELD_AT(Hchi,y);
              _F_DIR3(s1,chi1,chi2);
            }
            
            _algebra_project(f,s1);
            /*_print_avect(f); */
            #ifdef UPDATE_EO
            _algebra_vector_mul_add_assign_g(*_4FIELD_AT(force,x,mu),-dt*ratio->a[n+1]*(_REPR_NORM2/_FUND_NORM2),f);	
            #else
            _algebra_vector_mul_add_assign_g(*_4FIELD_AT(force,x,mu),dt*ratio->a[n+1]*(_REPR_NORM2/_FUND_NORM2),f);	
            #endif

#ifdef MEASURE_FORCERHMC
            double nsq;
            _algebra_vector_sqnorm_g(nsq,f);
            forcestat[0]+=sqrt(nsq);
            for(y=0;y<NG*NG-1;++y){
              if(forcestat[1]<fabs(*(((double*)&f)+y))) forcestat[1]=fabs(*(((double*)&f)+y));
            }
#endif
          } //directions for
        } //SITE_FOR
      } //PIECE_FOR
        
        
#ifdef MEASURE_FORCERHMC
      if(logger_getlevel("FORCE-STAT")>=10){
        global_sum(forcestat,1);
        global_max(forcestat+1,1);
          
        forcestat[0]*=ratio->a[n+1]*(_REPR_NORM2/_FUND_NORM2)/(4.*GLB_VOLUME);
        forcestat[1]*=ratio->a[n+1]*(_REPR_NORM2/_FUND_NORM2);
        lprintf("FORCE-STAT",10," force_rhmc: dt= %1.8e avr |force|= %1.8e maxforce= %1.8e mass= %f k= %d n= %d a= %1.8e b= %1.8e\n",dt,forcestat[0],forcestat[1],par->mass,k,n,ratio->a[n+1],ratio->b[n]);
        }
#endif
        
      } //For n ratio order
    } //For pseudofermions loop

  apply_BCs_on_momentum_field(force);

  
  #ifdef UPDATE_EO
  chi[0].type=&glattice;
  #endif
  
  #ifdef TIMING
  #ifdef TIMING_WITH_BARRIERS
  MPI_Barrier(GLB_COMM);
  #endif
  gettimeofday(&end,0);
  timeval_subtract(&etime,&end,&start);
  lprintf("TIMING",0,"force_rhmc %.6f s\n",1.*etime.tv_sec+1.e-6*etime.tv_usec);
  #endif
}

#undef _F_DIR0
#undef _F_DIR1
#undef _F_DIR2
#undef _F_DIR3
