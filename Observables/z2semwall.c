/***************************************************************************\
* Copyright (c) 2009, Agostino Patella, Antonio Rago                        *   
* All rights reserved.                                                      * 
\***************************************************************************/

#include "global.h"
#include "linear_algebra.h"
#include "inverters.h"
#include "suN.h"
#include "observables.h"
#include "dirac.h"
#include "utils.h"
#include "memory.h"
#include "update.h"
#include "error.h"
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "logger.h"
#include "io.h"
#include "random.h"
#include "communications.h"
#include "ranlux.h"


#define RESTRICTED_SET



#ifdef RESTRICTED_SET
  #define G1_CHANNEL
  #define PCAC_CHANNEL
#else
  #define ID_CHANNEL
  #define G0_CHANNEL
  #define G1_CHANNEL
  #define G2_CHANNEL
  #define G3_CHANNEL
  #define G0G5_CHANNEL
  #define G5G1_CHANNEL
  #define G5G2_CHANNEL
  #define G5G3_CHANNEL
  #define G0G1_CHANNEL
  #define G0G2_CHANNEL
  #define G0G3_CHANNEL
  #define G0G5G1_CHANNEL
  #define G0G5G2_CHANNEL
  #define G0G5G3_CHANNEL
  #define PCAC_CHANNEL
#endif



/*static void create_diluted_source(spinor_field *source, int tau, int beta) {*/
/*  int c[4];*/
/*  spinor_field_zero_f(source);*/
/*  if(COORD[0]==tau/T) {*/
/*    c[0]=tau%T;*/
/*    for(c[1]=0; c[1]<X; c[1]++)*/
/*    for(c[2]=0; c[2]<Y; c[2]++)*/
/*    for(c[3]=0; c[3]<Z; c[3]++)*/
/*      ranz2((double*)(&(_FIELD_AT(source,ipt(c[0],c[1],c[2],c[3])))->c[beta]),sizeof(suNf_vector)/sizeof(double));*/
/*  }*/
/*}*/

static void create_diluted_source_even(spinor_field *source, int tau, int beta) {
  int c[4];
  spinor_field_zero_f(source);
  if(COORD[0]==tau/T) {
    c[0]=tau%T;
    for(c[1]=0; c[1]<X; c[1]++)
    for(c[2]=0; c[2]<Y; c[2]++)
    for(c[3]=0; c[3]<Z; c[3]++)
      if(((tau+COORD[1]*X+c[1]+COORD[2]*Y+c[2]+COORD[3]*Z+c[3])&1)==0)
        ranz2((double*)(&(_FIELD_AT(source,ipt(c[0],c[1],c[2],c[3])))->c[beta]),sizeof(suNf_vector)/sizeof(double));
  }
}


static double hmass, hmass_pre;

static void D(spinor_field *out, spinor_field *in){
  Dphi(hmass,out,in);
}

static void D_pre(spinor_field *out, spinor_field *in){
  Dphi_eopre(hmass_pre,out,in);
}


/*enum { ID, G0, G1, G2, G3, G5, G0G1, G0G2, G0G3, G0G5, G5G1, G5G2, G5G3, G0G5G1, G0G5G2, G0G5G3 };*/

static void z2semwall_qprop_QMR_eo(void (*Gamma)(suNf_spinor*,suNf_spinor*), spinor_field *psi, spinor_field *eta, int nm, double *mass, double acc) {
	mshift_par QMR_par;
	double *shift;
  spinor_field *eta2, *resd, *resdn;
  _DECLARE_INT_ITERATOR(ix);
 	double norm;
	spinor_field qprop_mask;
	int i, cgiter=0;

#ifndef NDEBUG
	spinor_field *test=0;
 	spinor_field *test_e=0;
	test=alloc_spinor_field_f(2,&glattice);
	test_e=alloc_spinor_field_f(1,&glat_even);
#endif

	eta2=alloc_spinor_field_f(2*nm+1,&glat_even);
	resdn=eta2+1;
	resd=resdn+nm;


	/* set up inverters parameters */
	shift=(double*)malloc(sizeof(double)*(nm));
	hmass_pre=mass[0]; /* we can put any number here!!! */
	for(i=0;i<nm;++i){
		shift[i]=(4.+hmass_pre)*(4.+hmass_pre)-(4.+mass[i])*(4.+mass[i]);
	}
	QMR_par.n = nm;
	QMR_par.shift = shift;
	QMR_par.err2 = .5*acc;
	QMR_par.max_iter = 0;

	/* noisy background */
	gaussian_spinor_field(eta2);
	norm=sqrt(spinor_field_sqnorm_f(eta2));
	spinor_field_mul_f(eta2,1./norm,eta2);

	/* invert noise */
	for(i=0;i<QMR_par.n;++i){
		spinor_field_zero_f(&resdn[i]);
	}
	cgiter+=g5QMR_mshift(&QMR_par, &D_pre, eta2, resdn);

	/* add source */
  _MASTER_FOR(&glat_even,ix)
    (*Gamma)(_FIELD_AT(eta2,ix),_FIELD_AT(eta,ix));

  /* invert source */
	for(i=0;i<QMR_par.n;++i){
	  		spinor_field_zero_f(&resd[i]);
	}
	cgiter+=g5QMR_mshift(&QMR_par, &D_pre, eta2, resd);

	for(i=0;i<QMR_par.n;++i){

#ifndef NDEBUG
		/* this is a test of the solution */
		hmass_pre=mass[i];
		D_pre(test_e,&resd[i]);
		++cgiter;
		spinor_field_sub_f(test_e,test_e,eta2);
		norm=spinor_field_sqnorm_f(test_e);
		lprintf("PROPAGATOR",0,"g5QMR_eo residuum of source [%d] = %e\n",i,norm);
		hmass_pre=mass[0];
#endif /* NDEBUG */

    /* substract noise */
		spinor_field_sub_f(&resd[i],&resd[i],&resdn[i]);

		/* compute solution */
		qprop_mask=psi[i];
		qprop_mask.type=&glat_even;
		spinor_field_mul_f(&qprop_mask,(4.+mass[i]),&resd[i]);
		qprop_mask.type=&glat_odd;
		Dphi_(&qprop_mask,&resd[i]);
		spinor_field_minus_f(&qprop_mask,&qprop_mask);
		if(i&1) ++cgiter; /* count only half of calls. works because the number of sources is even */

#ifndef NDEBUG
		/* this is a test of the solution */
		hmass=mass[i];
		D(&test[0],&psi[i]);
		++cgiter;
		spinor_field_zero_f(&test[1]);
    _MASTER_FOR(&glat_even,ix)
      (*Gamma)(_FIELD_AT(&test[1],ix),_FIELD_AT(eta,ix));
		spinor_field_sub_f(&test[0],&test[0],&test[1]);
		norm=spinor_field_sqnorm_f(&test[0]);
		lprintf("PROPAGATOR",0,"g5QMR_eo residuum of source [%d] = %e\n",i,norm);
		hmass=mass[0];
#endif /* NDEBUG */
  }

	lprintf("PROPAGATOR",10,"QMR_eo MVM = %d\n",cgiter);
	
	free_spinor_field(eta2);
	free(shift);
#ifndef NDEBUG
	free_spinor_field(test_e);
	free_spinor_field(test);
#endif
}


void z2semwall_mesons(int conf, int nm, double *mass, double acc) {
  spinor_field *eta;
  spinor_field *psi0;
  spinor_field *psi;
  suNf_spinor sp;
  int ix, i;
  int beta, tau;
  int t,x,y,z;
  double ran, tmp;
  double corr[GLB_T*nm];
  
  eta=alloc_spinor_field_f(4,&glat_even);
  psi0=alloc_spinor_field_f(5*nm,&glattice);
  psi=psi0+4*nm;
  
  ranlxd(&ran,1);
  tau=(int)(ran*GLB_T);

  for(beta=0;beta<4;beta++) {
    create_diluted_source_even(&eta[beta], tau, beta);
    z2semwall_qprop_QMR_eo(&g5_eval_g5GammaDag_times_spinor,&psi0[beta*nm],&eta[beta],nm,mass,acc);
  }

  for(i=0; i<nm*GLB_T; i++)
    corr[i] = 0.;

  for(beta=0;beta<4;beta++) {
    for(i=0; i<nm; i++) {  
      for (t=0; t<T; t++) {
        for (x=0; x<X; x++) for (y=0; y<Y; y++) for (z=0; z<Z; z++) {
          ix=ipt(t,x,y,z);
          _spinor_prod_re_f(tmp,*_FIELD_AT(&psi0[beta*nm+i],ix),*_FIELD_AT(&psi0[beta*nm+i],ix));
          corr[(COORD[0]*T+t+GLB_T-tau)%GLB_T+i*GLB_T]+=tmp;
        }
      }
    }
  }

  global_sum(corr,GLB_T);

  
#define PRINT_CORR(name) \
  for(i=0; i<nm; i++) { \
    lprintf("MAIN",0,"conf #%d mass=%2.6f TRIPLET " #name "= ",conf,mass[i]); \
    for(t=0;t<GLB_T;++t) { \
      lprintf("MAIN",0,"%e ",corr[t]); \
    } \
    lprintf("MAIN",0,"\n"); \
    fflush(stdout); \
  } \

  PRINT_CORR(g5);


#define COMPUTE_CORR(name) \
  for(i=0; i<nm*GLB_T; i++) \
    corr[i] = 0.; \
  for(beta=0;beta<4;beta++) { \
    z2semwall_qprop_QMR_eo(& name##_eval_g5GammaDag_times_spinor,psi,&eta[beta],nm,mass,acc); \
    for(i=0; i<nm; i++) { \
      for (t=0; t<T; t++) { \
        for (x=0; x<X; x++) for (y=0; y<Y; y++) for (z=0; z<Z; z++) { \
          ix=ipt(t,x,y,z); \
          _spinor_zero_f(sp); \
          name##_eval_g5GammaDag_times_spinor(&sp,_FIELD_AT(&psi0[beta*nm+i],ix)); \
          _spinor_prod_re_f(tmp,*_FIELD_AT(&psi[i],ix),sp); \
          corr[(COORD[0]*T+t+GLB_T-tau)%GLB_T+i*GLB_T] += tmp; \
        } \
      } \
    } \
  } \
  global_sum(corr,GLB_T)

#ifdef ID_CHANNEL
  COMPUTE_CORR(id);
  PRINT_CORR(id);
#endif

#ifdef G0_CHANNEL
  COMPUTE_CORR(g0);
  PRINT_CORR(g0);
#endif
  
#ifdef G1_CHANNEL
  COMPUTE_CORR(g1);
  for(i=0; i<nm*GLB_T; i++)
    corr[i] = -corr[i];
  PRINT_CORR(g1);
#endif
   
#ifdef G2_CHANNEL
  COMPUTE_CORR(g2);
  for(i=0; i<nm*GLB_T; i++)
    corr[i] = -corr[i];
  PRINT_CORR(g2);
#endif
  
#ifdef G3_CHANNEL
  COMPUTE_CORR(g3);
  for(i=0; i<nm*GLB_T; i++)
    corr[i] = -corr[i];
  PRINT_CORR(g3);
#endif
  
#ifdef G0G5_CHANNEL
  COMPUTE_CORR(g0g5);
  for(i=0; i<nm*GLB_T; i++)
    corr[i] = -corr[i];
  PRINT_CORR(g0g5);
#endif
  
#ifdef G5G1_CHANNEL
  COMPUTE_CORR(g5g1);
  PRINT_CORR(g5g1);
#endif
  
#ifdef G5G2_CHANNEL
  COMPUTE_CORR(g5g2);
  PRINT_CORR(g5g2);
#endif
  
#ifdef G5G3_CHANNEL
  COMPUTE_CORR(g5g3);
  PRINT_CORR(g5g3);
#endif
  
#ifdef G0G1_CHANNEL
  COMPUTE_CORR(g0g1);
  for(i=0; i<nm*GLB_T; i++)
    corr[i] = -corr[i];
  PRINT_CORR(g0g1);
#endif
  
#ifdef G0G2_CHANNEL
  COMPUTE_CORR(g0g2);
  for(i=0; i<nm*GLB_T; i++)
    corr[i] = -corr[i];
  PRINT_CORR(g0g2);
#endif
  
#ifdef G0G3_CHANNEL
  COMPUTE_CORR(g0g3);
  for(i=0; i<nm*GLB_T; i++)
    corr[i] = -corr[i];
  PRINT_CORR(g0g3);
#endif
  
#ifdef G0G5G1_CHANNEL
  COMPUTE_CORR(g0g5g1);
  PRINT_CORR(g0g5g1);
#endif
  
#ifdef G0G5G2_CHANNEL
  COMPUTE_CORR(g0g5g2);
  PRINT_CORR(g0g5g2);
#endif
  
#ifdef G0G5G3_CHANNEL
  COMPUTE_CORR(g0g5g3);
  PRINT_CORR(g0g5g3);
#endif


#ifdef PCAC_CHANNEL
  for(i=0; i<nm*GLB_T; i++)
    corr[i] = 0.;
  for(beta=0;beta<4;beta++) {
    for(i=0; i<nm; i++) {
      for (t=0; t<T; t++) {
        for (x=0; x<X; x++) for (y=0; y<Y; y++) for (z=0; z<Z; z++) {
          ix=ipt(t,x,y,z);
          _spinor_zero_f(sp);
          g0g5_eval_g5GammaDag_times_spinor(&sp,_FIELD_AT(&psi0[beta*nm+i],ix));
          _spinor_prod_re_f(tmp,*_FIELD_AT(&psi0[beta*nm+i],ix),sp);
          corr[(COORD[0]*T+t+GLB_T-tau)%GLB_T+i*GLB_T] += tmp;
        }
      }
    }
  }
  global_sum(corr,GLB_T);
  PRINT_CORR(g5_g0g5_re);
#endif


}
