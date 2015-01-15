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
#include "meson_observables.h"

/*#define POINT_TO_ALL*/

/* #define RESTRICTED_SET */



#ifdef RESTRICTED_SET
  #define G1_CHANNEL
  #define PCAC_CHANNEL
  #define NCHANNELS 3
// enum { _g5=0, _g1, _g5_g0g5_re };
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
  #define NCHANNELS 17
//  enum { _id=0, _g0, _g1, _g2, _g3, _g5, _g0g1, _g0g2, _g0g3, _g0g5, _g5g1, _g5g2, _g5g3, _g0g5g1, _g0g5g2, _g0g5g3, _g5_g0g5_re };
#endif

static const int _g5_g0g5_re=NCHANNELS-1;


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


#ifdef POINT_TO_ALL
static void create_point_source_even(spinor_field *source, int beta, int b) {
  spinor_field_zero_f(source);
  if(COORD[0]==0 && COORD[1]==0 && COORD[2]==0 && COORD[3]==0) {
    int ix=ipt(0,0,0,0);
    _FIELD_AT(source,ix)->c[beta].c[b].re = 1.;
  }
}
#else
static void create_diluted_source_even(spinor_field *source, int tau, int beta) {
  int c[4];
  spinor_field_zero_f(source);
  if(COORD[0]==tau/T) {
    c[0]=tau%T;
    for(c[1]=0; c[1]<X; c[1]++)
    for(c[2]=0; c[2]<Y; c[2]++)
    for(c[3]=0; c[3]<Z; c[3]++)
      if(((tau+zerocoord[1]+c[1]+zerocoord[2]+c[2]+zerocoord[3]+c[3])&1)==0)
        ranz2((double*)(&(_FIELD_AT(source,ipt(c[0],c[1],c[2],c[3])))->c[beta]),sizeof(suNf_vector)/sizeof(double));
  }
  
#ifndef NDEBUG
  double norm=spinor_field_sqnorm_f(source);
  lprintf("Z2SEMWALL",0,"Source sqnorm (must be %.0f) = %e\n",(.5*NF)*GLB_VOL3,norm);
#endif
}
#endif



#ifndef NDEBUG
static double hmass;
#endif
static double hmass_pre;


static void D_pre(spinor_field *out, spinor_field *in){
  Dphi_eopre(hmass_pre,out,in);
}



static int init=0;
static mshift_par QMR_par;
static double *shift;
static double *mass;
static spinor_field *QMR_noise;
static spinor_field *QMR_resdn;
#ifndef NDEBUG
static spinor_field *test;
static spinor_field *test_e;
#endif
static spinor_field *eta2;
static spinor_field *resd;
static spinor_field *eta;
static spinor_field *psi0;
static spinor_field *psi;
static void z2semwall_qprop_init(int nm, double *m, double acc) {
  int i, cgiter=0;
  double norm;

  if(init==0) {


    shift=(double*)malloc(sizeof(double)*(nm));
    mass=(double*)malloc(sizeof(double)*(nm));
    hmass_pre=m[0]; /* we can put any number here!!! */
#ifndef NDEBUG
    hmass=m[0]; /* we can put any number here!!! */
#endif
    for(i=0;i<nm;++i){
      mass[i]=m[i];
      shift[i]=(4.+hmass_pre)*(4.+hmass_pre)-(4.+m[i])*(4.+m[i]);
    }
    QMR_par.n = nm;
    QMR_par.shift = shift;
    QMR_par.err2 = .5*acc;
    QMR_par.max_iter = 0;
  
    QMR_noise=alloc_spinor_field_f(nm+1,&glat_even);
    QMR_resdn=QMR_noise+1;

#ifndef NDEBUG
    test=alloc_spinor_field_f(2,&glattice);
    test_e=alloc_spinor_field_f(1,&glat_even);
#endif
 
    eta2=alloc_spinor_field_f(QMR_par.n+1,&glat_even);
    resd=eta2+1;

    eta=alloc_spinor_field_f(4,&glat_even);
    psi0=alloc_spinor_field_f(5*nm,&glattice);
    psi=psi0+4*nm;

    /* noisy background */
    gaussian_spinor_field(QMR_noise);
    norm=sqrt(spinor_field_sqnorm_f(QMR_noise));
    spinor_field_mul_f(QMR_noise,1./norm,QMR_noise);
  }

  /* invert noise */
  for(i=0;i<QMR_par.n;++i)
    spinor_field_zero_f(&QMR_resdn[i]);
  cgiter+=g5QMR_mshift(&QMR_par, &D_pre, QMR_noise, QMR_resdn);

#ifndef NDEBUG
  for(i=0;i<QMR_par.n;++i){
    hmass_pre=mass[i];
    D_pre(test_e,&QMR_resdn[i]);
    ++cgiter;
    spinor_field_sub_assign_f(test_e,QMR_noise);
    norm=spinor_field_sqnorm_f(test_e);
    lprintf("Z2SEMWALL",0,"g5QMR_eo residuum of source [%d] = %e\n",i,norm);
    hmass_pre=mass[0];
  }
#endif /* NDEBUG */

  lprintf("Z2SEMWALL",10,"QMR_eo MVM = %d\n",cgiter);

  init=1;
}

void z2semwall_qprop_free() {
  error(init==0,1,"z2semwall.c","z2semwall method not initialized!");

#ifndef NDEBUG
  free_spinor_field_f(test_e);
  free_spinor_field_f(test);
#endif /* NDEBUG */

  free_spinor_field_f(eta2);

  free_spinor_field_f(eta);
  free_spinor_field_f(psi0);


  free(shift);
  free(mass);
    
  free_spinor_field_f(QMR_noise);
  init=0;
}


/***************************************************************************\

Gamma is supposed to be one of the name##_eval_g5GammaDag_times_spinor
functions in mesons.c

psi = D^{-1} g5 Gamma^+ eta

\***************************************************************************/

static void z2semwall_qprop_QMR_eo(void (*Gamma)(suNf_spinor*,suNf_spinor*), spinor_field *psi, spinor_field *eta) {

  spinor_field qprop_mask;
  int i, cgiter=0;

  error(init==0,1,"z2semwall.c","z2semwall method not initialized!");

#ifndef NDEBUG
    double norm;
#endif /* NDEBUG */

  for(i=0;i<QMR_par.n;++i){
#ifndef NDEBUG
    /* this is a test of the solution */
    hmass_pre=mass[i];
    D_pre(test_e,&QMR_resdn[i]);
    ++cgiter;
    spinor_field_sub_assign_f(test_e,QMR_noise);
    norm=spinor_field_sqnorm_f(test_e);
    lprintf("Z2SEMWALL",0,"g5QMR_eo residuum of source QMR_noise [%d] = %e\n",i,norm);
    hmass_pre=mass[0];
#endif /* NDEBUG */
  }


  /* add source */
  _MASTER_FOR(&glat_even,ix)
    (*Gamma)(_FIELD_AT(eta2,ix),_FIELD_AT(eta,ix));
  spinor_field_add_assign_f(eta2,QMR_noise);


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
    spinor_field_sub_assign_f(test_e,eta2);
    norm=spinor_field_sqnorm_f(test_e);
    lprintf("Z2SEMWALL",0,"g5QMR_eo residuum of source Gamma(eta)+QMR_noise [%d] = %e\n",i,norm);
    hmass_pre=mass[0];
#endif /* NDEBUG */

    /* substract noise */
    spinor_field_sub_assign_f(&resd[i],&QMR_resdn[i]);
#ifndef NDEBUG
    spinor_field_sub_assign_f(eta2,QMR_noise);
#endif /* NDEBUG */


#ifndef NDEBUG
    /* this is a test of the solution */
    hmass_pre=mass[i];
    D_pre(test_e,&resd[i]);
    ++cgiter;
    spinor_field_sub_assign_f(test_e,eta2);
    norm=spinor_field_sqnorm_f(test_e);
    lprintf("Z2SEMWALL",0,"g5QMR_eo residuum of source Gamma(eta) [%d] = %e\n",i,norm);
    hmass_pre=mass[0];
#endif /* NDEBUG */

    /* compute solution */
    qprop_mask=psi[i];
    qprop_mask.type=&glat_even;
    /* qprop_mask.ptr=psi[i].ptr+glat_even.master_shift; */
    spinor_field_mul_f(&qprop_mask,(4.+mass[i]),&resd[i]);
    qprop_mask.type=&glat_odd;
    qprop_mask.ptr=psi[i].ptr+glat_odd.master_shift; 
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
    spinor_field_sub_assign_f(&test[0],&test[1]);
    norm=spinor_field_sqnorm_f(&test[0]);
    lprintf("Z2SEMWALL",0,"g5QMR_eo residuum of source Gamma(eta) (2) [%d] = %e\n",i,norm);
    hmass=mass[0];
#endif /* NDEBUG */
  }

  lprintf("Z2SEMWALL",10,"QMR_eo MVM = %d\n",cgiter);
}



/***************************************************************************\

Computes

          1
C(t) = - ---- sum < tr[ \bar{Gamma} D^{-1}(x,t;0) Gamma D^{-1}(0;x,t) ] >
          V3   x

following 0804.1501


ifdef POINT_TO_ALL

a,b = spin indices
i,j = color indices

eta^{(b,j)}_{a,i}(x,t) = \delta_{ab} \delta_{ij} \delta_{x0} \delta_{t0}
psi0^{(b,j)} = D^{-1} eta^{(b,j)}
psi^{(b,j)} = D^{-1} g5 Gamma^dag eta^{(b,j)}
corr[t] = - s/L^3 \sum_{x,b,j}
            [psi^{(b,j)}(x,t)]^dag g5 Gamma^dag psi0^{(b,j)}(x,t)
where s is the sign such that
s Gamma^dag = g0 Gamma^dag g0

else * POINT_TO_ALL *

a,b = spin indices
h = 1, ..., Nh
t0 = random timeslice

eta^{(b,h)}_{a,i}(x,t) = z2(i,x) \delta_{ab} \chi_EVEN(x) \delta_{t,t0}
where z2(i,x) are independent z2 random numbers

psi0^{(b,j)} = D^{-1} eta^{(b,j)}
psi^{(b,j)} = D^{-1} g5 Gamma^dag eta^{(b,j)}
corr[t] = - 2s/(L^6*Nh) \sum_{x,b,j}
            [psi^{(b,j)}(x,t)]^dag g5 Gamma^dag psi0^{(b,j)}(x,t)
where s is the sign such that
s Gamma^dag = g0 Gamma^dag g0

endif * POINT_TO_ALL *

\***************************************************************************/

void z2semwall_mesons(int conf, int nhits, int nm, double *m, double acc) {
  suNf_spinor sp;
  int ix, i, k, n;
  int beta, tau;
#ifndef POINT_TO_ALL
  double ran;
#endif
  int t,x,y,z;
  double tmp;
  double corr[NCHANNELS][GLB_T*nm];

  error(nhits<1,1,"z2semwall.c","Bad value for nhits!");

  z2semwall_qprop_init(nm, m, acc);

  for(i=0; i<nm*GLB_T; i++)
  for(k=0; k<NCHANNELS; k++)
    corr[k][i] = 0.;

#ifdef POINT_TO_ALL
  nhits=NF;
#endif

  for(n=0; n<nhits; n++) {
  
#ifdef POINT_TO_ALL
    tau=0;
#else
    do{
      ranlxd(&ran,1);
      tau=(int)(ran*GLB_T);
    } while(tau==GLB_T);
    bcast_int(&tau,1);
#endif

    for(beta=0;beta<4;beta++) {
#ifdef POINT_TO_ALL
      create_point_source_even(&eta[beta], beta, n);
#else
      create_diluted_source_even(&eta[beta], tau, beta);
#endif
      z2semwall_qprop_QMR_eo(&g5_eval_g5GammaDag_times_spinor,&psi0[beta*nm],&eta[beta]);
    }


    for(beta=0;beta<4;beta++)
      for(i=0; i<nm; i++)
        for (t=0; t<T; t++)
          for (x=0; x<X; x++)
            for (y=0; y<Y; y++) 
              for (z=0; z<Z; z++) {
                ix=ipt(t,x,y,z);
                _spinor_prod_re_f(tmp,*_FIELD_AT(&psi0[beta*nm+i],ix),*_FIELD_AT(&psi0[beta*nm+i],ix));
                corr[_g5][(zerocoord[0]+t+GLB_T-tau)%GLB_T+i*GLB_T]+=tmp;
              }
    

#define COMPUTE_CORR(name) \
    for(beta=0;beta<4;beta++) { \
      z2semwall_qprop_QMR_eo(& name##_eval_g5GammaDag_times_spinor,psi,&eta[beta]); \
      for(i=0; i<nm; i++) { \
        for (t=0; t<T; t++) { \
          for (x=0; x<X; x++) for (y=0; y<Y; y++) for (z=0; z<Z; z++) { \
            ix=ipt(t,x,y,z); \
            _spinor_zero_f(sp); \
            name##_eval_g5GammaDag_times_spinor(&sp,_FIELD_AT(&psi0[beta*nm+i],ix)); \
            _spinor_prod_re_f(tmp,*_FIELD_AT(&psi[i],ix),sp); \
            corr[ _##name ][(zerocoord[0]+t+GLB_T-tau)%GLB_T+i*GLB_T] += tmp; \
          } \
        } \
      } \
    }

#ifdef ID_CHANNEL
    COMPUTE_CORR(id);
#endif

#ifdef G0_CHANNEL
    COMPUTE_CORR(g0);
#endif
  
#ifdef G1_CHANNEL
    COMPUTE_CORR(g1);
#endif
   
#ifdef G2_CHANNEL
    COMPUTE_CORR(g2);
#endif
  
#ifdef G3_CHANNEL
    COMPUTE_CORR(g3);
#endif
  
#ifdef G0G5_CHANNEL
    COMPUTE_CORR(g0g5);
#endif
  
#ifdef G5G1_CHANNEL
    COMPUTE_CORR(g5g1);
#endif
  
#ifdef G5G2_CHANNEL
    COMPUTE_CORR(g5g2);
#endif
  
#ifdef G5G3_CHANNEL
    COMPUTE_CORR(g5g3);
#endif
  
#ifdef G0G1_CHANNEL
    COMPUTE_CORR(g0g1);
#endif
  
#ifdef G0G2_CHANNEL
    COMPUTE_CORR(g0g2);
#endif
  
#ifdef G0G3_CHANNEL
    COMPUTE_CORR(g0g3);
#endif
  
#ifdef G0G5G1_CHANNEL
    COMPUTE_CORR(g0g5g1);
#endif
  
#ifdef G0G5G2_CHANNEL
    COMPUTE_CORR(g0g5g2);
#endif
  
#ifdef G0G5G3_CHANNEL
    COMPUTE_CORR(g0g5g3);
#endif

#ifdef PCAC_CHANNEL
     for(beta=0;beta<4;beta++) {
      for(i=0; i<nm; i++) {
        for (t=0; t<T; t++) {
          for (x=0; x<X; x++) for (y=0; y<Y; y++) for (z=0; z<Z; z++) {
            ix=ipt(t,x,y,z);
            _spinor_zero_f(sp);
            g0g5_eval_g5GammaDag_times_spinor(&sp,_FIELD_AT(&psi0[beta*nm+i],ix));
            _spinor_prod_re_f(tmp,*_FIELD_AT(&psi0[beta*nm+i],ix),sp);
            corr[_g5_g0g5_re][(zerocoord[0]+t+GLB_T-tau)%GLB_T+i*GLB_T] += tmp;
          }
        }
      }
    }
#endif

  }


  for(k=0; k<NCHANNELS; k++) {
    global_sum(corr[k],GLB_T*nm);
    for(i=0; i<nm*GLB_T; i++)
#ifdef POINT_TO_ALL
      corr[k][i] *= -1./GLB_VOL3;
#else
      corr[k][i] *= -((2./nhits)/GLB_VOL3)/GLB_VOL3;
#endif
  }

  
#define PRINT_CORR(name) \
  for(i=0; i<nm; i++) { \
    lprintf("MAIN",0,"conf #%d mass=%2.6f TRIPLET " #name "= ",conf,mass[i]); \
    for(t=0;t<GLB_T;++t) { \
      lprintf("MAIN",0,"%e ",corr[ _##name ][t+i*GLB_T]); \
    } \
    lprintf("MAIN",0,"\n"); \
    fflush(stdout); \
  } \

  for(i=0; i<nm*GLB_T; i++)
    corr[_g5][i] *= -1.;
  PRINT_CORR(g5);

#ifdef ID_CHANNEL
  PRINT_CORR(id);
#endif

#ifdef G0_CHANNEL
  PRINT_CORR(g0);
#endif
  
#ifdef G1_CHANNEL
  for(i=0; i<nm*GLB_T; i++)
    corr[_g1][i] *= -1.;
  PRINT_CORR(g1);
#endif
   
#ifdef G2_CHANNEL
  for(i=0; i<nm*GLB_T; i++)
    corr[_g2][i] *= -1.;
  PRINT_CORR(g2);
#endif
  
#ifdef G3_CHANNEL
  for(i=0; i<nm*GLB_T; i++)
    corr[_g3][i] *= -1.;
  PRINT_CORR(g3);
#endif
  
#ifdef G0G5_CHANNEL
  for(i=0; i<nm*GLB_T; i++)
    corr[_g0g5][i] *= -1.;
  PRINT_CORR(g0g5);
#endif
  
#ifdef G5G1_CHANNEL
  PRINT_CORR(g5g1);
#endif
  
#ifdef G5G2_CHANNEL
  PRINT_CORR(g5g2);
#endif
  
#ifdef G5G3_CHANNEL
  PRINT_CORR(g5g3);
#endif
  
#ifdef G0G1_CHANNEL
  for(i=0; i<nm*GLB_T; i++)
    corr[_g0g1][i] *= -1.;
  PRINT_CORR(g0g1);
#endif
  
#ifdef G0G2_CHANNEL
  for(i=0; i<nm*GLB_T; i++)
    corr[_g0g2][i] *= -1.;
  PRINT_CORR(g0g2);
#endif
  
#ifdef G0G3_CHANNEL
  for(i=0; i<nm*GLB_T; i++)
    corr[_g0g3][i] *= -1.;
  PRINT_CORR(g0g3);
#endif
  
#ifdef G0G5G1_CHANNEL
  PRINT_CORR(g0g5g1);
#endif
  
#ifdef G0G5G2_CHANNEL
  PRINT_CORR(g0g5g2);
#endif
  
#ifdef G0G5G3_CHANNEL
  PRINT_CORR(g0g5g3);
#endif

#ifdef PCAC_CHANNEL
  PRINT_CORR(g5_g0g5_re);
#endif

}
