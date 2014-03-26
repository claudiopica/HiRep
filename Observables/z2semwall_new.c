/***************************************************************************\
* Copyright (c) 2009, Agostino Patella, Antonio Rago                        *
* All rights reserved.                                                      * 
*                                                                           *
*   Modified by Rudy Arthur, Ari Hietanen                                   *
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
#include "gamma_spinor.h"
#include "meson_observables.h"




/* #define RESTRICTED_SET */

/*  Introduces small gaussian noise for source field 
    Gives smaller errors with one extra inversion */
#define GAUSSIAN_NOISE 

#ifdef RESTRICTED_SET
  #define G1_CHANNEL
  #define PCAC_CHANNEL
  #define NCHANNELS 3
  enum { _g5=0, _g1, _g5_g0g5_re };
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
  #define ID_DISC_CHANNEL
  #define G5_DISC_CHANNEL
  #define NCHANNELS 19
  enum {_g5_g0g5_re=NCHANNELS-3, _id_disc,_g5_disc };
#endif


static void create_diluted_source_equal_even(spinor_field *source, int tau) {
  int c[4];
  suNf_vector *v1,*v2;
  int i;
  for (i=0;i<4;++i){
    spinor_field_zero_f(&source[i]);
  }
  
  if(COORD[0]==tau/T) {// Check that tau is in this thread.
    c[0]=tau%T;
    for(c[1]=0; c[1]<X; c[1]++) for(c[2]=0; c[2]<Y; c[2]++)  for(c[3]=0; c[3]<Z; c[3]++){
	  if(((tau+zerocoord[1]+c[1]+zerocoord[2]+c[2]+zerocoord[3]+c[3])&1)==0){
	    v1 = &((_FIELD_AT(&source[0],ipt(c[0],c[1],c[2],c[3])))->c[0]);
	    ranz2((double*)(v1),sizeof(suNf_vector)/sizeof(double)); // Make new sources
	    for (i=1;i<4;++i){
	      v2 = &((_FIELD_AT(&source[i],ipt(c[0],c[1],c[2],c[3])))->c[i]); //Copy previous index.
	      *v2 = *v1;
	    }
	  }
	}
  }
}


static double hmass_pre;


static void D_pre(spinor_field *out, spinor_field *in){
  Dphi_eopre(hmass_pre,out,in);
}

static int init=0;
static mshift_par QMR_par;
static double *shift;
static double *mass;
#ifdef GAUSSIAN_NOISE
static spinor_field *QMR_noise;
static spinor_field *QMR_resdn;
#endif
static spinor_field *resd;
static spinor_field *eta;
static spinor_field *eta2;
static spinor_field *psi0;
static spinor_field *psi;



static void z2semwall_qprop_init(int nm, double *m, double acc) {
  int i, cgiter=0;
#ifdef GAUSSIAN_NOISE
  double norm;
#endif

  if(init==0) {


    shift=(double*)malloc(sizeof(double)*(nm));
    mass=(double*)malloc(sizeof(double)*(nm));
    hmass_pre=m[0]; /* we can put any number here!!! */
    for(i=0;i<nm;++i){
      mass[i]=m[i];
      shift[i]=(4.+hmass_pre)*(4.+hmass_pre)-(4.+m[i])*(4.+m[i]);
    }
    QMR_par.n = nm;
    QMR_par.shift = shift;
    QMR_par.err2 = .5*acc;
    QMR_par.max_iter = 0;

    resd=alloc_spinor_field_f(QMR_par.n,&glat_even);

    eta=alloc_spinor_field_f(5,&glat_even);
    eta2=eta+4;
    psi0=alloc_spinor_field_f(5*nm,&glattice);
    psi=psi0+4*nm;

#ifdef GAUSSIAN_NOISE
    QMR_noise=alloc_spinor_field_f(nm+1,&glat_even);
    QMR_resdn=QMR_noise+1;
    /* noisy background */
    gaussian_spinor_field(QMR_noise);
    norm=sqrt(spinor_field_sqnorm_f(QMR_noise));
    spinor_field_mul_f(QMR_noise,1./norm,QMR_noise);
#endif
  }
#ifdef GAUSSIAN_NOISE
  /* invert noise */
  for(i=0;i<QMR_par.n;++i) spinor_field_zero_f(&QMR_resdn[i]);
  cgiter+=g5QMR_mshift(&QMR_par, &D_pre, QMR_noise, QMR_resdn);
#endif

  lprintf("Z2SEMWALL",10,"QMR_eo MVM = %d\n",cgiter);

  init=1;
}

void z2semwall_qprop_free_new() {
  error(init==0,1,"z2semwall.c","z2semwall method not initialized!");

  free_spinor_field_f(eta);
  free_spinor_field_f(psi0);
  free_spinor_field_f(resd);

  free(shift);
  free(mass);
  
#ifdef GAUSSIAN_NOISE
  free_spinor_field_f(QMR_noise);
#endif
  init=0;
}


/***************************************************************************\

psi = D^{-1} eta

\***************************************************************************/


static void z2semwall_qprop_QMR_eo(spinor_field *psi, spinor_field *eta) {
  spinor_field qprop_mask;
  int i, cgiter=0;
  lprintf("ZSEMWALL",0,"%g\n",spinor_field_sqnorm_f(eta));
  error(init==0,1,"z2semwall.c","z2semwall method not initialized!");

  /* add source */
#ifdef GAUSSIAN_NOISE
  spinor_field_add_f(eta2,eta,QMR_noise);
#else
  spinor_field_copy_f(eta2,eta);
#endif


  /* invert source */
  for(i=0;i<QMR_par.n;++i){
    spinor_field_zero_f(&resd[i]);
  }
  cgiter+=g5QMR_mshift(&QMR_par, &D_pre, eta2, resd);

  for(i=0;i<QMR_par.n;++i){
#ifdef GAUSSIAN_NOISE
    spinor_field_sub_assign_f(&resd[i],&QMR_resdn[i]);
#endif
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
  }
  
  lprintf("ZSEMWALL",0,"%g\n",spinor_field_sqnorm_f(eta));
  lprintf("Z2SEMWALL NEW",10,"QMR_eo MVM = %d\n",cgiter);
}


/***************************************************************************\

Computes

          1
C(t) = - ---- sum < tr[ \bar{Gamma} D^{-1}(x,t;0) Gamma D^{-1}(0;x,t) ] >
          V3   x

following 0804.1501

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


\***************************************************************************/

void z2semwall_mesons_new(int conf, int nhits, int nm, double *m, double acc) {
  int ix, i, k, n;
  int beta, tau;
  double ran;
  int t,x,y,z;
  double tmp;
  double corr[NCHANNELS][GLB_T*nm];
  suNf_spinor sp1,sp2;
  suNf_spinor *psi0p;
  int slices[GLB_T];
  int counter = 0;
  int itmp;
  
  error(nhits<1,1,"z2semwall_new.c","Bad value for nhits!");

  z2semwall_qprop_init(nm, m, acc);

  for(i=0; i<nm*GLB_T; i++)
    for(k=0; k<NCHANNELS; k++)
      corr[k][i] = 0.;

  for(n=0; n<nhits; n++) {
    /* Random timeslize not previously chosen */
    if (counter == 0){
      for (i=0;i<GLB_T;++i){
	slices[i]=i;
      }
      counter=GLB_T;
    }
    do{
      ranlxd(&ran,1);
      itmp=(int)(ran*counter);    
    } while(itmp==counter);
    counter--;
    tau = slices[itmp];
    slices[itmp]=slices[counter];
    slices[counter]=tau;

    bcast_int(&tau,1);

    create_diluted_source_equal_even(eta, tau);

    for(beta=0;beta<4;beta++) {
      z2semwall_qprop_QMR_eo(&psi0[beta*nm],&eta[beta]);
    }

    for(i=0; i<nm; i++) {						
      for (t=0; t<T; t++) {						
	for (x=0; x<X; x++) for (y=0; y<Y; y++) for (z=0; z<Z; z++) { 
	      ix=ipt(t,x,y,z);					
	      for (beta=0;beta<4;beta++){
		_spinor_prod_re_f(tmp,*_FIELD_AT(&psi0[beta*nm+i],ix),*_FIELD_AT(&psi0[beta*nm+i],ix));
		corr[_g5][(zerocoord[0]+t+GLB_T-tau)%GLB_T+i*GLB_T]+=tmp;
	    }

#ifdef ID_CHANNEL
	      psi0p = _FIELD_AT(&psi0[i],ix);
	      _spinor_g5_f(sp1,*psi0p);
	      sp2=*psi0p;
	      _spinor_prod_re_f(tmp,sp1,sp2);
	      corr[_id][(zerocoord[0]+t+GLB_T-tau)%GLB_T+i*GLB_T] += tmp; 

	      psi0p = _FIELD_AT(&psi0[nm+i],ix);
	      _spinor_g5_f(sp1,*psi0p);
	      sp2=*psi0p;
	      _spinor_prod_re_f(tmp,sp1,sp2);	
	      corr[_id][(zerocoord[0]+t+GLB_T-tau)%GLB_T+i*GLB_T] += tmp;

	      psi0p = _FIELD_AT(&psi0[2*nm+i],ix);
	      _spinor_g5_f(sp1,*psi0p);
	      _spinor_minus_f(sp2,*psi0p);
	      _spinor_prod_re_f(tmp,sp1,sp2);	
	      corr[_id][(zerocoord[0]+t+GLB_T-tau)%GLB_T+i*GLB_T] += tmp; 

	      psi0p = _FIELD_AT(&psi0[3*nm+i],ix);
	      _spinor_g5_f(sp1,*psi0p);
	      _spinor_minus_f(sp2,*psi0p);
	      _spinor_prod_re_f(tmp,sp1,sp2);	
	      corr[_id][(zerocoord[0]+t+GLB_T-tau)%GLB_T+i*GLB_T] += tmp; 
#endif
#ifdef G0_CHANNEL
	      psi0p = _FIELD_AT(&psi0[i],ix);
	      _spinor_g5g0_f(sp1,*psi0p);
	      psi0p = _FIELD_AT(&psi0[2*nm+i],ix);
	      sp2 = *psi0p;
	      _spinor_prod_re_f(tmp,sp1,sp2);	
	      corr[_g0][(zerocoord[0]+t+GLB_T-tau)%GLB_T+i*GLB_T] += tmp; 

	      psi0p = _FIELD_AT(&psi0[nm+i],ix);
	      _spinor_g5g0_f(sp1,*psi0p);
	      psi0p = _FIELD_AT(&psi0[3*nm+i],ix);
	      sp2 = *psi0p;
	      _spinor_prod_re_f(tmp,sp1,sp2);	
	      corr[_g0][(zerocoord[0]+t+GLB_T-tau)%GLB_T+i*GLB_T] += tmp;

	      psi0p = _FIELD_AT(&psi0[2*nm+i],ix);
	      _spinor_g5g0_f(sp1,*psi0p);
	      psi0p = _FIELD_AT(&psi0[0*nm+i],ix);
	      _spinor_minus_f(sp2,*psi0p);
	      _spinor_prod_re_f(tmp,sp1,sp2);	
	      corr[_g0][(zerocoord[0]+t+GLB_T-tau)%GLB_T+i*GLB_T] += tmp; 

	      psi0p = _FIELD_AT(&psi0[3*nm+i],ix);
	      _spinor_g5g0_f(sp1,*psi0p);
	      psi0p = _FIELD_AT(&psi0[1*nm+i],ix);
	      _spinor_minus_f(sp2,*psi0p);
	      _spinor_prod_re_f(tmp,sp1,sp2);	
	      corr[_g0][(zerocoord[0]+t+GLB_T-tau)%GLB_T+i*GLB_T] += tmp; 
#endif

#ifdef G1_CHANNEL
	      psi0p = _FIELD_AT(&psi0[i],ix);
	      _spinor_g5g1_f(sp1,*psi0p);
	      psi0p = _FIELD_AT(&psi0[3*nm+i],ix);
	      _spinor_i_plus_f(sp2,*psi0p);
	      _spinor_prod_re_f(tmp,sp1,sp2);	
	      corr[_g1][(zerocoord[0]+t+GLB_T-tau)%GLB_T+i*GLB_T] += tmp; 

	      psi0p = _FIELD_AT(&psi0[nm+i],ix);
	      _spinor_g5g1_f(sp1,*psi0p);
	      psi0p = _FIELD_AT(&psi0[2*nm+i],ix);
	      _spinor_i_plus_f(sp2,*psi0p);
	      _spinor_prod_re_f(tmp,sp1,sp2);	
	      corr[_g1][(zerocoord[0]+t+GLB_T-tau)%GLB_T+i*GLB_T] += tmp;

	      psi0p = _FIELD_AT(&psi0[2*nm+i],ix);
	      _spinor_g5g1_f(sp1,*psi0p);
	      psi0p = _FIELD_AT(&psi0[1*nm+i],ix);
	      _spinor_i_plus_f(sp2,*psi0p);
	      _spinor_prod_re_f(tmp,sp1,sp2);	
	      corr[_g1][(zerocoord[0]+t+GLB_T-tau)%GLB_T+i*GLB_T] += tmp; 

	      psi0p = _FIELD_AT(&psi0[3*nm+i],ix);
	      _spinor_g5g1_f(sp1,*psi0p);
	      psi0p = _FIELD_AT(&psi0[0*nm+i],ix);
	      _spinor_i_plus_f(sp2,*psi0p);
	      _spinor_prod_re_f(tmp,sp1,sp2);	
	      corr[_g1][(zerocoord[0]+t+GLB_T-tau)%GLB_T+i*GLB_T] += tmp; 	      
#endif	      
#ifdef G2_CHANNEL
	      psi0p = _FIELD_AT(&psi0[i],ix);
	      _spinor_g5g2_f(sp1,*psi0p);
	      psi0p = _FIELD_AT(&psi0[3*nm+i],ix);
	      _spinor_minus_f(sp2,*psi0p);
	      _spinor_prod_re_f(tmp,sp1,sp2);	
	      corr[_g2][(zerocoord[0]+t+GLB_T-tau)%GLB_T+i*GLB_T] += tmp; 

	      psi0p = _FIELD_AT(&psi0[nm+i],ix);
	      _spinor_g5g2_f(sp1,*psi0p);
	      psi0p = _FIELD_AT(&psi0[2*nm+i],ix);
	      sp2 = *psi0p;
	      _spinor_prod_re_f(tmp,sp1,sp2);	
	      corr[_g2][(zerocoord[0]+t+GLB_T-tau)%GLB_T+i*GLB_T] += tmp;

	      psi0p = _FIELD_AT(&psi0[2*nm+i],ix);
	      _spinor_g5g2_f(sp1,*psi0p);
	      psi0p = _FIELD_AT(&psi0[1*nm+i],ix);
	      _spinor_minus_f(sp2,*psi0p);
	      _spinor_prod_re_f(tmp,sp1,sp2);	
	      corr[_g2][(zerocoord[0]+t+GLB_T-tau)%GLB_T+i*GLB_T] += tmp; 

	      psi0p = _FIELD_AT(&psi0[3*nm+i],ix);
	      _spinor_g5g2_f(sp1,*psi0p);
	      psi0p = _FIELD_AT(&psi0[0*nm+i],ix);
	      sp2 = *psi0p;
	      _spinor_prod_re_f(tmp,sp1,sp2);	
	      corr[_g2][(zerocoord[0]+t+GLB_T-tau)%GLB_T+i*GLB_T] += tmp; 
#endif

#ifdef G3_CHANNEL
	      psi0p = _FIELD_AT(&psi0[i],ix);
	      _spinor_g5g3_f(sp1,*psi0p);
	      psi0p = _FIELD_AT(&psi0[2*nm+i],ix);
	      _spinor_i_plus_f(sp2,*psi0p);
	      _spinor_prod_re_f(tmp,sp1,sp2);	
	      corr[_g3][(zerocoord[0]+t+GLB_T-tau)%GLB_T+i*GLB_T] += tmp; 

	      psi0p = _FIELD_AT(&psi0[nm+i],ix);
	      _spinor_g5g3_f(sp1,*psi0p);
	      psi0p = _FIELD_AT(&psi0[3*nm+i],ix);
	      _spinor_i_minus_f(sp2,*psi0p);
	      _spinor_prod_re_f(tmp,sp1,sp2);	
	      corr[_g3][(zerocoord[0]+t+GLB_T-tau)%GLB_T+i*GLB_T] += tmp;

	      psi0p = _FIELD_AT(&psi0[2*nm+i],ix);
	      _spinor_g5g3_f(sp1,*psi0p);
	      psi0p = _FIELD_AT(&psi0[0*nm+i],ix);
	      _spinor_i_plus_f(sp2,*psi0p);
	      _spinor_prod_re_f(tmp,sp1,sp2);	
	      corr[_g3][(zerocoord[0]+t+GLB_T-tau)%GLB_T+i*GLB_T] += tmp; 

	      psi0p = _FIELD_AT(&psi0[3*nm+i],ix);
	      _spinor_g5g3_f(sp1,*psi0p);
	      psi0p = _FIELD_AT(&psi0[1*nm+i],ix);
	      _spinor_i_minus_f(sp2,*psi0p);
	      _spinor_prod_re_f(tmp,sp1,sp2);	
	      corr[_g3][(zerocoord[0]+t+GLB_T-tau)%GLB_T+i*GLB_T] += tmp; 
#endif

#ifdef G0G5_CHANNEL
	      psi0p = _FIELD_AT(&psi0[i],ix);
	      _spinor_g0_f(sp1,*psi0p);
	      psi0p = _FIELD_AT(&psi0[2*nm+i],ix);
	      sp2 = *psi0p;
	      _spinor_prod_re_f(tmp,sp1,sp2);	
	      corr[_g0g5][(zerocoord[0]+t+GLB_T-tau)%GLB_T+i*GLB_T] += tmp; 

	      psi0p = _FIELD_AT(&psi0[nm+i],ix);
	      _spinor_g0_f(sp1,*psi0p);
	      psi0p = _FIELD_AT(&psi0[3*nm+i],ix);
	      sp2 = *psi0p;
	      _spinor_prod_re_f(tmp,sp1,sp2);	
	      corr[_g0g5][(zerocoord[0]+t+GLB_T-tau)%GLB_T+i*GLB_T] += tmp;

	      psi0p = _FIELD_AT(&psi0[2*nm+i],ix);
	      _spinor_g0_f(sp1,*psi0p);
	      psi0p = _FIELD_AT(&psi0[0*nm+i],ix);
	      sp2 = *psi0p;
	      _spinor_prod_re_f(tmp,sp1,sp2);	
	      corr[_g0g5][(zerocoord[0]+t+GLB_T-tau)%GLB_T+i*GLB_T] += tmp; 

	      psi0p = _FIELD_AT(&psi0[3*nm+i],ix);
	      _spinor_g0_f(sp1,*psi0p);
	      psi0p = _FIELD_AT(&psi0[1*nm+i],ix);
	      sp2 = *psi0p;
	      _spinor_prod_re_f(tmp,sp1,sp2);	
	      corr[_g0g5][(zerocoord[0]+t+GLB_T-tau)%GLB_T+i*GLB_T] += tmp; 
#endif

#ifdef G5G1_CHANNEL
	      psi0p = _FIELD_AT(&psi0[i],ix);
	      _spinor_g1_f(sp1,*psi0p);
	      psi0p = _FIELD_AT(&psi0[3*nm+i],ix);
	      _spinor_i_plus_f(sp2,*psi0p);
	      _spinor_prod_re_f(tmp,sp1,sp2);	
	      corr[_g5g1][(zerocoord[0]+t+GLB_T-tau)%GLB_T+i*GLB_T] += tmp; 

	      psi0p = _FIELD_AT(&psi0[nm+i],ix);
	      _spinor_g1_f(sp1,*psi0p);
	      psi0p = _FIELD_AT(&psi0[2*nm+i],ix);
	      _spinor_i_plus_f(sp2,*psi0p);
	      _spinor_prod_re_f(tmp,sp1,sp2);	
	      corr[_g5g1][(zerocoord[0]+t+GLB_T-tau)%GLB_T+i*GLB_T] += tmp;

	      psi0p = _FIELD_AT(&psi0[2*nm+i],ix);
	      _spinor_g1_f(sp1,*psi0p);
	      psi0p = _FIELD_AT(&psi0[1*nm+i],ix);
	      _spinor_i_minus_f(sp2,*psi0p);
	      _spinor_prod_re_f(tmp,sp1,sp2);	
	      corr[_g5g1][(zerocoord[0]+t+GLB_T-tau)%GLB_T+i*GLB_T] += tmp; 

	      psi0p = _FIELD_AT(&psi0[3*nm+i],ix);
	      _spinor_g1_f(sp1,*psi0p);
	      psi0p = _FIELD_AT(&psi0[0*nm+i],ix);
	      _spinor_i_minus_f(sp2,*psi0p);
	      _spinor_prod_re_f(tmp,sp1,sp2);	
	      corr[_g5g1][(zerocoord[0]+t+GLB_T-tau)%GLB_T+i*GLB_T] += tmp; 	      
#endif	      
#ifdef G5G2_CHANNEL
	      psi0p = _FIELD_AT(&psi0[i],ix);
	      _spinor_g2_f(sp1,*psi0p);
	      psi0p = _FIELD_AT(&psi0[3*nm+i],ix);
	      _spinor_minus_f(sp2,*psi0p);
	      _spinor_prod_re_f(tmp,sp1,sp2);	
	      corr[_g5g2][(zerocoord[0]+t+GLB_T-tau)%GLB_T+i*GLB_T] += tmp; 

	      psi0p = _FIELD_AT(&psi0[nm+i],ix);
	      _spinor_g2_f(sp1,*psi0p);
	      psi0p = _FIELD_AT(&psi0[2*nm+i],ix);
	      sp2 = *psi0p;
	      _spinor_prod_re_f(tmp,sp1,sp2);	
	      corr[_g5g2][(zerocoord[0]+t+GLB_T-tau)%GLB_T+i*GLB_T] += tmp;

	      psi0p = _FIELD_AT(&psi0[2*nm+i],ix);
	      _spinor_g2_f(sp1,*psi0p);
	      psi0p = _FIELD_AT(&psi0[1*nm+i],ix);
	      sp2 = *psi0p;
	      _spinor_prod_re_f(tmp,sp1,sp2);	
	      corr[_g5g2][(zerocoord[0]+t+GLB_T-tau)%GLB_T+i*GLB_T] += tmp; 

	      psi0p = _FIELD_AT(&psi0[3*nm+i],ix);
	      _spinor_g2_f(sp1,*psi0p);
	      psi0p = _FIELD_AT(&psi0[0*nm+i],ix);
	      _spinor_minus_f(sp2,*psi0p);
	      _spinor_prod_re_f(tmp,sp1,sp2);	
	      corr[_g5g2][(zerocoord[0]+t+GLB_T-tau)%GLB_T+i*GLB_T] += tmp; 
#endif

#ifdef G5G3_CHANNEL
	      psi0p = _FIELD_AT(&psi0[i],ix);
	      _spinor_g3_f(sp1,*psi0p);
	      psi0p = _FIELD_AT(&psi0[2*nm+i],ix);
	      _spinor_i_plus_f(sp2,*psi0p);
	      _spinor_prod_re_f(tmp,sp1,sp2);	
	      corr[_g5g3][(zerocoord[0]+t+GLB_T-tau)%GLB_T+i*GLB_T] += tmp; 

	      psi0p = _FIELD_AT(&psi0[nm+i],ix);
	      _spinor_g3_f(sp1,*psi0p);
	      psi0p = _FIELD_AT(&psi0[3*nm+i],ix);
	      _spinor_i_minus_f(sp2,*psi0p);
	      _spinor_prod_re_f(tmp,sp1,sp2);	
	      corr[_g5g3][(zerocoord[0]+t+GLB_T-tau)%GLB_T+i*GLB_T] += tmp;

	      psi0p = _FIELD_AT(&psi0[2*nm+i],ix);
	      _spinor_g3_f(sp1,*psi0p);
	      psi0p = _FIELD_AT(&psi0[0*nm+i],ix);
	      _spinor_i_minus_f(sp2,*psi0p);
	      _spinor_prod_re_f(tmp,sp1,sp2);	
	      corr[_g5g3][(zerocoord[0]+t+GLB_T-tau)%GLB_T+i*GLB_T] += tmp; 

	      psi0p = _FIELD_AT(&psi0[3*nm+i],ix);
	      _spinor_g3_f(sp1,*psi0p);
	      psi0p = _FIELD_AT(&psi0[1*nm+i],ix);
	      _spinor_i_plus_f(sp2,*psi0p);
	      _spinor_prod_re_f(tmp,sp1,sp2);	
	      corr[_g5g3][(zerocoord[0]+t+GLB_T-tau)%GLB_T+i*GLB_T] += tmp; 
#endif

#ifdef G0G1_CHANNEL
	      psi0p = _FIELD_AT(&psi0[i],ix);
	      _spinor_g5g0g1_f(sp1,*psi0p);
	      psi0p = _FIELD_AT(&psi0[1*nm+i],ix);
	      _spinor_i_plus_f(sp2,*psi0p);
	      _spinor_prod_re_f(tmp,sp1,sp2);	
	      corr[_g0g1][(zerocoord[0]+t+GLB_T-tau)%GLB_T+i*GLB_T] += tmp; 

	      psi0p = _FIELD_AT(&psi0[nm+i],ix);
	      _spinor_g5g0g1_f(sp1,*psi0p);
	      psi0p = _FIELD_AT(&psi0[0*nm+i],ix);
	      _spinor_i_plus_f(sp2,*psi0p);
	      _spinor_prod_re_f(tmp,sp1,sp2);	
	      corr[_g0g1][(zerocoord[0]+t+GLB_T-tau)%GLB_T+i*GLB_T] += tmp;

	      psi0p = _FIELD_AT(&psi0[2*nm+i],ix);
	      _spinor_g5g0g1_f(sp1,*psi0p);
	      psi0p = _FIELD_AT(&psi0[3*nm+i],ix);
	      _spinor_i_plus_f(sp2,*psi0p);
	      _spinor_prod_re_f(tmp,sp1,sp2);	
	      corr[_g0g1][(zerocoord[0]+t+GLB_T-tau)%GLB_T+i*GLB_T] += tmp; 

	      psi0p = _FIELD_AT(&psi0[3*nm+i],ix);
	      _spinor_g5g0g1_f(sp1,*psi0p);
	      psi0p = _FIELD_AT(&psi0[2*nm+i],ix);
	      _spinor_i_plus_f(sp2,*psi0p);
	      _spinor_prod_re_f(tmp,sp1,sp2);	
	      corr[_g0g1][(zerocoord[0]+t+GLB_T-tau)%GLB_T+i*GLB_T] += tmp; 
#endif

#ifdef G0G2_CHANNEL
	      psi0p = _FIELD_AT(&psi0[i],ix);
	      _spinor_g5g0g2_f(sp1,*psi0p);
	      psi0p = _FIELD_AT(&psi0[1*nm+i],ix);
	      _spinor_minus_f(sp2,*psi0p);
	      _spinor_prod_re_f(tmp,sp1,sp2);	
	      corr[_g0g2][(zerocoord[0]+t+GLB_T-tau)%GLB_T+i*GLB_T] += tmp; 

	      psi0p = _FIELD_AT(&psi0[nm+i],ix);
	      _spinor_g5g0g2_f(sp1,*psi0p);
	      psi0p = _FIELD_AT(&psi0[0*nm+i],ix);
	      sp2 = *psi0p;
	      _spinor_prod_re_f(tmp,sp1,sp2);	
	      corr[_g0g2][(zerocoord[0]+t+GLB_T-tau)%GLB_T+i*GLB_T] += tmp;

	      psi0p = _FIELD_AT(&psi0[2*nm+i],ix);
	      _spinor_g5g0g2_f(sp1,*psi0p);
	      psi0p = _FIELD_AT(&psi0[3*nm+i],ix);
	      _spinor_minus_f(sp2,*psi0p);
	      _spinor_prod_re_f(tmp,sp1,sp2);	
	      corr[_g0g2][(zerocoord[0]+t+GLB_T-tau)%GLB_T+i*GLB_T] += tmp; 

	      psi0p = _FIELD_AT(&psi0[3*nm+i],ix);
	      _spinor_g5g0g2_f(sp1,*psi0p);
	      psi0p = _FIELD_AT(&psi0[2*nm+i],ix);
	      sp2 = *psi0p;
	      _spinor_prod_re_f(tmp,sp1,sp2);	
	      corr[_g0g2][(zerocoord[0]+t+GLB_T-tau)%GLB_T+i*GLB_T] += tmp; 	      
#endif
#ifdef G0G3_CHANNEL
	      psi0p = _FIELD_AT(&psi0[i],ix);
	      _spinor_g5g0g3_f(sp1,*psi0p);
	      psi0p = _FIELD_AT(&psi0[0*nm+i],ix);
	      _spinor_i_plus_f(sp2,*psi0p);
	      _spinor_prod_re_f(tmp,sp1,sp2);	
	      corr[_g0g3][(zerocoord[0]+t+GLB_T-tau)%GLB_T+i*GLB_T] += tmp; 

	      psi0p = _FIELD_AT(&psi0[nm+i],ix);
	      _spinor_g5g0g3_f(sp1,*psi0p);
	      psi0p = _FIELD_AT(&psi0[1*nm+i],ix);
	      _spinor_i_minus_f(sp2,*psi0p);
	      _spinor_prod_re_f(tmp,sp1,sp2);	
	      corr[_g0g3][(zerocoord[0]+t+GLB_T-tau)%GLB_T+i*GLB_T] += tmp;

	      psi0p = _FIELD_AT(&psi0[2*nm+i],ix);
	      _spinor_g5g0g3_f(sp1,*psi0p);
	      psi0p = _FIELD_AT(&psi0[2*nm+i],ix);
	      _spinor_i_plus_f(sp2,*psi0p);
	      _spinor_prod_re_f(tmp,sp1,sp2);	
	      corr[_g0g3][(zerocoord[0]+t+GLB_T-tau)%GLB_T+i*GLB_T] += tmp; 

	      psi0p = _FIELD_AT(&psi0[3*nm+i],ix);
	      _spinor_g5g0g3_f(sp1,*psi0p);
	      psi0p = _FIELD_AT(&psi0[3*nm+i],ix);
	      _spinor_i_minus_f(sp2,*psi0p);
	      _spinor_prod_re_f(tmp,sp1,sp2);	
	      corr[_g0g3][(zerocoord[0]+t+GLB_T-tau)%GLB_T+i*GLB_T] += tmp; 	      
#endif
#ifdef G0G5G1_CHANNEL
	      psi0p = _FIELD_AT(&psi0[i],ix);
	      _spinor_g0g1_f(sp1,*psi0p);
	      psi0p = _FIELD_AT(&psi0[1*nm+i],ix);
	      _spinor_i_minus_f(sp2,*psi0p);
	      _spinor_prod_re_f(tmp,sp1,sp2);	
	      corr[_g0g5g1][(zerocoord[0]+t+GLB_T-tau)%GLB_T+i*GLB_T] += tmp; 

	      psi0p = _FIELD_AT(&psi0[nm+i],ix);
	      _spinor_g0g1_f(sp1,*psi0p);
	      psi0p = _FIELD_AT(&psi0[0*nm+i],ix);
	      _spinor_i_minus_f(sp2,*psi0p);
	      _spinor_prod_re_f(tmp,sp1,sp2);	
	      corr[_g0g5g1][(zerocoord[0]+t+GLB_T-tau)%GLB_T+i*GLB_T] += tmp;

	      psi0p = _FIELD_AT(&psi0[2*nm+i],ix);
	      _spinor_g0g1_f(sp1,*psi0p);
	      psi0p = _FIELD_AT(&psi0[3*nm+i],ix);
	      _spinor_i_plus_f(sp2,*psi0p);
	      _spinor_prod_re_f(tmp,sp1,sp2);	
	      corr[_g0g5g1][(zerocoord[0]+t+GLB_T-tau)%GLB_T+i*GLB_T] += tmp; 

	      psi0p = _FIELD_AT(&psi0[3*nm+i],ix);
	      _spinor_g0g1_f(sp1,*psi0p);
	      psi0p = _FIELD_AT(&psi0[2*nm+i],ix);
	      _spinor_i_plus_f(sp2,*psi0p);
	      _spinor_prod_re_f(tmp,sp1,sp2);	
	      corr[_g0g5g1][(zerocoord[0]+t+GLB_T-tau)%GLB_T+i*GLB_T] += tmp; 
#endif
#ifdef G0G5G2_CHANNEL
	      psi0p = _FIELD_AT(&psi0[i],ix);
	      _spinor_g0g2_f(sp1,*psi0p);
	      psi0p = _FIELD_AT(&psi0[1*nm+i],ix);
	      sp2 = *psi0p;
	      _spinor_prod_re_f(tmp,sp1,sp2);	
	      corr[_g0g5g2][(zerocoord[0]+t+GLB_T-tau)%GLB_T+i*GLB_T] += tmp; 

	      psi0p = _FIELD_AT(&psi0[nm+i],ix);
	      _spinor_g0g2_f(sp1,*psi0p);
	      psi0p = _FIELD_AT(&psi0[0*nm+i],ix);
	      _spinor_minus_f(sp2,*psi0p);
	      _spinor_prod_re_f(tmp,sp1,sp2);	
	      corr[_g0g5g2][(zerocoord[0]+t+GLB_T-tau)%GLB_T+i*GLB_T] += tmp;

	      psi0p = _FIELD_AT(&psi0[2*nm+i],ix);
	      _spinor_g0g2_f(sp1,*psi0p);
	      psi0p = _FIELD_AT(&psi0[3*nm+i],ix);
	      _spinor_minus_f(sp2,*psi0p);
	      _spinor_prod_re_f(tmp,sp1,sp2);	
	      corr[_g0g5g2][(zerocoord[0]+t+GLB_T-tau)%GLB_T+i*GLB_T] += tmp; 

	      psi0p = _FIELD_AT(&psi0[3*nm+i],ix);
	      _spinor_g0g2_f(sp1,*psi0p);
	      psi0p = _FIELD_AT(&psi0[2*nm+i],ix);
	      sp2 = *psi0p;
	      _spinor_prod_re_f(tmp,sp1,sp2);	
	      corr[_g0g5g2][(zerocoord[0]+t+GLB_T-tau)%GLB_T+i*GLB_T] += tmp; 
#endif
  
#ifdef G0G5G3_CHANNEL
	      psi0p = _FIELD_AT(&psi0[i],ix);
	      _spinor_g0g3_f(sp1,*psi0p);
	      psi0p = _FIELD_AT(&psi0[0*nm+i],ix);
	      _spinor_i_minus_f(sp2,*psi0p);
	      _spinor_prod_re_f(tmp,sp1,sp2);	
	      corr[_g0g5g3][(zerocoord[0]+t+GLB_T-tau)%GLB_T+i*GLB_T] += tmp; 

	      psi0p = _FIELD_AT(&psi0[nm+i],ix);
	      _spinor_g0g3_f(sp1,*psi0p);
	      psi0p = _FIELD_AT(&psi0[1*nm+i],ix);
	      _spinor_i_plus_f(sp2,*psi0p);
	      _spinor_prod_re_f(tmp,sp1,sp2);	
	      corr[_g0g5g3][(zerocoord[0]+t+GLB_T-tau)%GLB_T+i*GLB_T] += tmp;

	      psi0p = _FIELD_AT(&psi0[2*nm+i],ix);
	      _spinor_g0g3_f(sp1,*psi0p);
	      psi0p = _FIELD_AT(&psi0[2*nm+i],ix);
	      _spinor_i_plus_f(sp2,*psi0p);
	      _spinor_prod_re_f(tmp,sp1,sp2);	
	      corr[_g0g5g3][(zerocoord[0]+t+GLB_T-tau)%GLB_T+i*GLB_T] += tmp; 

	      psi0p = _FIELD_AT(&psi0[3*nm+i],ix);
	      _spinor_g0g3_f(sp1,*psi0p);
	      psi0p = _FIELD_AT(&psi0[3*nm+i],ix);
	      _spinor_i_minus_f(sp2,*psi0p);
	      _spinor_prod_re_f(tmp,sp1,sp2);	
	      corr[_g0g5g3][(zerocoord[0]+t+GLB_T-tau)%GLB_T+i*GLB_T] += tmp; 	      

#endif

#ifdef PCAC_CHANNEL
	      for (beta = 0;beta<4;beta++){
		psi0p = _FIELD_AT(&psi0[beta*nm+i],ix);
		_spinor_g0_f(sp1,*psi0p);
		_spinor_minus_f(sp2,*psi0p);
		_spinor_prod_re_f(tmp,sp1,sp2);	
		corr[_g5_g0g5_re][(zerocoord[0]+t+GLB_T-tau)%GLB_T+i*GLB_T] += tmp; 
	      }
#endif
	      if (t+zerocoord[0]==tau){
#ifdef ID_DISC_CHANNEL
		for (beta=0;beta<4;beta++){
		  _spinor_prod_re_f(tmp,*_FIELD_AT(&psi0[beta*nm+i],ix),*_FIELD_AT(&eta[beta],ix));
		  corr[_id_disc][tau+i*GLB_T]+=tmp;
		}
#endif

#ifdef G5_DISC_CHANNEL
		for (beta=0;beta<4;beta++){
		  _spinor_g5_prod_re_f(tmp,*_FIELD_AT(&psi0[beta*nm+i],ix),*_FIELD_AT(&eta[beta],ix));
		  corr[_g5_disc][tau+i*GLB_T]+=tmp;
		}
#endif
	      }
	    }
      }
    }
  }

  
  for(k=0; k<NCHANNELS; k++) {
    global_sum(corr[k],GLB_T*nm);
    for(i=0; i<nm*GLB_T; i++)
      corr[k][i] *= -((2./nhits)/GLB_VOL3)/GLB_VOL3;
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
  /*  for(i=0; i<nm*GLB_T; i++)
      corr[_g1][i] *= -1.;*/
  PRINT_CORR(g1);
#endif

#ifdef G2_CHANNEL
  /*  for(i=0; i<nm*GLB_T; i++)
      corr[_g2][i] *= -1.;*/
  PRINT_CORR(g2);
#endif
  
  
#ifdef G3_CHANNEL
  /*  for(i=0; i<nm*GLB_T; i++)
    corr[_g3][i] *= -1.;*/
  PRINT_CORR(g3);
#endif

#ifdef G0G5_CHANNEL
  /*  for(i=0; i<nm*GLB_T; i++)
      corr[_g0g5][i] *= -1.;*/
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
  /*  for(i=0; i<nm*GLB_T; i++)
      corr[_g0g1][i] *= -1.;*/
  PRINT_CORR(g0g1);
#endif

#ifdef G0G2_CHANNEL
  /*  for(i=0; i<nm*GLB_T; i++)
      corr[_g0g2][i] *= -1.;*/
  PRINT_CORR(g0g2);
#endif
  
#ifdef G0G3_CHANNEL
  /*  for(i=0; i<nm*GLB_T; i++)
      corr[_g0g3][i] *= -1.;*/
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

#ifdef ID_DISC_CHANNEL
  PRINT_CORR(id_disc);
#endif

#ifdef G5_DISC_CHANNEL
  PRINT_CORR(g5_disc);
#endif

}
#undef GAUSSIAN_NOISE
