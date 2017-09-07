/***************************************************************************\
* Copyright (c) 2013 Rudy Arthur, Ari Hietanen, Jarno Rantaharju            *
*                                                                           *
*                                                                           *
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
#include "spin_matrix.h"
#include "propagator.h"
#include <string.h>
#include "meson_observables.h"
#define PI 3.141592653589793238462643383279502884197

meson_observable *triplet_discon_correlators = NULL; //For four fermion interaction

//If ind2=-2 disconnected correlator
static void add_meson_observable(meson_observable** mop, gamma_ind ind1, gamma_ind ind2, char* channel_name, char* channel_type,double sign){
  meson_observable* motmp,*mo;
  if (*mop==NULL){
    *mop = (meson_observable*) (malloc(sizeof(meson_observable)));
    mo=*mop;
    motmp = mo;
  }
  else{
    mo = *mop;
    motmp=mo;
    while (motmp->next!=NULL) motmp=motmp->next;
    motmp->next = (meson_observable*) (malloc(sizeof(meson_observable)));
    motmp = motmp->next;
  }
  motmp->ind1 = ind1;
  motmp->ind2 = ind2;
  strcpy(motmp->channel_name,channel_name);
  strcpy(motmp->channel_type,channel_type);
  motmp->sign = sign;
  motmp->corr_size = -1;
  motmp->corr_re = NULL;
  motmp->corr_im = NULL;
  motmp->next = NULL;
}

void init_triplet_discon_correlators(){
  int i;
  for (i=0;i<NGAMMA_IND;++i){
    char name[100];
    sprintf(name,"%s_disc",meson_channel_names[i]);
    add_meson_observable(&triplet_discon_correlators,i,_NOGAMMA,name,"TRIPLET",-1);
  }
}

static void free_observables_core(  meson_observable* mo){
  meson_observable* motmp;
  while (mo!=NULL){
    motmp=mo;
    mo=mo->next;
    free(motmp->corr_re);
    free(motmp->corr_im);
    free(motmp);
  }
}

void free_triplet_discon_observables(){
  free_observables_core(triplet_discon_correlators);
  discon_correlators=NULL;
}

static void op_spinor(suNf_spinor* out, suNf_spinor* in, gamma_ind i){
  switch (i){
  case _g5: _spinor_g5_f( (*out),(*in) ); break;
  case _id: *out=*in; break;
  case _g0: _spinor_g0_f( (*out),(*in) ); break;
  case _g1: _spinor_g1_f( (*out),(*in) ); break;
  case _g2: _spinor_g2_f( (*out),(*in) ); break;
  case _g3: _spinor_g3_f( (*out),(*in) ); break;
  case _g0g1: _spinor_g0g1_f(*out,*in); break;
  case _g0g2: _spinor_g0g2_f(*out,*in); break;
  case _g0g3: _spinor_g0g3_f(*out,*in); break;
  case _g0g5: _spinor_g0g5_f(*out,*in); break;
  case _g5g1: _spinor_g5g1_f(*out,*in); break;
  case _g5g2: _spinor_g5g2_f(*out,*in); break;
  case _g5g3: _spinor_g5g3_f(*out,*in); break;
  case _g0g5g1: _spinor_g5g0g1_f(*out,*in); break;
  case _g0g5g2: _spinor_g5g0g2_f(*out,*in); break;
  case _g0g5g3: _spinor_g5g0g3_f(*out,*in); break;
  default: break;
  }
}



#define corr_ind(px,py,pz,n_mom,tc,nm,cm) ((px)*(n_mom)*(n_mom)*(lt)*(nm)+(py)*(n_mom)*(lt)*(nm)+(pz)*(lt)*(nm)+ ((cm)*(lt)) +(tc))
static void measure_mesons_disconnected_core(spinor_field* psi0, spinor_field* psi1, spinor_field* eta, meson_observable* mo, int offset,int lt){
  int ix,t,x,y,z,tc;
  complex tr;
  //suNf_spin_matrix sma,smb, sm1,sm2,smtmp1,smtmp2,sm_src;
  meson_observable* motmp=mo;
  suNf_spinor sma, smtmp1,sm_src;
  _spinor_zero_f(smtmp1);
  lprintf("measure_mesons_disconnected_core",50,"Measuring channels: ");
  while (motmp!=NULL){
    lprintf("measure_mesons_disconnected_core",50," %s",motmp->channel_name);
    motmp=motmp->next;
  }
  lprintf("measure_mesons_disconnected_core",50,"\n");

	  for (t=0; t<T; t++) {
	    tc = (zerocoord[0]+t+GLB_T )%GLB_T+offset;
	    for (x=0; x<X; x++) for (y=0; y<Y; y++) for (z=0; z<Z; z++) { 
		  ix=ipt(t,x,y,z);
		  sma=*_FIELD_AT(psi0,ix);
        sm_src=*_FIELD_AT(eta,ix);
		  motmp=mo;
		  while (motmp!=NULL){
		    op_spinor(&smtmp1,&sma,motmp->ind1);
		    _spinor_prod_f(tr,smtmp1,sm_src);
		    motmp->corr_re[corr_ind(0,0,0,1,tc,1,0)] += motmp->sign*(tr.re);
		    motmp->corr_im[corr_ind(0,0,0,1,tc,1,0)] += motmp->sign*(tr.im);
		    motmp=motmp->next;
		  }
            }
	  }
      
  lprintf("measure_mesons_core",50,"Measuring DONE! ");
}

static void init_corrs(int nm, int n_mom, meson_observable* mo){
  int i,n_mom_tot,size;
  meson_observable* motmp=mo;
  if (n_mom>=GLB_X){
    n_mom=GLB_X-1;
    lprintf("measure_mesons_with_momenta",0,"Reduced n_mom to %d (no more momenta accessible with this lattice size) ",n_mom);
  }
  n_mom_tot = n_mom*n_mom*n_mom;
  size = GLB_T*nm*n_mom_tot;
  while (motmp!=NULL){
    if (size!=motmp->corr_size){
      free(motmp->corr_re);
      free(motmp->corr_im);
      motmp->corr_re=malloc(sizeof(double)*size);
      motmp->corr_im=malloc(sizeof(double)*size);    
      motmp->corr_size=size;
      for(i=0; i<motmp->corr_size; i++){
	motmp->corr_re[i] = 0.;
	motmp->corr_im[i] = 0.;
      }
    }
    motmp=motmp->next;
  }
}


void measure_mesons_disconnected(meson_observable* mo,spinor_field* psi0, spinor_field* eta ){
  init_corrs(1,1,mo);
  lprintf("measure_mesons",50,"measure default mesons");
  measure_mesons_disconnected_core(psi0, psi0, eta, mo, 0,GLB_T);
}


