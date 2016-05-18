/***************************************************************************\
* Copyright (c) 2013 Rudy Arthur, Ari Hietanen                              *
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

meson_observable *meson_correlators  = NULL;
meson_observable *discon_correlators = NULL;
meson_observable *cvc_correlators = NULL;

char* meson_channel_names[NGAMMA_IND]={"g5","id","g0","g1","g2","g3","g0g1","g0g2","g0g3","g0g5","g5g1","g5g2","g5g3","g0g5g1","g0g5g2","g0g5g3"};

static int vector_gammas[4] = {_g0, _g1, _g2, _g3};

static double determine_sign(gamma_ind ind){
  if (ind == _g5 || ind == _g0 || ind == _g1 || ind == _g3 || ind == _g0g1 || ind == _g0g3 || ind == _g0g5 || ind == _g5g1 || ind==_g5g3 || ind==_g0g5g2){
    return -1.;
  }
  else {
    return 1.;
  }
}

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

void init_meson_correlators(int meas_offdiag){
  int i,j;
  char name[100];
  for (i=0;i<NGAMMA_IND;++i){
    add_meson_observable(&meson_correlators,i,i,meson_channel_names[i],"TRIPLET",determine_sign(i));
  }
  add_meson_observable(&meson_correlators,_g5,_g0g5,"g5_g0g5","TRIPLET",determine_sign(_g0g5));
  if (meas_offdiag){
    for (i=0;i<NGAMMA_IND;++i){
      for (j=0;j<NGAMMA_IND;++j){
	if ( i!=j && !(i==_g5 && j == _g0g5) ){
	  sprintf(name,"%s_%s",meson_channel_names[i],meson_channel_names[j]);
	  add_meson_observable(&meson_correlators,i,j,name,"TRIPLET",determine_sign(j));
	}
      }
    }  
  }
}


void init_discon_correlators(){
  int i;
  for (i=0;i<NGAMMA_IND;++i){
    char name[100];
    sprintf(name,"%s_disc",meson_channel_names[i]);
    add_meson_observable(&discon_correlators,i,_NOGAMMA,name,"SINGLET",-1);
  }
}


void init_cvc_correlators(){
  int i;
  for (i=0;i<4;++i){
    char name[100];
    sprintf(name,"%s_CL",meson_channel_names[ vector_gammas[i] ]);
    add_meson_observable(&cvc_correlators,vector_gammas[i],i,name,"TRIPLET", determine_sign( vector_gammas[i] ));
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

void free_meson_observables(){
  free_observables_core(meson_correlators);
  meson_correlators=NULL;
  free_observables_core(discon_correlators);
  discon_correlators=NULL;
}


static void spinmatrix_op(suNf_spin_matrix *out, suNf_spin_matrix *in, gamma_ind i){
  switch (i){
  case _g5: _spinmatrix_g5(*out,*in); break;
  case _id: *out=*in; break;
  case _g0: _spinmatrix_g0(*out,*in); break;
  case _g1: _spinmatrix_g1(*out,*in); break;
  case _g2: _spinmatrix_g2(*out,*in); break;
  case _g3: _spinmatrix_g3(*out,*in); break;
  case _g0g1: _spinmatrix_g0g1(*out,*in); break;
  case _g0g2: _spinmatrix_g0g2(*out,*in); break;
  case _g0g3: _spinmatrix_g0g3(*out,*in); break;
  case _g0g5: _spinmatrix_g0g5(*out,*in); break;
  case _g5g1: _spinmatrix_g5g1(*out,*in); break;
  case _g5g2: _spinmatrix_g5g2(*out,*in); break;
  case _g5g3: _spinmatrix_g5g3(*out,*in); break;
  case _g0g5g1: _spinmatrix_g5g0g1(*out,*in); break;
  case _g0g5g2: _spinmatrix_g5g0g2(*out,*in); break;
  case _g0g5g3: _spinmatrix_g5g0g3(*out,*in); break;
  default: break;
  }
}

static void op_spinmatrix(suNf_spin_matrix* out, suNf_spin_matrix* in, gamma_ind i){
  switch (i){
  case _g5: _g5_spinmatrix(*out,*in); break;
  case _id: *out=*in; break;
  case _g0: _g0_spinmatrix(*out,*in); break;
  case _g1: _g1_spinmatrix(*out,*in); break;
  case _g2: _g2_spinmatrix(*out,*in); break;
  case _g3: _g3_spinmatrix(*out,*in); break;
  case _g0g1: _g0g1_spinmatrix(*out,*in); break;
  case _g0g2: _g0g2_spinmatrix(*out,*in); break;
  case _g0g3: _g0g3_spinmatrix(*out,*in); break;
  case _g0g5: _g0g5_spinmatrix(*out,*in); break;
  case _g5g1: _g5g1_spinmatrix(*out,*in); break;
  case _g5g2: _g5g2_spinmatrix(*out,*in); break;
  case _g5g3: _g5g3_spinmatrix(*out,*in); break;
  case _g0g5g1: _g5g0g1_spinmatrix(*out,*in); break;
  case _g0g5g2: _g5g0g2_spinmatrix(*out,*in); break;
  case _g0g5g3: _g5g0g3_spinmatrix(*out,*in); break;
  default: break;
  }
}




static void op_propagator(suNf_propagator* out, suNf_propagator* in, gamma_ind i){
  switch (i){
  case _g5: _g5_propagator( (*out),(*in) ); break;
  case _id: *out=*in; break;
  case _g0: _g0_propagator( (*out),(*in) ); break;
  case _g1: _g1_propagator( (*out),(*in) ); break;
  case _g2: _g2_propagator( (*out),(*in) ); break;
  case _g3: _g3_propagator( (*out),(*in) ); break;
  default: break;
  }
}




#define corr_ind(px,py,pz,n_mom,tc,nm,cm) ((px)*(n_mom)*(n_mom)*(lt)*(nm)+(py)*(n_mom)*(lt)*(nm)+(pz)*(lt)*(nm)+ ((cm)*(lt)) +(tc))

/* spinor_fields* are 4xnm arrays of spinor_field ordered([color][spinor])*/
void measure_mesons_core(spinor_field* psi0, spinor_field* psi1, spinor_field* eta, meson_observable* mo, int nm, int tau, int n_mom, int offset,int lt){
  int i,ix,t,x,y,z,beta,px,py,pz,tc;
  double pdotx,cpdotx,spdotx;
  complex tr;
  suNf_spin_matrix sma,smb, sm1,sm2,smtmp1,smtmp2,sm_src;
  meson_observable* motmp=mo;
  _spinmatrix_zero(smtmp1);
  _spinmatrix_zero(smtmp2);
  lprintf("measure_mesons_core",50,"Measuring channels: ");
  while (motmp!=NULL){
    lprintf("",50," %s",motmp->channel_name);
    motmp=motmp->next;
  }
  lprintf("",50,"\n");
  for (px=0;px<n_mom;++px) for (py=0;py<n_mom;++py) for (pz=0;pz<n_mom;++pz){
	for(i=0; i<nm; i++) {
	  for (t=0; t<T; t++) {	 
	    tc = (zerocoord[0]+t+GLB_T-tau)%GLB_T+offset;
	    for (x=0; x<X; x++) for (y=0; y<Y; y++) for (z=0; z<Z; z++) { 
		  ix=ipt(t,x,y,z);					
		  pdotx = 2.0*PI*(((double) px)*(x+zerocoord[1])/GLB_X + ((double) py)*(y+zerocoord[2])/GLB_Y + ((double) pz)*(z+zerocoord[3])/GLB_Z);
		  cpdotx=cos(pdotx);
		  spdotx=sin(pdotx);
		  for (beta=0;beta<4;beta++){ 
		    _spinmatrix_assign_row(sma, *_FIELD_AT(&psi0[beta*nm+i],ix), beta);
		    _spinmatrix_assign_row(smb, *_FIELD_AT(&psi1[beta*nm+i],ix), beta);
		    _spinmatrix_assign_row(sm_src, *_FIELD_AT(&eta[beta],ix), beta);
		  }
		  motmp=mo;
		  while (motmp!=NULL){
		    if (motmp->ind2!=_NOGAMMA){
		      spinmatrix_op(&smtmp1,&sma,motmp->ind2);
		      _spinmatrix_g5(sm1,smtmp1);
		      op_spinmatrix(&smtmp2,&smb,motmp->ind1);
		      _g5_spinmatrix(sm2,smtmp2);
		      _spinmatrix_mul_trace(tr,sm1,sm2);
		    }
		    else{
		      spinmatrix_op(&smtmp1,&sma,motmp->ind1);
		      _spinmatrix_mul_trace(tr,smtmp1,sm_src);
		    }
		    motmp->corr_re[corr_ind(px,py,pz,n_mom,tc,nm,i)] += motmp->sign*(tr.re*cpdotx+tr.im*spdotx);
		    motmp->corr_im[corr_ind(px,py,pz,n_mom,tc,nm,i)] += motmp->sign*(tr.im*cpdotx-tr.re*spdotx);
		    motmp=motmp->next;
		  }
		}
	  }
	}
      }
  lprintf("measure_mesons_core",50,"Measuring DONE! ");
}



static void measure_conserved_core(spinor_field* psi0, spinor_field* psi1, spinor_field* eta, meson_observable* mo, int nm, int tau, int n_mom, int offset,int lt){

  int i,ix,t,x,y,z,beta,px,py,pz,tc,a,ixmu;
  double pdotx,cpdotx,spdotx;
  complex tr;
  suNf_propagator sp0,sp1,Usp,spf,sptmp1,sptmp2,sptmp3,spleft,spdag;
  suNf *u1;
  meson_observable* motmp=mo;

  lprintf("measure_conserved_core",50,"Measuring channels: ");
  while (motmp!=NULL){
    lprintf("",50," %s",motmp->channel_name);
    motmp=motmp->next;
  }
  lprintf("",50,"\n");
  for (px=0;px<n_mom;++px) for (py=0;py<n_mom;++py) for (pz=0;pz<n_mom;++pz){
	for(i=0; i<nm; i++) {
	  for (t=0; t<T; t++) {	 
	    tc = (zerocoord[0]+t+GLB_T-tau)%GLB_T+i*GLB_T+offset;
	    for (x=0; x<X; x++) for (y=0; y<Y; y++) for (z=0; z<Z; z++) {
		  ix=ipt(t,x,y,z);

					
		  pdotx = 2.0*PI*(((double) px)*(x+zerocoord[1])/GLB_X + ((double) py)*(y+zerocoord[2])/GLB_Y + ((double) pz)*(z+zerocoord[3])/GLB_Z);
		  cpdotx=cos(pdotx);
		  spdotx=sin(pdotx);
		  for (a=0;a<NF;++a){
		    for (beta=0;beta<4;beta++){ 
		      _propagator_assign(sp0, *_FIELD_AT(&psi0[a*4*nm+beta*nm+i],ix),a,beta);
		      _propagator_assign(sp1, *_FIELD_AT(&psi0[a*4*nm+beta*nm+i],ix),a,beta);
		    }
		  }
                  _propagator_dagger(spdag,sp1);

		  motmp=mo;
		  while (motmp!=NULL){

			  ixmu = iup(ix,motmp->ind2);
			  for (a=0;a<NF;++a){
			    for (beta=0;beta<4;beta++){ 
			      _propagator_assign(spf, *_FIELD_AT(&psi0[a*4*nm+beta*nm+i],ixmu), a,beta); //S(x+mu, y)
			    }
			  }
			  u1 = _4FIELD_AT(u_gauge_f,ix,motmp->ind2);

			  // Tr [ g5 (1+g_mu) U^(x) S(x,0) g5 g_mu S^(x+mu, y) ]
			  _suNf_inverse_prop_multiply(Usp,*u1,sp0); //U^(x) S(x,0) 
			  sptmp1=Usp;
			  op_propagator(&sptmp2,&sptmp1,motmp->ind1); //g_mu U^(x) S(x,0) 
			  _propagator_add(sptmp1,sptmp1,sptmp2); //(1+g_mu) U^(x) S(x,0) 
			  _g5_propagator(sptmp2,sptmp1); //g5 (1+g_mu) U^(x) S(x,0) 
			  _propagator_dagger(sptmp1,spf); //S^(x+mu, 0)
			  _g5_propagator(sptmp3,sptmp1); //g5 g_mu S^(x+mu, 0)
			  op_propagator(&sptmp1,&sptmp3,motmp->ind1); //g_mu S^(x+mu, 0)
			  _propagator_mul(spleft,sptmp2,sptmp1);	
			    

			  // Tr [ g5 (1+g_mu) U(x) S(x+mu,0) g5 g_mu S^(x, y) ]
			  _suNf_prop_multiply(Usp,*u1,spf); //U(x) S(x+mu,0)
			  sptmp1=Usp;
			  op_propagator(&sptmp2,&sptmp1,motmp->ind1);//g_mu U(x) S(x,0) 
			  _propagator_sub(sptmp1,sptmp1,sptmp2);//(1-g_mu) U(x) S(x,0) 
			  _g5_propagator(sptmp2,sptmp1);//g5(1-g_mu) U(x) S(x,0) 
			  _g5_propagator(sptmp1,spdag); //g5 S^(x, 0)
			  op_propagator(&sptmp3,&sptmp1,motmp->ind1); //g_mu g5 S^(x, 0)
			  _propagator_mul(sptmp1,sptmp2,sptmp3);

			  _propagator_sub(sptmp2,spleft,sptmp1);
			  _propagator_trace(tr,sptmp2);    

			  motmp->corr_re[corr_ind(px,py,pz,n_mom,tc,nm,i)] += 0.5*motmp->sign*(tr.re*cpdotx+tr.im*spdotx);
			  motmp->corr_im[corr_ind(px,py,pz,n_mom,tc,nm,i)] += 0.5*motmp->sign*(tr.im*cpdotx-tr.re*spdotx);
			  motmp=motmp->next;

		  } //END CORRELATOR LOOP
		} //END SPATIAL LOOP
	  } //END T LOOP
	} //END MASS LOOP
      } //END MOMENTUM LOOP
  lprintf("measure_formfactor_core",50,"Measuring DONE! ");
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

static void zero_corrs(meson_observable* mo){
  meson_observable* motmp=mo;
  int i;
  while (motmp!=NULL){
    for(i=0; i<motmp->corr_size; i++){
      motmp->corr_re[i] = 0.;
      motmp->corr_im[i] = 0.;
    }
    motmp=motmp->next;
  }
}
void measure_diquarks(meson_observable* mo,spinor_field* psi0, spinor_field* psi1,spinor_field* eta, int nm, int tau){
  init_corrs(nm,1,mo);
  lprintf("measure_mesons",50,"measure default diquarks");
  measure_mesons_core(psi0, psi1, eta, mo,nm, tau, 1, 0,GLB_T);
//  measure_mesons_core(psi1, psi0, eta, mo,nm, tau, 1, 0,GLB_T);
}


void measure_mesons(meson_observable* mo,spinor_field* psi0, spinor_field* eta, int nm, int tau){
  init_corrs(nm,1,mo);
  lprintf("measure_mesons",50,"measure default mesons");
  measure_mesons_core(psi0, psi0, eta, mo,nm, tau, 1, 0,GLB_T);
}

void measure_conserved_currents(meson_observable* mo,spinor_field* psi0, spinor_field* eta, int nm, int tau){
  init_corrs(nm,1,mo);
  lprintf("measure_mesons",50,"measure conserved vector currents");
  measure_conserved_core(psi0, psi0, eta, mo,nm, tau, 1, 0,GLB_T);
}

void measure_mesons_ext(meson_observable* mo,spinor_field* psi0, spinor_field* eta, int nm, int tau,int begin){
  init_corrs(nm*3,1,mo);
  lprintf("measure_mesons",50,"measure extended mesons segment %d",begin);
  measure_mesons_core(psi0, psi0,eta, mo,nm, tau, 1, GLB_T*begin, 3*GLB_T);
}


void measure_point_mesons_momenta(meson_observable* mo,spinor_field* psi0, spinor_field* eta, int nm, int tau, int n_mom){
  init_corrs(nm,n_mom,mo);
  lprintf("measure_mesons",50,"measure point mesons with momenta");
  measure_mesons_core(psi0, psi0, eta, mo, nm, tau, n_mom, 0, GLB_T);
}

void measure_point_mesons_momenta_ext(meson_observable* mo,spinor_field* psi0, spinor_field* eta, int nm, int tau, int n_mom, int begin){
  init_corrs(3*nm,n_mom,mo);
  lprintf("measure_mesons",50,"measure extended point mesons with moment");
  measure_mesons_core(psi0, psi0, eta, mo, nm, tau, n_mom, GLB_T*begin, 3*GLB_T);
}

static void print_corr_core(meson_observable* mo,int lt,int conf, int nm, double* mass, char* label, int n_mom){
  int i,t,px,py,pz;
  for (px=0;px<n_mom;++px) for (py=0;py<n_mom;++py) for (pz=0;pz<n_mom;++pz){ 
	for(i=0; i<nm; i++) {
	  /*Real*/
	  if (n_mom>1){
	    lprintf("MAIN",0,"conf #%d mass=%2.6f %s %s %s_re momentum(%d,%d,%d)= ",conf,mass[i],label,mo->channel_type,mo->channel_name,px,py,pz); 
	  }
	  else{
	    if (mo->ind1==mo->ind2){	    
	      lprintf("MAIN",0,"conf #%d mass=%2.6f %s %s %s= ",conf,mass[i],label,mo->channel_type,mo->channel_name); //To be compatible with the old output
	    }
	    else{
	      lprintf("MAIN",0,"conf #%d mass=%2.6f %s %s %s_re= ",conf,mass[i],label,mo->channel_type,mo->channel_name); 
	    }
	  }
	  for(t=0;t<lt;++t) lprintf("MAIN",0,"%e ",mo->corr_re[corr_ind(px,py,pz,n_mom,t,nm,i)]);
	  lprintf("MAIN",0,"\n");
	  /*Imaginary */
	  if (n_mom>1){
	    lprintf("MAIN",0,"conf #%d mass=%2.6f %s %s %s_im momentum(%d,%d,%d)= ",conf,mass[i],label,mo->channel_type,mo->channel_name,px,py,pz); 
	  }
	  else{
	      lprintf("MAIN",0,"conf #%d mass=%2.6f %s %s %s_im= ",conf,mass[i],label,mo->channel_type,mo->channel_name); 
	  }
	  for(t=0;t<lt;++t) lprintf("MAIN",0,"%e ",mo->corr_im[corr_ind(px,py,pz,n_mom,t,nm,i)]);
	  lprintf("MAIN",0,"\n"); 
	}
      }
}


static void print_corr(meson_observable* mo,int lt,int conf, int nm, double* mass, char* label, int n_mom){
  meson_observable* motmp=mo;
  while (motmp!=NULL){  
    print_corr_core(motmp,lt,conf,nm,mass,label,n_mom);
    motmp=motmp->next;
  }
}

static void do_global_sum(meson_observable* mo, double norm){
  meson_observable* motmp=mo;
  int i;
  while (motmp!=NULL){
      global_sum(motmp->corr_re,motmp->corr_size);
      global_sum(motmp->corr_im,motmp->corr_size);
      for(i=0; i<motmp->corr_size; i++){
	motmp->corr_re[i] *= norm;
	motmp->corr_im[i] *= norm;
      }
    motmp=motmp->next;
  }
}

void print_mesons(meson_observable* mo,double norm, int conf, int nm, double* mass, int lt, int n_mom, char* label){
  do_global_sum(mo,-((1./norm)/GLB_VOL3));
  print_corr(mo,lt,conf,nm,mass,label,n_mom);
  zero_corrs(mo);
}
