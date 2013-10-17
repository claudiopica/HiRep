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
#define PI 3.141592653589793238462643383279502884197

enum { _g5= 0,_id, _g0, _g1, _g2, _g3,  _g0g1, _g0g2, _g0g3, _g0g5, _g5g1, _g5g2, _g5g3, _g0g5g1, _g0g5g2, _g0g5g3, _g5_g0g5_re, _id_disc,_g5_disc,_pipig_ff_re,_pipig_ff_im,_pipig_conserved_ff,NCHANNELS };

char* channel_names[NCHANNELS]={"g5","id","g0","g1","g2","g3","g0g1","g0g2","g0g3","g0g5","g5g1","g5g2","g5g3","g0g5g1","g0g5g2","g0g5g3","g5_g0g5_re","id_disc", "g5_disc","pipig_ff_re","pipig_ff_im","pipig_conserved_ff"};

char* channel_types[NCHANNELS]={"TRIPLET","TRIPLET","TRIPLET","TRIPLET","TRIPLET","TRIPLET","TRIPLET","TRIPLET","TRIPLET","TRIPLET","TRIPLET","TRIPLET","TRIPLET","TRIPLET","TRIPLET","TRIPLET","TRIPLET","SINGLET","SINGLET","FORM_FACTOR","FORM_FACTOR","FORM_FACTOR"};

int measure_channels[NCHANNELS]={1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0};

static double* corr[NCHANNELS];

static int init = 0;

#define corr_ind(px,py,pz,n_mom,tc,nm) ((px)*(n_mom)*(n_mom)*(lt)*(nm)+(py)*(n_mom)*(lt)*(nm)+(pz)*(lt)*(nm) +(tc))

/* spinor_fields* are 4xnm arrays of spinor_field ordered([color][spinor])*/
static void measure_mesons_core(spinor_field* psi0, spinor_field* psi1, spinor_field* eta, int nm, int tau, int n_mom, int offset,int lt){
  int i,ix,t,x,y,z,beta,px,py,pz,tc;
  double pdotx,cpdotx,spdotx;
  complex tr;
  suNf_spin_matrix sma,smb, sm1,sm2;
  lprintf("measure_mesons_core",50,"Measuring channels: ");
  for (i=0;i<NCHANNELS-1;++i){
    if (measure_channels[i]) lprintf("",50," %s",channel_names[i]);
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
		  for (beta=0;beta<4;beta++){ 
		    _spinmatrix_assign_row(sma, *_FIELD_AT(&psi0[beta*nm+i],ix), beta);
		    _spinmatrix_assign_row(smb, *_FIELD_AT(&psi1[beta*nm+i],ix), beta); 
		  }
		  if (measure_channels[_g5]){
		    _spinmatrix_mul_trace(tr, sma, smb);
		    corr[_g5][corr_ind(px,py,pz,n_mom,tc,nm)]+= -cpdotx*tr.re;
		    //		    corr_im[_g5][corr_ind(px,py,pz,n_mom,tc,nm)]+= spdotx*tr.re;
		  }
		  if(measure_channels[_id]){
		    _g5_spinmatrix(sm1, sma);
		    _spinmatrix_g5(sm2, smb);
		    _spinmatrix_mul_trace(tr, sm1, sm2);
		    corr[_id][corr_ind(px,py,pz,n_mom,tc,nm)] += cpdotx*tr.re;
		    //		    corr_im[_id][corr_ind(px,py,pz,n_mom,tc,nm)] += -spdotx*tr.re;
		  }
		  if(measure_channels[_g0]){
		    _g5g0_spinmatrix(sm1, sma);
		    _spinmatrix_g5g0(sm2, smb);
		    _spinmatrix_mul_trace(tr, sm1, sm2);
		    corr[_g0][corr_ind(px,py,pz,n_mom,tc,nm)] += -cpdotx*tr.re;
		    //		    corr_im[_g0][corr_ind(px,py,pz,n_mom,tc,nm)] += spdotx*tr.re;
		  }
		  if(measure_channels[_g1]){
		    _g5g1_spinmatrix(sm1, sma);
		    _spinmatrix_g5g1(sm2, smb);
		    _spinmatrix_mul_trace(tr, sm1, sm2);
		    corr[_g1][corr_ind(px,py,pz,n_mom,tc,nm)] += -cpdotx*tr.re;
		    //		    corr_im[_g1][corr_ind(px,py,pz,n_mom,tc,nm)] += spdotx*tr.re;
		  }      
		  if(measure_channels[_g2]){
		    _g5g2_spinmatrix(sm1, sma);
		    _spinmatrix_g5g2(sm2, smb);
		    _spinmatrix_mul_trace(tr, sm2, sm1);
		    corr[_g2][corr_ind(px,py,pz,n_mom,tc,nm)] += cpdotx*tr.re;
		    //		    corr_im[_g2][corr_ind(px,py,pz,n_mom,tc,nm)] += -spdotx*tr.re;
		  }
		  if(measure_channels[_g3]){
		    _g5g3_spinmatrix(sm1, sma);
		    _spinmatrix_g5g3(sm2, smb);
		    _spinmatrix_mul_trace(tr, sm1, sm2);
		    corr[_g3][corr_ind(px,py,pz,n_mom,tc,nm)] += -cpdotx*tr.re;
		    //		    corr_im[_g3][corr_ind(px,py,pz,n_mom,tc,nm)] += spdotx*tr.re;
		  }
		  if(measure_channels[_g0g5]){
		    _g0_spinmatrix(sm1, sma);
		    _spinmatrix_g0(sm2, smb);
		    _spinmatrix_mul_trace(tr, sm1, sm2);
		    corr[_g0g5][corr_ind(px,py,pz,n_mom,tc,nm)] += -cpdotx*tr.re;
		    //		    corr_im[_g0g5][corr_ind(px,py,pz,n_mom,tc,nm)] += spdotx*tr.re;
		  }
		  if(measure_channels[_g5g1]){
		    _g1_spinmatrix(sm1, sma);
		    _spinmatrix_g1(sm2, smb);
		    _spinmatrix_mul_trace(tr, sm1, sm2);
		    corr[_g5g1][corr_ind(px,py,pz,n_mom,tc,nm)] += -cpdotx*tr.re;
		    //		    corr_im[_g5g1][corr_ind(px,py,pz,n_mom,tc,nm)] += spdotx*tr.re;
		  }	      
		  if(measure_channels[_g5g2]){
		    _g2_spinmatrix(sm1, sma);
		    _spinmatrix_g2(sm2, smb);
		    _spinmatrix_mul_trace(tr, sm1, sm2);
		    corr[_g5g2][corr_ind(px,py,pz,n_mom,tc,nm)] += cpdotx*tr.re;
		    //		    corr_im[_g5g2][corr_ind(px,py,pz,n_mom,tc,nm)] += -spdotx*tr.re;
		  }
		  if(measure_channels[_g5g3]){
		    _g3_spinmatrix(sm1, sma);
		    _spinmatrix_g3(sm2, smb);
		    _spinmatrix_mul_trace(tr, sm1, sm2);
		    corr[_g5g3][corr_ind(px,py,pz,n_mom,tc,nm)] += -cpdotx*tr.re;
		    //		    corr_im[_g5g3][corr_ind(px,py,pz,n_mom,tc,nm)] += spdotx*tr.re;
		  }
		  if(measure_channels[_g0g1]){
		    _g5g0g1_spinmatrix(sm1, sma);
		    _spinmatrix_g5g0g1(sm2, smb);
		    _spinmatrix_mul_trace(tr, sm1, sm2);
		    corr[_g0g1][corr_ind(px,py,pz,n_mom,tc,nm)] += -cpdotx*tr.re;
		    //		    corr_im[_g0g1][corr_ind(px,py,pz,n_mom,tc,nm)] += spdotx*tr.re;
		  }
		  if(measure_channels[_g0g2]){
		    _g5g0g2_spinmatrix(sm1, sma);
		    _spinmatrix_g5g0g2(sm2, smb);
		    _spinmatrix_mul_trace(tr, sm1, sm2);
		    corr[_g0g2][corr_ind(px,py,pz,n_mom,tc,nm)] += cpdotx*tr.re;
		    //		    corr_im[_g0g2][corr_ind(px,py,pz,n_mom,tc,nm)] += -spdotx*tr.re;
		  }
		  if(measure_channels[_g0g3]){
		    _g5g0g3_spinmatrix(sm1, sma);
		    _spinmatrix_g5g0g3(sm2, smb);
		    _spinmatrix_mul_trace(tr, sm1, sm2);
		    corr[_g0g3][corr_ind(px,py,pz,n_mom,tc,nm)] += -cpdotx*tr.re;
		    //		    corr_im[_g0g3][corr_ind(px,py,pz,n_mom,tc,nm)] += spdotx*tr.re;
		  }
		  if(measure_channels[_g0g5g1]){
		    _g0g1_spinmatrix(sm1, sma);
		    _spinmatrix_g0g1(sm2, smb);
		    _spinmatrix_mul_trace(tr, sm1, sm2);
		    corr[_g0g5g1][corr_ind(px,py,pz,n_mom,tc,nm)] += +cpdotx*tr.re;
		    //		    corr_im[_g0g5g1][corr_ind(px,py,pz,n_mom,tc,nm)] += -spdotx*tr.re;
		  }
		  if(measure_channels[_g0g5g2]){
		    _g0g2_spinmatrix(sm1, sma);
		    _spinmatrix_g0g2(sm2, smb);
		    _spinmatrix_mul_trace(tr, sm1, sm2);
		    corr[_g0g5g2][corr_ind(px,py,pz,n_mom,tc,nm)] += -cpdotx*tr.re;
		    //		    corr_im[_g0g5g2][corr_ind(px,py,pz,n_mom,tc,nm)] += +spdotx*tr.re;
		  }
		  if(measure_channels[_g0g5g3]){
		    _g0g3_spinmatrix(sm1, sma);
		    _spinmatrix_g0g3(sm2, smb);
		    _spinmatrix_mul_trace(tr, sm1, sm2);
		    corr[_g0g5g3][corr_ind(px,py,pz,n_mom,tc,nm)] += +cpdotx*tr.re;
		    //		    corr_im[_g0g5g3][corr_ind(px,py,pz,n_mom,tc,nm)] += -spdotx*tr.re;
		  }
		  if(measure_channels[_g5_g0g5_re]){
		    _spinmatrix_g0(sm2, smb);
		    _spinmatrix_mul_trace(tr, sma, sm2);
		    corr[_g5_g0g5_re][corr_ind(px,py,pz,n_mom,tc,nm)] += cpdotx*tr.re+spdotx*tr.im;
		    //		    corr_im[_g5_g0g5_re][corr_ind(px,py,pz,n_mom,tc,nm)] += -spdotx*tr.re;
		  }
		  if(measure_channels[_id_disc]){
		    for (beta=0;beta<4;beta++){
		      _spinor_prod_re_f(tr.re,*_FIELD_AT(&psi0[beta*nm+i],ix),*_FIELD_AT(&eta[beta],ix));
		      corr[_id_disc][tau+i*lt+offset] += -tr.re;
		    }
		  }
		  if(measure_channels[_g5_disc]){
		    for (beta=0;beta<4;beta++){
		      _spinor_g5_prod_re_f(tr.re,*_FIELD_AT(&psi0[beta*nm+i],ix),*_FIELD_AT(&eta[beta],ix));
		      corr[_g5_disc][tau+i*lt+offset]+=tr.re;
		    }
		  }
		  if(measure_channels[_pipig_ff_re]){
		    _g5g0_spinmatrix(sm1, sma); //g5 g0 psi1
		    _spinmatrix_mul_trace(tr, sm1, smb); // g5 g0 psi1 psi0^dagger
		    corr[_pipig_ff_re][corr_ind(px,py,pz,n_mom,tc,nm)] += +cpdotx*tr.re+spdotx*tr.im;;
		    corr[_pipig_ff_im][corr_ind(px,py,pz,n_mom,tc,nm)] += -spdotx*tr.re+cpdotx*tr.im;
		  }
		}
	  }
	}
      }
  lprintf("measure_mesons_core",50,"Measuring DONE! ");
}

/* spinor_fields* are nmx4xNF arrays of spinor_field ordered([nm][color][spinor])*/
static void measure_nonlocal_core(spinor_field* psi0, spinor_field* psi1, spinor_field* eta, int nm, int tau, int n_mom, int offset,int lt){
  int i,ix,t,x,y,z,beta,px,py,pz,tc,a;
  double pdotx,cpdotx,spdotx;
  complex tr;
  suNf_propagator sp0,sp1,Usp,spf,sptmp1,sptmp2;
  lprintf("measure_mesons_core",50,"Measuring channels: ");
  for (i=NCHANNELS-1;i<NCHANNELS;++i){
    if (measure_channels[i]) lprintf("",50," %s",channel_names[i]);
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
		      _propagator_assign(sp1, *_FIELD_AT(&psi1[a*4*nm+beta*nm+i],ix),a,beta);
		    }
		  }
		  if (measure_channels[_pipig_conserved_ff]){
		    suNf *u1 = _4FIELD_AT(u_gauge_f,ix,0);
		    int ixmu = iup(ix,0);
		    int alpha,a;
		    /* Tr [ 1/2 S^(0,y+\mu) g_5(1+g_0) U^(y) S(y,x) */
		    _suNf_inverse_prop_multiply(Usp,*u1,sp0);
		    sptmp1=Usp;
		    _g0_propagator(stmp2,stmp1);
		    _propagator_add(stmp1,stmp1,stmp2);
		    _g5_propagator(stmp2,stmp1);
		    for (a=0;a<NF;a++) for (beta=0;beta<4;beta++){ 
			_propagator_assign(spf, *_FIELD_AT(&psi1[4*a*nm+beta*nm+i],ixmu), a,beta); 
		    }
		    _tr_propagator_mul_propagator(tr,stmp2,sfp);
						  
		    corr[_pipig_conserved_ff][corr_ind(px,py,pz,n_mom,tc,nm)] += +0.5*(cpdotx*tr.re+spdotx*tr.im);
		    /* Tr [ -1/2 S^(0,y) g_5(1-g_0) U(y) S(y+mu,x) */
		    for (a=0;a<NF;a++) for (beta=0;beta<4;beta++){ 
			_propagator_assign(spf, *_FIELD_AT(&psi0[4*a*nm+beta*nm+i],ixmu), a,beta); 
		      }
		    _suNf_prop_multiply(Usp,*u1,spf);
		    sptmp1=Usp;
		    _g0_propagator(stmp2,stmp1);
		    _propagator_sub(stmp1,stmp1,stmp2);
		    _g5_propagator(stmp2,stmp1);
		    _tr_propagator_mul_propagator(tr,stmp2,sp1);		    
		    corr[_pipig_conserved_ff][corr_ind(px,py,pz,n_mom,tc,nm)] += -0.5*(cpdotx*tr.re+spdotx*tr.im);
		  }
		}
	  }
	}
      }
  lprintf("measure_nonlocal_core",50,"Measuring DONE! ");
}		  


static void init_corrs(int nm, int n_mom){
  int i,k,n_mom_tot,size;
  static int size_old=-1;
  if (n_mom>=GLB_X){
    n_mom=GLB_X-1;
    lprintf("measure_mesons_with_momenta",0,"Reduced n_mom to %d (no more momenta accessible with this lattice size) ",n_mom);
  }
  n_mom_tot = n_mom*n_mom*n_mom;
  size = GLB_T*nm*n_mom_tot;
  if (size_old!=size){
    if (size_old!=-1){
      for (i=0;i<NCHANNELS;++i){
	free(corr[i]);
      }
    }
    for (i=0;i<NCHANNELS;++i){
      corr[i]=(double*) malloc(sizeof(double)*size);
    }
    size_old=size;
  }
 
  if (!init){
    for(k=0; k<NCHANNELS; k++){
      for(i=0; i<nm*GLB_T*n_mom_tot; i++){
	corr[k][i] = 0.;
      }
    }
    init = 1;
  }
}


void measure_mesons(spinor_field* psi0, spinor_field* eta, int nm, int tau){
  init_corrs(nm,1);
  lprintf("measure_mesons",50,"measure default mesons");
  measure_mesons_core(psi0, psi0, eta, nm, tau, 1, 0,GLB_T);
}

void measure_mesons_ext(spinor_field* psi0, spinor_field* eta, int nm, int tau,int begin){
  init_corrs(nm*3,1);
  lprintf("measure_mesons",50,"measure extended mesons segment %d",begin);
  measure_mesons_core(psi0, psi0,eta, nm, tau, 1, GLB_T*begin, 3*GLB_T);
}


void measure_point_mesons_momenta(spinor_field* psi0, spinor_field* eta, int nm, int tau, int n_mom){
  init_corrs(nm,n_mom);
  lprintf("measure_mesons",50,"measure point mesons with momenta");
   measure_mesons_core(psi0, psi0, eta, nm, tau, n_mom, 0, GLB_T);
}

void measure_point_mesons_momenta_ext(spinor_field* psi0, spinor_field* eta, int nm, int tau, int n_mom, int begin){
  init_corrs(3*nm,n_mom);
  lprintf("measure_mesons",50,"measure extended point mesons with moment");
  measure_mesons_core(psi0, psi0, eta, nm, tau, n_mom, GLB_T*begin, 3*GLB_T);
}

/* psi0 sequential propagator from tf, psi point to all propagator from ti */
void measure_formfactors(spinor_field* psi0, spinor_field* psi1, spinor_field* eta, int nm, int ti, int tf, int n_mom){
  int i;
  for (i=0;i<NCHANNELS-1;++i){
    measure_channels[i]=0;
  }
  init_corrs(nm,n_mom);

  lprintf("MEASURE_FORMFACTORS",50,"psi0 = sequential prop, psi1 = point prop");

  measure_channels[_g5]=1;
  measure_mesons_core(psi1, psi1, eta, nm, ti, n_mom, 0, GLB_T);
  measure_channels[_g5]=0;
  measure_channels[_pipig_ff_re]=1;
  measure_channels[_pipig_ff_im]=1;
  measure_channels[_pipig_conserved_ff]=1;
  measure_mesons_core(psi0, psi1, eta, nm, ti, n_mom, 0, GLB_T);
  measure_channels[_g5]=1;
}

void measure_formfactors_ext(spinor_field* psi0, spinor_field* psi1, spinor_field* eta, int nm, int ti, int tf, int n_mom, int begin){
  int i;
  for (i=0;i<NCHANNELS-1;++i){
    measure_channels[i]=0;
  }
  init_corrs(2*nm,n_mom);
  if (begin){
   measure_channels[_g5]=1;
   measure_mesons_core(psi1, psi1, eta, nm, ti, n_mom, 0, 2*GLB_T);
   measure_channels[_g5]=0;
   measure_channels[_pipig_ff_re]=1;
   measure_channels[_pipig_ff_im]=1;
   measure_mesons_core(psi0, psi1, eta, nm, ti, n_mom, 0, 2*GLB_T);
   measure_channels[_g5]=1;
  }
  else{
   measure_channels[_g5]=1;
   measure_mesons_core(psi1, psi1, eta, nm, ti, n_mom, GLB_T, 2*GLB_T);
   measure_channels[_g5]=0;
   measure_channels[_pipig_ff_re]=1;
   measure_channels[_pipig_ff_im]=1;
   measure_mesons_core(psi0, psi1, eta, nm, ti, n_mom, GLB_T, 2*GLB_T);
   measure_channels[_g5]=1;
  }
}


void measure_formfactors_conserved(spinor_field* psi0, spinor_field* psi1, spinor_field* eta, int nm, int ti, int tf, int n_mom, int begin){
  init_corrs(nm,n_mom);
  measure_channels[_pipig_conserved_ff]=1;
  measure_nonlocal_core(psi1, psi1, eta, nm, ti, n_mom, 0, GLB_T);    
}


static void print_corr_core(int channel,int lt,int conf, int nm, double* mass, char* label, int n_mom){
  int i,t,px,py,pz;
  for (px=0;px<n_mom;++px) for (py=0;py<n_mom;++py) for (pz=0;pz<n_mom;++pz){ 
	/*	if (is_complex){
	  for(i=0; i<nm; i++) { 
	    if (n_mom>1){
	      lprintf("MAIN",0,"conf #%d mass=%2.6f %s %s %s momentum(%d,%d,%d)= ",conf,mass[i],label,channel_types[channel],channel_names[channel],px,py,pz);
	    }else{
	      lprintf("MAIN",0,"conf #%d mass=%2.6f %s %s %s= ",conf,mass[i],label,channel_types[channel],channel_names[channel]);
	    }
	    for(t=0;t<lt;++t) { 
	      lprintf("MAIN",0,"( %e %e ) ",corr[ channel ][corr_ind(px,py,pz,n_mom,t+i*lt,nm)],corr_im[channel][corr_ind(px,py,pz,n_mom,t+i*lt,nm)]); 
	    } 
	    lprintf("MAIN",0,"\n"); 
	    fflush(stdout); 
	  } 
	}
	else{*/
	  for(i=0; i<nm; i++) {
	    if (n_mom>1){
	      lprintf("MAIN",0,"conf #%d mass=%2.6f %s %s %s momentum(%d,%d,%d)= ",conf,mass[i],label,channel_types[channel],channel_names[channel],px,py,pz); 
	    }
	    else{
	      lprintf("MAIN",0,"conf #%d mass=%2.6f %s %s %s= ",conf,mass[i],label,channel_types[channel],channel_names[channel]); 
            }
	    for(t=0;t<lt;++t) lprintf("MAIN",0,"%e ",corr[channel][corr_ind(px,py,pz,n_mom,t+i*lt,nm)]);
	    lprintf("MAIN",0,"\n"); 
	    fflush(stdout); 
	    //	  } 
	}
      }
}

static void print_corr(int lt,int conf, int nm, double* mass, char* label, int n_mom){
  int i;
  for(i = 0; i<NCHANNELS; ++i){
    if(measure_channels[i]){
      print_corr_core(i,lt,conf,nm,mass,label,n_mom);
    }
  }
}

void print_mesons(double norm, int conf, int nm, double* mass, char* label){
  int i,k;
  if (init !=0 ){
    for(k=0; k<NCHANNELS; k++) {
      global_sum(corr[k],GLB_T*nm);
      for(i=0; i<nm*GLB_T; i++)
	corr[k][i] *= -((1./norm)/GLB_VOL3);
    }
  }
  init = 0;
  print_corr(GLB_T,conf,nm,mass,label,1);
}

void print_mesons_ext(double norm, int conf, int nm, double* mass, char* label){
  int i,k;
  if (init !=0 ){
    for(k=0; k<NCHANNELS; k++) {
      global_sum(corr[k],3*GLB_T*nm);
      for(i=0; i<3*nm*GLB_T; i++)
	corr[k][i] *= -((1./norm)/GLB_VOL3);
    }
  }
  init = 0;
  print_corr(3*GLB_T,conf,nm,mass,label,1);
}


void print_mesons_momenta(int conf, int nm, double* mass, int n_mom, char* label){
  int i,k;
  if (init){
    for(k=0; k<NCHANNELS; k++) {
      global_sum(corr[k],GLB_T*nm*n_mom*n_mom*n_mom);
      for(i=0; i<nm*GLB_T*n_mom*n_mom*n_mom; i++){
	corr[k][i] *= -((1.)/GLB_VOL3);
      }
    }
  }
  init = 0;
  print_corr(GLB_T,conf,nm,mass,label,n_mom);
}


void print_mesons_momenta_ext(int conf, int nm, double* mass, int n_mom, char* label){
  int i,k;
  if (init){
    for(k=0; k<NCHANNELS; k++) {
      global_sum(corr[k],3*GLB_T*nm*n_mom*n_mom*n_mom);
      for(i=0; i<3*nm*GLB_T*n_mom*n_mom*n_mom; i++){
	corr[k][i] *= -((1.)/GLB_VOL3);
      }
    }
  }
  init = 0;
  print_corr(3*GLB_T,conf,nm,mass,label,n_mom);
}

void print_formfactor(int conf, int nm, double* mass, int n_mom, char* label, int tf){
  int i,k;
  if (init){
    for(k=0; k<NCHANNELS; k++) {
      global_sum(corr[k],GLB_T*nm*n_mom*n_mom*n_mom);
      //      global_sum(corr_im[k],GLB_T*nm*n_mom*n_mom*n_mom);
      for(i=0; i<nm*GLB_T*n_mom*n_mom*n_mom; i++){
	corr[k][i] *= -((1.)/GLB_VOL3);
	//	corr_im[k][i] *= -((1.)/GLB_VOL3);
      }
    }
  }
  init = 0;
  char str[256]; sprintf(str, "%s tf %d", label, tf); //Modify label
  print_corr(GLB_T,conf,nm,mass,str,n_mom);
}

void print_formfactor_ext(int conf, int nm, double* mass, int n_mom, char* label, int tf){
  int i,k;
  if (init){
    for(k=0; k<NCHANNELS; k++) {
      global_sum(corr[k],2*GLB_T*nm*n_mom*n_mom*n_mom);
      for(i=0; i<2*nm*GLB_T*n_mom*n_mom*n_mom; i++){
        corr[k][i] *= -((1.)/GLB_VOL3);
      }
    }
  }
  init = 0;
  char str[256]; sprintf(str, "%s tf %d", label, tf); //Modify label
  print_corr(2*GLB_T,conf,nm,mass,str,n_mom);
}


#undef corr_ind
