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

enum { _g5_ff= 0,_pipig_ff_re,_pipig_ff_im,_pipig_conserved_ff,NCHANNELS_FF };

static char* ff_channel_names[NCHANNELS_FF]={"g5_ff","pipig_ff_re","pipig_ff_im","pipig_conserved_ff"};

static double* corr[NCHANNELS_FF];

static int init = 0;

#define corr_ind(px,py,pz,n_mom,tc,nm) ((px)*(n_mom)*(n_mom)*(lt)*(nm)+(py)*(n_mom)*(lt)*(nm)+(pz)*(lt)*(nm) +(tc))

/* spinor_fields* are nmx4xNF arrays of spinor_field ordered([nm][color][spinor])*/
static void measure_formfactor_core(spinor_field* psi0, spinor_field* psi1, spinor_field* eta, int nm, int tau, int tf, int n_mom, int offset,int lt, int* pt){
  int i,ix,t,x,y,z,beta,px,py,pz,tc,a,ixmu;
  double pdotx,cpdotx,spdotx;
  complex tr;
  suNf_propagator sp0,sp1,Usp,spf,sptmp1,sptmp2,sptmp3,spdag;
  suNf *u1;




  lprintf("measure_mesons_core",50,"Measuring channels: ");
  for (i=NCHANNELS_FF-1;i<NCHANNELS_FF;++i){
    lprintf("",50," %s",ff_channel_names[i]);
  }
  lprintf("",50,"\n");

  for (px=0;px<n_mom;++px) for (py=0;py<n_mom;++py) for (pz=0;pz<n_mom;++pz){
	for(i=0; i<nm; i++) {
	  for (t=0; t<T; t++) {	 
	    tc = (zerocoord[0]+t+GLB_T-tau)%GLB_T+i*GLB_T+offset;
	    for (x=0; x<X; x++) for (y=0; y<Y; y++) for (z=0; z<Z; z++) {
		  ix=ipt(t,x,y,z);					
		  pdotx = 2.0*PI*(((double) px)*(x+zerocoord[1]-pt[1])/GLB_X + ((double) py)*(y+zerocoord[2]-pt[2])/GLB_Y + ((double) pz)*(z+zerocoord[3]-pt[3])/GLB_Z);
		  cpdotx=cos(pdotx);
		  spdotx=sin(pdotx);
		  for (a=0;a<NF;++a){
		    for (beta=0;beta<4;beta++){ 
		      _propagator_assign(sp0, *_FIELD_AT(&psi0[a*4*nm+beta*nm+i],ix),a,beta);
		      _propagator_assign(sp1, *_FIELD_AT(&psi1[a*4*nm+beta*nm+i],ix),a,beta);
		    }
		  }
                  _propagator_dagger(spdag,sp1);

		  //Pion
		  _propagator_muldag_trace(tr,sp1,sp1);
		  corr[_g5_ff][corr_ind(px,py,pz,n_mom,tc,nm)] += (cpdotx*tr.re + spdotx*tr.im);		  

		  // Local vector current (pipig)
		  //Sequential source  
		  _g5g0_propagator(sptmp2, sp0);
		  _propagator_mul(sptmp1, sptmp2, spdag); // g5 g0 G(x,0) S(y,0)^dagger
		  _propagator_trace(tr,sptmp1);
		  corr[_pipig_ff_re][corr_ind(px,py,pz,n_mom,tc,nm)] += +0.5*(cpdotx*tr.re + spdotx*tr.im);
		  corr[_pipig_ff_im][corr_ind(px,py,pz,n_mom,tc,nm)] += +0.5*(spdotx*tr.re - cpdotx*tr.im);

		  //2nd contraction
		  _g5g0_propagator(sptmp1, sp1);
		  _propagator_dagger(sptmp2,sp0); //G(y,0)^dagger
		  _propagator_mul(spf, sptmp2, sptmp1); // G(y,0)^dagger g5 g0 S(y,0)
		  _propagator_trace(tr,spf);
		  corr[_pipig_ff_re][corr_ind(px,py,pz,n_mom,tc,nm)] += -0.5*(cpdotx*tr.re + spdotx*tr.im);
		  corr[_pipig_ff_im][corr_ind(px,py,pz,n_mom,tc,nm)] += -0.5*(spdotx*tr.re - cpdotx*tr.im);
		  
		  //Conserved vector current

		  u1 = _4FIELD_AT(u_gauge_f,ix,0);
		  ixmu = iup(ix,0);
		  
		  // Tr [ 1/2 S^(0,y+\mu) g_5(1+g_0) U^(y) S(y,x) 
		  _suNf_inverse_prop_multiply(Usp,*u1,sp0);
		  sptmp1=Usp;
		  _g0_propagator(sptmp2,sptmp1);
		  _propagator_add(sptmp1,sptmp1,sptmp2);
		  _g5_propagator(sptmp2,sptmp1);
		  for (a=0;a<NF;a++) for (beta=0;beta<4;beta++){ 
		      _propagator_assign(spf, *_FIELD_AT(&psi1[4*a*nm+beta*nm+i],ixmu), a,beta); 
		    }
		  _propagator_dagger(sptmp1,spf);
		  _propagator_mul(sptmp3,sptmp1,sptmp2);	
		    

		  // Tr [ -1/2 S^(0,y) g_5(1-g_0) U(y) S(y+mu,x) 
		  for (a=0;a<NF;a++) for (beta=0;beta<4;beta++){ 
		      _propagator_assign(spf, *_FIELD_AT(&psi0[4*a*nm+beta*nm+i],ixmu), a,beta); 
		    }
		  _suNf_prop_multiply(Usp,*u1,spf);
		  sptmp1=Usp;
		  _g0_propagator(sptmp2,sptmp1);
		  _propagator_sub(sptmp1,sptmp1,sptmp2);
		  _g5_propagator(sptmp2,sptmp1);
		  _propagator_mul(sptmp1,spdag,sptmp2);	
		  _propagator_sub(sptmp2,sptmp3,sptmp1);
		  _propagator_trace(tr,sptmp2);    
		  corr[_pipig_conserved_ff][corr_ind(px,py,pz,n_mom,tc,nm)] += 0.5*(cpdotx*tr.re + spdotx*tr.im);
		} //END SPATIAL LOOP
	  } //END T LOOP
	} //END MASS LOOP
      } //END MOMENTUM LOOP
  lprintf("measure_formfactor_core",50,"Measuring DONE! ");
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
      for (i=0;i<NCHANNELS_FF;++i){
	free(corr[i]);
      }
    }
    for (i=0;i<NCHANNELS_FF;++i){
      corr[i]=(double*) malloc(sizeof(double)*size);
    }
    size_old=size;
  }
  if (!init){
    for(k=0; k<NCHANNELS_FF; k++){
      for(i=0; i<nm*GLB_T*n_mom_tot; i++){
	corr[k][i] = 0.;
      }
    }
    init = 1;
  }
}


/* psi0 sequential propagator from tf, psi point to all propagator from ti */
void measure_formfactors(spinor_field* psi0, spinor_field* psi1, spinor_field* eta, int nm, int ti, int tf, int n_mom, int *pt){
  init_corrs(nm,n_mom);
  lprintf("MEASURE_FORMFACTORS",50,"psi0 = sequential prop, psi1 = point prop");
  measure_formfactor_core(psi0, psi1, eta, nm, ti, tf, n_mom, 0, GLB_T,pt);  
}

/*void measure_formfactors_ext(spinor_field* psi0, spinor_field* psi1, spinor_field* eta, int nm, int ti, int tf, int n_mom, int begin){
  init_corrs(3*nm,n_mom);
  lprintf("MEASURE_FORMFACTORS",50,"psi0 = sequential prop, psi1 = point prop");
  measure_formfactor_core(psi0, psi1, eta, nm, ti, tf, n_mom, GLB_T*begin, 3*GLB_T);  
  }*/

static void print_corr_core(int channel,int lt,int conf, int nm, double* mass, char* label, int n_mom){
  int i,t,px,py,pz;
  for (px=0;px<n_mom;++px) for (py=0;py<n_mom;++py) for (pz=0;pz<n_mom;++pz){ 
	  for(i=0; i<nm; i++) {
	    if (n_mom>1){
	      lprintf("MAIN",0,"conf #%d mass=%2.6f %s FORM_FACTOR %s momentum(%d,%d,%d)= ",conf,mass[i],label,ff_channel_names[channel],px,py,pz); 
	    }
	    else{
	      lprintf("MAIN",0,"conf #%d mass=%2.6f %s FORM_FACTOR %s= ",conf,mass[i],label,ff_channel_names[channel]); 
            }
	    for(t=0;t<lt;++t) lprintf("MAIN",0,"%e ",corr[channel][corr_ind(px,py,pz,n_mom,t+i*lt,nm)]);
	    lprintf("MAIN",0,"\n"); 
	    fflush(stdout); 
	  }
      }
}

static void print_corr(int lt,int conf, int nm, double* mass, char* label, int n_mom){
  int i;
  for(i = 0; i<NCHANNELS_FF; ++i){
    print_corr_core(i,lt,conf,nm,mass,label,n_mom); 
  }
}


void print_formfactor(int conf, int nm, double* mass, int n_mom, char* label, int tf){
  int i,k;
  if (init){
    for(k=0; k<NCHANNELS_FF; k++) {
      global_sum(corr[k],GLB_T*nm*n_mom*n_mom*n_mom);
      for(i=0; i<nm*GLB_T*n_mom*n_mom*n_mom; i++){
	corr[k][i] *= -((1.)/GLB_VOL3);
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
    for(k=0; k<NCHANNELS_FF; k++) {
      global_sum(corr[k],3*GLB_T*nm*n_mom*n_mom*n_mom);
      for(i=0; i<3*nm*GLB_T*n_mom*n_mom*n_mom; i++){
        corr[k][i] *= -((1.)/GLB_VOL3);
      }
    }
  }
  init = 0;
  char str[256]; sprintf(str, "%s tf %d", label, tf); //Modify label
  print_corr(3*GLB_T,conf,nm,mass,str,n_mom);
}


