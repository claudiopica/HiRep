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
  for (i=0;i<NCHANNELS;++i){
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
		  }
		  if(measure_channels[_id]){
		    _g5_spinmatrix(sm1, sma);
		    _spinmatrix_g5(sm2, smb);
		    _spinmatrix_mul_trace(tr, sm1, sm2);
		    corr[_id][corr_ind(px,py,pz,n_mom,tc,nm)] += cpdotx*tr.re;
		  }
		  if(measure_channels[_g0]){
		    _g5g0_spinmatrix(sm1, sma);
		    _spinmatrix_g5g0(sm2, smb);
		    _spinmatrix_mul_trace(tr, sm1, sm2);
		    corr[_g0][corr_ind(px,py,pz,n_mom,tc,nm)] += -cpdotx*tr.re;	
		  }
		  if(measure_channels[_g1]){
		    _g5g1_spinmatrix(sm1, sma);
		    _spinmatrix_g5g1(sm2, smb);
		    _spinmatrix_mul_trace(tr, sm1, sm2);
		    corr[_g1][corr_ind(px,py,pz,n_mom,tc,nm)] += -cpdotx*tr.re;
		  }      
		  if(measure_channels[_g2]){
		    _g5g2_spinmatrix(sm1, sma);
		    _spinmatrix_g5g2(sm2, smb);
		    _spinmatrix_mul_trace(tr, sm2, sm1);
		    corr[_g2][corr_ind(px,py,pz,n_mom,tc,nm)] += cpdotx*tr.re;
		  }
		  if(measure_channels[_g3]){
		    _g5g3_spinmatrix(sm1, sma);
		    _spinmatrix_g5g3(sm2, smb);
		    _spinmatrix_mul_trace(tr, sm1, sm2);
		    corr[_g3][corr_ind(px,py,pz,n_mom,tc,nm)] += -cpdotx*tr.re;
		  }
		  if(measure_channels[_g0g5]){
		    _g0_spinmatrix(sm1, sma);
		    _spinmatrix_g0(sm2, smb);
		    _spinmatrix_mul_trace(tr, sm1, sm2);
		    corr[_g0g5][corr_ind(px,py,pz,n_mom,tc,nm)] += -cpdotx*tr.re;
		  }
		  if(measure_channels[_g5g1]){
		    _g1_spinmatrix(sm1, sma);
		    _spinmatrix_g1(sm2, smb);
		    _spinmatrix_mul_trace(tr, sm1, sm2);
		    corr[_g5g1][corr_ind(px,py,pz,n_mom,tc,nm)] += -cpdotx*tr.re;
		  }	      
		  if(measure_channels[_g5g2]){
		    _g2_spinmatrix(sm1, sma);
		    _spinmatrix_g2(sm2, smb);
		    _spinmatrix_mul_trace(tr, sm1, sm2);
		    corr[_g5g2][corr_ind(px,py,pz,n_mom,tc,nm)] += cpdotx*tr.re;
		  }
		  if(measure_channels[_g5g3]){
		    _g3_spinmatrix(sm1, sma);
		    _spinmatrix_g3(sm2, smb);
		    _spinmatrix_mul_trace(tr, sm1, sm2);
		    corr[_g5g3][corr_ind(px,py,pz,n_mom,tc,nm)] += -cpdotx*tr.re;
		  }
		  if(measure_channels[_g0g1]){
		    _g5g0g1_spinmatrix(sm1, sma);
		    _spinmatrix_g5g0g1(sm2, smb);
		    _spinmatrix_mul_trace(tr, sm1, sm2);
		    corr[_g0g1][corr_ind(px,py,pz,n_mom,tc,nm)] += -cpdotx*tr.re;
		  }
		  if(measure_channels[_g0g2]){
		    _g5g0g2_spinmatrix(sm1, sma);
		    _spinmatrix_g5g0g2(sm2, smb);
		    _spinmatrix_mul_trace(tr, sm1, sm2);
		    corr[_g0g2][corr_ind(px,py,pz,n_mom,tc,nm)] += cpdotx*tr.re;
		  }
		  if(measure_channels[_g0g3]){
		    _g5g0g3_spinmatrix(sm1, sma);
		    _spinmatrix_g5g0g3(sm2, smb);
		    _spinmatrix_mul_trace(tr, sm1, sm2);
		    corr[_g0g3][corr_ind(px,py,pz,n_mom,tc,nm)] += -cpdotx*tr.re;
		  }
		  if(measure_channels[_g0g5g1]){
		    _g0g1_spinmatrix(sm1, sma);
		    _spinmatrix_g0g1(sm2, smb);
		    _spinmatrix_mul_trace(tr, sm1, sm2);
		    corr[_g0g5g1][corr_ind(px,py,pz,n_mom,tc,nm)] += +cpdotx*tr.re;
		  }
		  if(measure_channels[_g0g5g2]){
		    _g0g2_spinmatrix(sm1, sma);
		    _spinmatrix_g0g2(sm2, smb);
		    _spinmatrix_mul_trace(tr, sm1, sm2);
		    corr[_g0g5g2][corr_ind(px,py,pz,n_mom,tc,nm)] += -cpdotx*tr.re;
		  }
		  if(measure_channels[_g0g5g3]){
		    _g0g3_spinmatrix(sm1, sma);
		    _spinmatrix_g0g3(sm2, smb);
		    _spinmatrix_mul_trace(tr, sm1, sm2);
		    corr[_g0g5g3][corr_ind(px,py,pz,n_mom,tc,nm)] += +cpdotx*tr.re;
		  }
		  if(measure_channels[_g5_g0g5_re]){
		    _spinmatrix_g0(sm2, smb);
		    _spinmatrix_mul_trace(tr, sma, sm2);
		    corr[_g5_g0g5_re][corr_ind(px,py,pz,n_mom,tc,nm)] += cpdotx*tr.re+spdotx*tr.im;
		  }
		  if(measure_channels[_id_disc]){
		    for (beta=0;beta<4;beta++){
		      _spinor_prod_re_f(tr.re,*_FIELD_AT(&psi0[beta*nm+i],ix),*_FIELD_AT(&eta[beta],ix));
		      corr[_id_disc][zerocoord[0]+t+i*lt+offset] += -tr.re;
		    }
		  }
		  if(measure_channels[_g5_disc]){
		    for (beta=0;beta<4;beta++){
		      _spinor_g5_prod_re_f(tr.re,*_FIELD_AT(&psi0[beta*nm+i],ix),*_FIELD_AT(&eta[beta],ix));
		      corr[_g5_disc][zerocoord[0]+t+i*lt+offset] += tr.re;
		    }
		  }
		} //END SPATIAL LOOP
	  } //END TIME LOOP
	} //END MASS LOOP
      } //END MOMENTUM LOOP
  lprintf("measure_mesons_core",50,"Measuring DONE! ");
}

/* spinor_fields* are 4xnm arrays of spinor_field ordered([color][spinor])*/
static void measure_discon_noise_core(spinor_field* psi0, spinor_field* psi1, spinor_field* eta, int nm, int n_mom, int offset,int lt){
  int i,ix,t,x,y,z,beta,px,py,pz,tc,tau,tm;
  double pdotx,cpdotx,spdotx;
  complex tr;
  suNf_spin_matrix sma,smb,smc, sm1,sm2;
  double **bracket=(double**) malloc(sizeof(double*)*NCHANNELS); //b(t) = sum_vec{x} <eta(x)| Gamma | psi(x) >
  for(i=0;i<NCHANNELS;i++){ 
	bracket[i]=(double*) malloc(sizeof(double)*lt); 
	for(t=0;t<lt;t++){ bracket[i][t] = 0; }
  }

  lprintf("measure_discon_noise_core",50,"Measuring channels: ");
  for (i=0;i<NCHANNELS;++i){
    if (measure_channels[i]) lprintf("",50," %s",channel_names[i]);
  }
  lprintf("",50,"\n");
  for (px=0;px<n_mom;++px) for (py=0;py<n_mom;++py) for (pz=0;pz<n_mom;++pz){
	for(i=0; i<nm; i++) {
	  for (t=0; t<T; t++) {	 
	    tc = (zerocoord[0]+t)+i*GLB_T+offset;
	    for (x=0; x<X; x++) for (y=0; y<Y; y++) for (z=0; z<Z; z++) { 
		  ix=ipt(t,x,y,z);					
		  pdotx = 2.0*PI*(((double) px)*(x+zerocoord[1])/GLB_X + ((double) py)*(y+zerocoord[2])/GLB_Y + ((double) pz)*(z+zerocoord[3])/GLB_Z);
		  cpdotx=cos(pdotx);
		  spdotx=sin(pdotx);
  
		  for (beta=0;beta<4;beta++){ 
		    _spinmatrix_assign_row(sma, *_FIELD_AT(&psi0[beta*nm+i],ix), beta);
		    _spinmatrix_assign_row(smb, *_FIELD_AT(&psi1[beta*nm+i],ix), beta); 
		    _spinmatrix_assign_row(smc, *_FIELD_AT(&eta[beta],ix), beta); 
		  }
		  if (measure_channels[_g5]){
		    //_spinmatrix_mul_trace(tr, sma, smb);
		    //corr[_g5][corr_ind(px,py,pz,n_mom,tc,nm)]+= -cpdotx*tr.re;
		    _g5_spinmatrix(sm1, sma);
		    _spinmatrix_mul_trace(tr, smc, sm1);
		    bracket[_g5][(zerocoord[0]+t)] += -tr.re;
		  }
		  /*if(measure_channels[_id]){
		    _g5_spinmatrix(sm1, sma);
		    _spinmatrix_g5(sm2, smb);
		    _spinmatrix_mul_trace(tr, sm1, sm2);
		    corr[_id][corr_ind(px,py,pz,n_mom,tc,nm)] += cpdotx*tr.re;
		  }
		  if(measure_channels[_g0]){
		    _g5g0_spinmatrix(sm1, sma);
		    _spinmatrix_g5g0(sm2, smb);
		    _spinmatrix_mul_trace(tr, sm1, sm2);
		    corr[_g0][corr_ind(px,py,pz,n_mom,tc,nm)] += -cpdotx*tr.re;	
		  }
		  if(measure_channels[_g1]){
		    _g5g1_spinmatrix(sm1, sma);
		    _spinmatrix_g5g1(sm2, smb);
		    _spinmatrix_mul_trace(tr, sm1, sm2);
		    corr[_g1][corr_ind(px,py,pz,n_mom,tc,nm)] += -cpdotx*tr.re;
		  }      
		  if(measure_channels[_g2]){
		    _g5g2_spinmatrix(sm1, sma);
		    _spinmatrix_g5g2(sm2, smb);
		    _spinmatrix_mul_trace(tr, sm2, sm1);
		    corr[_g2][corr_ind(px,py,pz,n_mom,tc,nm)] += cpdotx*tr.re;
		  }
		  if(measure_channels[_g3]){
		    _g5g3_spinmatrix(sm1, sma);
		    _spinmatrix_g5g3(sm2, smb);
		    _spinmatrix_mul_trace(tr, sm1, sm2);
		    corr[_g3][corr_ind(px,py,pz,n_mom,tc,nm)] += -cpdotx*tr.re;
		  }
		  if(measure_channels[_g0g5]){
		    _g0_spinmatrix(sm1, sma);
		    _spinmatrix_g0(sm2, smb);
		    _spinmatrix_mul_trace(tr, sm1, sm2);
		    corr[_g0g5][corr_ind(px,py,pz,n_mom,tc,nm)] += -cpdotx*tr.re;
		  }
		  if(measure_channels[_g5g1]){
		    _g1_spinmatrix(sm1, sma);
		    _spinmatrix_g1(sm2, smb);
		    _spinmatrix_mul_trace(tr, sm1, sm2);
		    corr[_g5g1][corr_ind(px,py,pz,n_mom,tc,nm)] += -cpdotx*tr.re;
		  }	      
		  if(measure_channels[_g5g2]){
		    _g2_spinmatrix(sm1, sma);
		    _spinmatrix_g2(sm2, smb);
		    _spinmatrix_mul_trace(tr, sm1, sm2);
		    corr[_g5g2][corr_ind(px,py,pz,n_mom,tc,nm)] += cpdotx*tr.re;
		  }
		  if(measure_channels[_g5g3]){
		    _g3_spinmatrix(sm1, sma);
		    _spinmatrix_g3(sm2, smb);
		    _spinmatrix_mul_trace(tr, sm1, sm2);
		    corr[_g5g3][corr_ind(px,py,pz,n_mom,tc,nm)] += -cpdotx*tr.re;
		  }
		  if(measure_channels[_g0g1]){
		    _g5g0g1_spinmatrix(sm1, sma);
		    _spinmatrix_g5g0g1(sm2, smb);
		    _spinmatrix_mul_trace(tr, sm1, sm2);
		    corr[_g0g1][corr_ind(px,py,pz,n_mom,tc,nm)] += -cpdotx*tr.re;
		  }
		  if(measure_channels[_g0g2]){
		    _g5g0g2_spinmatrix(sm1, sma);
		    _spinmatrix_g5g0g2(sm2, smb);
		    _spinmatrix_mul_trace(tr, sm1, sm2);
		    corr[_g0g2][corr_ind(px,py,pz,n_mom,tc,nm)] += cpdotx*tr.re;
		  }
		  if(measure_channels[_g0g3]){
		    _g5g0g3_spinmatrix(sm1, sma);
		    _spinmatrix_g5g0g3(sm2, smb);
		    _spinmatrix_mul_trace(tr, sm1, sm2);
		    corr[_g0g3][corr_ind(px,py,pz,n_mom,tc,nm)] += -cpdotx*tr.re;
		  }
		  if(measure_channels[_g0g5g1]){
		    _g0g1_spinmatrix(sm1, sma);
		    _spinmatrix_g0g1(sm2, smb);
		    _spinmatrix_mul_trace(tr, sm1, sm2);
		    corr[_g0g5g1][corr_ind(px,py,pz,n_mom,tc,nm)] += +cpdotx*tr.re;
		  }
		  if(measure_channels[_g0g5g2]){
		    _g0g2_spinmatrix(sm1, sma);
		    _spinmatrix_g0g2(sm2, smb);
		    _spinmatrix_mul_trace(tr, sm1, sm2);
		    corr[_g0g5g2][corr_ind(px,py,pz,n_mom,tc,nm)] += -cpdotx*tr.re;
		  }
		  if(measure_channels[_g0g5g3]){
		    _g0g3_spinmatrix(sm1, sma);
		    _spinmatrix_g0g3(sm2, smb);
		    _spinmatrix_mul_trace(tr, sm1, sm2);
		    corr[_g0g5g3][corr_ind(px,py,pz,n_mom,tc,nm)] += +cpdotx*tr.re;
		  }
		  if(measure_channels[_g5_g0g5_re]){
		    _spinmatrix_g0(sm2, smb);
		    _spinmatrix_mul_trace(tr, sma, sm2);
		    corr[_g5_g0g5_re][corr_ind(px,py,pz,n_mom,tc,nm)] += cpdotx*tr.re+spdotx*tr.im;
		  }
		  if(measure_channels[_id_disc]){
		    for (beta=0;beta<4;beta++){
		      _spinor_prod_re_f(tr.re,*_FIELD_AT(&psi0[beta*nm+i],ix),*_FIELD_AT(&eta[beta],ix));
		      corr[_id_disc][zerocoord[0]+t+i*lt+offset] += -tr.re;
		    }
		  }
		  if(measure_channels[_g5_disc]){
		    for (beta=0;beta<4;beta++){
		      _spinor_g5_prod_re_f(tr.re,*_FIELD_AT(&psi0[beta*nm+i],ix),*_FIELD_AT(&eta[beta],ix));
		      corr[_g5_disc][zerocoord[0]+t+i*lt+offset] += tr.re;
		    }
		  }*/
		} //END SPATIAL LOOP
	  } //END TIME LOOP
	    double *bucket=(double*) malloc(sizeof(double)*lt); 
	    #ifdef WITH_MPI
	     MPI_Gather(bracket[ _g5 ] + zerocoord[0], T, MPI_DOUBLE, bucket, T, MPI_DOUBLE, 0, GLB_COMM);
	    #else
	     for(t=0;t<T;++t){ bucket[t] = bracket[ _g5 ][ t ]; }
	    #endif

            for(t=0;t<lt;++t){ 
		for(tau=0;tau<lt;++tau){ 
       		 tm = (t+tau)%lt;
		 corr[_g5][t+i*lt+offset] += bucket[tau] * bucket[tm];
	        }
            }

	} //END MASS LOOP
      } //END MOMENTUM LOOP
  lprintf("measure_discon_noise_core",50,"Measuring DONE! ");
}


static double hmass_pre;
static void D_pre(spinor_field *out, spinor_field *in){  Dphi_eopre(hmass_pre,out,in); }
static void Ddag_pre(spinor_field *out, spinor_field *in, spinor_field *ttmp){ 
      spinor_field_copy_f( ttmp, in ); 
      spinor_field_g5_assign_f( ttmp ); 
      H_pre(out, ttmp);
}
/* spinor_fields* are 4xnm arrays of spinor_field ordered([color][spinor])*/
static void measure_discon_lma_core(int Nev, int nm, int n_mom, int offset,int lt, double *m){
  int i,ix,t,x,y,z,beta,px,py,pz,tc,tau,n;
  double pdotx,cpdotx,spdotx;
  complex tr;

  spinor_field* vi; vi = alloc_spinor_field_f(1,&glat_even);
  spinor_field* Dvi; Dvi = alloc_spinor_field_f(1,&glat_even);
  double li;

  lprintf("measure_discon_lma_core",10,"Measuring channels: ");
  for (i=0;i<NCHANNELS;++i){
    if (measure_channels[i]) lprintf("",10," %s",channel_names[i]);
  }
  lprintf("",10,"\n");

  hmass_pre = m[0];

for(n=0;n<Nev;n++){

  copy_evec(n, vi, &li);
  D_pre(Dvi, vi);
  //spinor_field_sub_assign_f(tmp1,Dvi);

  //lprintf("DISCON",0,"%g %g %g\n",hmass_pre, spinor_field_sqnorm_f(vi), spinor_field_sqnorm_f(tmp1) );
  for (px=0;px<n_mom;++px) for (py=0;py<n_mom;++py) for (pz=0;pz<n_mom;++pz){
	for(i=0; i<nm; i++) {
	  for (t=0; t<T; t++) {	 
	    tc = (zerocoord[0]+t)+i*GLB_T+offset;
	    for (x=0; x<X; x++) for (y=0; y<Y; y++) for (z=0; z<Z; z++) { 
		  ix=ipt(t,x,y,z);					
		  pdotx = 2.0*PI*(((double) px)*(x+zerocoord[1])/GLB_X + ((double) py)*(y+zerocoord[2])/GLB_Y + ((double) pz)*(z+zerocoord[3])/GLB_Z);
		  cpdotx=cos(pdotx);
		  spdotx=sin(pdotx);
  

		  
		  if(measure_channels[_id_disc]){
		      _spinor_prod_re_f(tr.re,*_FIELD_AT(Dvi,ix),*_FIELD_AT(vi,ix));
		      corr[_id_disc][zerocoord[0]+t+i*lt+offset] +=  -(4. + hmass_pre) * tr.re/li;
		  }
		  if(measure_channels[_g5_disc]){
		      _spinor_g5_prod_re_f(tr.re,*_FIELD_AT(Dvi,ix),*_FIELD_AT(vi,ix));
		      corr[_g5_disc][zerocoord[0]+t+i*lt+offset] +=  (4. + hmass_pre) *tr.re/li;
		  }

		
		

		} //END SPATIAL LOOP
	  } //END TIME LOOP
	} //END MASS LOOP
      } //END MOMENTUM LOOP
    }// END NEV LOOP
  free_spinor_field_f(vi);
  free_spinor_field_f(Dvi);
  lprintf("measure_discon_noise_core",50,"Measuring DONE! ");
}

/* spinor_fields* are nmx4xNF arrays of spinor_field ordered([nm][color][spinor])*/
static void measure_formfactor_core(spinor_field* psi0, spinor_field* psi1, spinor_field* eta, int nm, int tau, int tf, int n_mom, int offset,int lt){
  int i,ix,t,x,y,z,beta,px,py,pz,tc,a;
  double pdotx,cpdotx,spdotx;
  complex tr;
  suNf_propagator sp0,sp1,Usp,spf,sptmp1,sptmp2,sptmp3,spdag;

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
                  _propagator_dagger(spdag,sp1);

		  if(measure_channels[_pipig_ff_re]){
		    //Sequential source
		    //_g0_propagator(sptmp1, sp0); _g5_propagator(sptmp2, sptmp1);//g5 g0 G(y,0)
 		    _g5g0_propagator(sptmp2, sp0);
		    _propagator_mul(sptmp1, sptmp2, spdag); // g5 g0 G(x,0) S(y,0)^dagger
		    _propagator_trace(tr,sptmp1);
		    corr[_pipig_ff_re][corr_ind(px,py,pz,n_mom,tc,nm)] += +0.5*(cpdotx*tr.re + spdotx*tr.im);
		    corr[_pipig_ff_im][corr_ind(px,py,pz,n_mom,tc,nm)] += +0.5*(spdotx*tr.re - cpdotx*tr.im);

		    //2nd contraction
	            //_g0_propagator(sptmp2, sp1); _g5_propagator(sptmp1, sptmp2);//g5 g0 S(y,0)
 		    _g5g0_propagator(sptmp1, sp1);
                    _propagator_dagger(sptmp2,sp0); //G(y,0)^dagger
		    _propagator_mul(spf, sptmp2, sptmp1); // G(y,0)^dagger g5 g0 S(y,0)
		    _propagator_trace(tr,spf);
		    corr[_pipig_ff_re][corr_ind(px,py,pz,n_mom,tc,nm)] += -0.5*(cpdotx*tr.re + spdotx*tr.im);
		    corr[_pipig_ff_im][corr_ind(px,py,pz,n_mom,tc,nm)] += -0.5*(spdotx*tr.re - cpdotx*tr.im);
		  }
		  if (measure_channels[_pipig_conserved_ff]){
		    suNf *u1 = _4FIELD_AT(u_gauge_f,ix,0);
		    int ixmu = iup(ix,0);
		    int a;

		    // Tr [ 1/2 S^(0,y+\mu) g_5(1+g_0) U^(y) S(y,x) 
		    _suNf_inverse_prop_multiply(Usp,*u1,sp0);
		    sptmp1=Usp;
		    _g0_propagator(sptmp2,sptmp1);
		    _propagator_add(sptmp1,sptmp1,sptmp2);
		    _g5_propagator(sptmp2,sptmp1);
		    for (a=0;a<NF;a++) for (beta=0;beta<4;beta++){ 
			_propagator_assign(spf, *_FIELD_AT(&psi1[4*a*nm+beta*nm+i],ixmu), a,beta); 
		    }
		    //_propagator_muldag_trace(tr,sptmp2,spf);
		    _propagator_dagger(sptmp1,spf);
		    _propagator_mul(sptmp3,sptmp1,sptmp2);	
		    //_propagator_trace(tr,sptmp1);		  
		    //corr[_pipig_conserved_ff][corr_ind(px,py,pz,n_mom,tc,nm)] += +0.5*(cpdotx*tr.re+spdotx*tr.im);

		    // Tr [ -1/2 S^(0,y) g_5(1-g_0) U(y) S(y+mu,x) 
		    for (a=0;a<NF;a++) for (beta=0;beta<4;beta++){ 
			_propagator_assign(spf, *_FIELD_AT(&psi0[4*a*nm+beta*nm+i],ixmu), a,beta); 
		      }
		    _suNf_prop_multiply(Usp,*u1,spf);
		    sptmp1=Usp;
		    _g0_propagator(sptmp2,sptmp1);
		    _propagator_sub(sptmp1,sptmp1,sptmp2);
		    _g5_propagator(sptmp2,sptmp1);
		    //_propagator_muldag_trace(tr,sptmp2,sp1);		
		    _propagator_mul(sptmp1,spdag,sptmp2);	
		    //_propagator_trace(tr,sptmp1);    
		    //corr[_pipig_conserved_ff][corr_ind(px,py,pz,n_mom,tc,nm)] += -0.5*(cpdotx*tr.re+spdotx*tr.im);

		    _propagator_sub(sptmp2,sptmp3,sptmp1);
		    _propagator_trace(tr,sptmp2);    
		    corr[_pipig_conserved_ff][corr_ind(px,py,pz,n_mom,tc,nm)] += 0.5*(cpdotx*tr.re + spdotx*tr.im);
		  }
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

void measure_discon(spinor_field* psi0, spinor_field* eta, int nm, int tau){
  /*int i;
  int mtmp[NCHANNELS];

  for (i=0;i<NCHANNELS;++i){
    mtmp[i] = measure_channels[i];
    measure_channels[i]=0;
  }*/
  init_corrs(nm,1);
  measure_channels[_g5_disc]=1;
  measure_channels[_id_disc]=1;
  lprintf("measure_mesons",50,"measure default mesons");
  measure_mesons_core(psi0, psi0, eta, nm, tau, 1, 0,GLB_T);
  /*for (i=0;i<NCHANNELS;++i){
    measure_channels[i]=mtmp[i];
  }*/
}

void measure_discon_noise(spinor_field* psi0, spinor_field* eta, int nm, int tau){
  /*int i;
  int mtmp[NCHANNELS];

  for (i=0;i<NCHANNELS;++i){
    mtmp[i] = measure_channels[i];
    measure_channels[i]=0;
  }*/
  init_corrs(nm,1);
  lprintf("measure_mesons",50,"measure default mesons");
  measure_discon_noise_core(psi0, psi0, eta, nm, 1, 0, GLB_T);
  /*for (i=0;i<NCHANNELS;++i){
    measure_channels[i]=mtmp[i];
  }*/
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
  for (i=0;i<NCHANNELS;++i){
    measure_channels[i]=0;
  }
  init_corrs(nm,n_mom);

  lprintf("MEASURE_FORMFACTORS",50,"psi0 = sequential prop, psi1 = point prop");

  measure_channels[_g5]=1;
  for(i=0;i<NF;i++){
	  measure_mesons_core(psi1 + 4*i, psi1 + 4*i, eta + 4*i, nm, ti, n_mom, 0, GLB_T);
  }
  measure_channels[_pipig_ff_re]=1;
  measure_channels[_pipig_ff_im]=1;
  measure_channels[_pipig_conserved_ff]=1;
  measure_formfactor_core(psi0, psi1, eta, nm, ti, tf, n_mom, 0, GLB_T);  
}

void measure_formfactors_ext(spinor_field* psi0, spinor_field* psi1, spinor_field* eta, int nm, int ti, int tf, int n_mom, int begin){
  int i;
  for (i=0;i<NCHANNELS;++i){
    measure_channels[i]=0;
  }
  init_corrs(3*nm,n_mom);

  lprintf("MEASURE_FORMFACTORS",50,"psi0 = sequential prop, psi1 = point prop");

  measure_channels[_g5]=1;
  for(i=0;i<NF;i++){
	  measure_mesons_core(psi1 + 4*i, psi1 + 4*i, eta + 4*i, nm, ti, n_mom, GLB_T*begin, 3*GLB_T);
  }
  measure_channels[_pipig_ff_re]=1;
  measure_channels[_pipig_ff_im]=1;
  measure_channels[_pipig_conserved_ff]=1;
  measure_formfactor_core(psi0, psi1, eta, nm, ti, tf, n_mom, GLB_T*begin, 3*GLB_T);  
}

static void print_corr_core(int channel,int lt,int conf, int nm, double* mass, char* label, int n_mom){
  int i,t,px,py,pz;
  for (px=0;px<n_mom;++px) for (py=0;py<n_mom;++py) for (pz=0;pz<n_mom;++pz){ 
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
	}
      }
}
void print_discon_core(int channel,int lt,int conf, int nm, double* mass, char* label, int n_mom){

  int i,t,tau,tm,px,py,pz;
  double *Dcorr = (double*) malloc(sizeof(double)*lt); 

  for (px=0;px<n_mom;++px) for (py=0;py<n_mom;++py) for (pz=0;pz<n_mom;++pz){ 
	  for(i=0; i<nm; i++) {

	    for(t=0;t<lt;++t){ 
	        Dcorr[t] = 0;
		for(tau=0;tau<lt;++tau){ 
       		 tm = (t+tau)%lt;
		 Dcorr[t] += corr[channel][corr_ind(px,py,pz,n_mom,tau+i*lt,nm)] * corr[channel][corr_ind(px,py,pz,n_mom,tm+i*lt,nm)];
	        }
            }
	    lprintf("MAIN",0,"conf #%d mass=%2.6f %s %s %s= ",conf,mass[i],label,channel_types[channel],channel_names[channel]); 
            
	    for(t=0;t<lt;++t) lprintf("MAIN",0,"%e ",Dcorr[t]);
	    lprintf("MAIN",0,"\n"); 
	    fflush(stdout); 
	}
      }

  free(Dcorr);
}

static void print_corr(int lt,int conf, int nm, double* mass, char* label, int n_mom){
  int i;
  char str1[256]; sprintf(str1, "LCORR %s", label); //Modify label
  for(i = 0; i<NCHANNELS; ++i){
    if(measure_channels[i]){

  	if( i == _id_disc || i == _g5_disc ){ 
		print_corr_core(i,lt,conf,nm,mass,str1,n_mom); }
	else{ print_corr_core(i,lt,conf,nm,mass,label,n_mom); }
      
    }
  }
  char str[256]; sprintf(str, "DCORR %s", label); //Modify label
  if(measure_channels[_id_disc]){ print_discon_core(_id_disc,lt,conf,nm,mass,str,n_mom); }
  if(measure_channels[_g5_disc]){ print_discon_core(_g5_disc,lt,conf,nm,mass,str,n_mom); }
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
	corr[k][i] *= -(1./GLB_VOL3);
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
    for(k=0; k<NCHANNELS; k++) {
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


