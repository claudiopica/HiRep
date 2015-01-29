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

enum { _Sin=0, _Sout, _g5r, _idr, _g0r, _g1r, _g2r, _g3r,  _g0g1r, _g0g2r, _g0g3r, _g5g0r, _g5g1r, _g5g2r, _g5g3r, _g5g0g1r, _g5g0g2r, _g5g0g3r, _cg0r, _cg1r, _cg2r, _cg3r, NCHANNELSR };

char* tr_channel_names[NCHANNELSR]={"Sin", "Sout", "g5","id","g0","g1","g2","g3","g0g1","g0g2","g0g3","g5g0","g5g1","g5g2","g5g3","g5g0g1","g5g0g2","g5g0g3","cg0","cg1","cg2","cg3"};

int tr_measure_channels[NCHANNELSR]={1,1, 1,1,1,1, 1,1,1,1, 1,1,1,1, 1,1,1,1, 1,1,1,1};

static suNf_propagator* tr_corr[NCHANNELSR];

static int init = 0;

static void measure_renormalization_core(spinor_field* psi_in, spinor_field* psi_out, int nm, int pt_in, int px_in, int py_in, int pz_in, int pt_out, int px_out, int py_out, int pz_out){
  int i, ix,t,x,y,z,a,beta;
  double pinx, poutx;
  complex eipinx, eipoutx;
  suNf_propagator Sin,Sout,Sout_dag,Stmp1,Stmp2,Stmp3, Stmp4, Sf;

  for(i=0; i<nm; i++) {
    for (t=0; t<T; t++) for (x=0; x<X; x++) for (y=0; y<Y; y++) for (z=0; z<Z; z++) { 
	  ix=ipt(t,x,y,z);					
	  pinx = 2.0*PI*(((double) pt_in)*(t+zerocoord[0])/GLB_T + ((double) px_in)*(x+zerocoord[1])/GLB_X + ((double) py_in)*(y+zerocoord[2])/GLB_Y + ((double) pz_in)*(z+zerocoord[3])/GLB_Z);
	  poutx = 2.0*PI*(((double) pt_out)*(t+zerocoord[0])/GLB_T + ((double) px_out)*(x+zerocoord[1])/GLB_X + ((double) py_out)*(y+zerocoord[2])/GLB_Y + ((double) pz_out)*(z+zerocoord[3])/GLB_Z);

	  eipinx.re = cos(pinx); eipinx.im = -sin(pinx);
	  eipoutx.re = cos(poutx); eipoutx.im = -sin(poutx);

	  for (a=0;a<NF;++a){
	    for (beta=0;beta<4;beta++){ 
	      _propagator_assign(Sin, *_FIELD_AT(&psi_in[a*4+beta],ix),a,beta);
	      _propagator_assign(Sout, *_FIELD_AT(&psi_out[a*4+beta],ix),a,beta);
	    }
	  }

	  _propagator_mulc_assign(Sin, eipinx);
	  _propagator_mulc_assign(Sout, eipoutx);
	  _propagator_dagger(Sout_dag,Sout);

	  _propagator_add(tr_corr[_Sin][i],tr_corr[_Sin][i],Sin); 
	  _propagator_add(tr_corr[_Sout][i],tr_corr[_Sout][i],Sout); 

	  //g5
          _g5_propagator(Stmp1, Sout_dag);
	  _propagator_mul(Stmp2, Stmp1, Sin);
	  _propagator_add(tr_corr[_g5r][i],tr_corr[_g5r][i],Stmp2); 

	  //id
          _g5_propagator(Stmp1, Sout_dag);
          _g5_propagator(Stmp2, Sin);
	  _propagator_mul(Stmp3, Stmp1, Stmp2);
	  _propagator_add(tr_corr[_idr][i],tr_corr[_idr][i],Stmp3); 

	  //g0
          _g5_propagator(Stmp1, Sout_dag);
          _g5g0_propagator(Stmp2, Sin);
	  _propagator_mul(Stmp3, Stmp1, Stmp2);
	  _propagator_add(tr_corr[_g0r][i],tr_corr[_g0r][i],Stmp3);
	  //g1
          _g5_propagator(Stmp1, Sout_dag);
          _g5g1_propagator(Stmp2, Sin);
	  _propagator_mul(Stmp3, Stmp1, Stmp2);
	  _propagator_add(tr_corr[_g1r][i],tr_corr[_g1r][i],Stmp3);
	  //g2
          _g5_propagator(Stmp1, Sout_dag);
          _g5g2_propagator(Stmp2, Sin);
	  _propagator_mul(Stmp3, Stmp1, Stmp2);
	  _propagator_add(tr_corr[_g2r][i],tr_corr[_g2r][i],Stmp3); 
	  //g3
          _g5_propagator(Stmp1, Sout_dag);
          _g5g3_propagator(Stmp2, Sin);
	  _propagator_mul(Stmp3, Stmp1, Stmp2);
	  _propagator_add(tr_corr[_g3r][i],tr_corr[_g3r][i],Stmp3); 

	  //g5g0
          _g5_propagator(Stmp1, Sout_dag);
          _g0_propagator(Stmp2, Sin);
	  _propagator_mul(Stmp3, Stmp1, Stmp2);
	  _propagator_add(tr_corr[_g5g0r][i],tr_corr[_g5g0r][i],Stmp3);
	  //g5g1
          _g5_propagator(Stmp1, Sout_dag);
          _g1_propagator(Stmp2, Sin);
	  _propagator_mul(Stmp3, Stmp1, Stmp2);
	  _propagator_add(tr_corr[_g5g1r][i],tr_corr[_g5g1r][i],Stmp3);
	  //g5g2
          _g5_propagator(Stmp1, Sout_dag);
          _g2_propagator(Stmp2, Sin);
	  _propagator_mul(Stmp3, Stmp1, Stmp2);
	  _propagator_add(tr_corr[_g5g2r][i],tr_corr[_g5g2r][i],Stmp3); 
	  //g5g3
          _g5_propagator(Stmp1, Sout_dag);
          _g3_propagator(Stmp2, Sin);
	  _propagator_mul(Stmp3, Stmp1, Stmp2);
	  _propagator_add(tr_corr[_g5g3r][i],tr_corr[_g5g3r][i],Stmp3); 

	  //g0g1
          _g5_propagator(Stmp1, Sout_dag);
          _g5g0g1_propagator(Stmp2, Sin);
	  _propagator_mul(Stmp3, Stmp1, Stmp2);
	  _propagator_add(tr_corr[_g0g1r][i],tr_corr[_g0g1r][i],Stmp3);
	  //g0g2
          _g5_propagator(Stmp1, Sout_dag);
          _g5g0g2_propagator(Stmp2, Sin);
	  _propagator_mul(Stmp3, Stmp1, Stmp2);
	  _propagator_add(tr_corr[_g0g2r][i],tr_corr[_g0g2r][i],Stmp3); 
	  //g0g3
          _g5_propagator(Stmp1, Sout_dag);
          _g5g0g3_propagator(Stmp2, Sin);
	  _propagator_mul(Stmp3, Stmp1, Stmp2);
	  _propagator_add(tr_corr[_g0g3r][i],tr_corr[_g0g3r][i],Stmp3); 

	  //g5g0g1
          _g5_propagator(Stmp1, Sout_dag);
          _g0g1_propagator(Stmp2, Sin);
	  _propagator_mul(Stmp3, Stmp1, Stmp2);
	  _propagator_add(tr_corr[_g5g0g1r][i],tr_corr[_g5g0g1r][i],Stmp3); 
	  //g5g0g2
          _g5_propagator(Stmp1, Sout_dag);
          _g0g2_propagator(Stmp2, Sin);
	  _propagator_mul(Stmp3, Stmp1, Stmp2);
	  _propagator_add(tr_corr[_g5g0g2r][i],tr_corr[_g5g0g2r][i],Stmp3); 
	  //g5g0g3
          _g5_propagator(Stmp1, Sout_dag);
          _g0g3_propagator(Stmp2, Sin);
	  _propagator_mul(Stmp3, Stmp1, Stmp2);
	  _propagator_add(tr_corr[_g5g0g3r][i],tr_corr[_g5g0g3r][i],Stmp3); 

//cg0
	  suNf *u1 = _4FIELD_AT(u_gauge_f,ix,0); //U_mu (x)
	  int ixmu = iup(ix,0); //x + mu

          // g5 S^(p, z+mu) g5 (1 + g0) U^(z) S(p,z)
	  _suNf_inverse_prop_multiply(Stmp1,*u1,Sin); //U^(z) S(p,z)
	  Stmp3 = Stmp1; //U^(z) S(p,z)
	  _g0_propagator(Stmp2, Stmp1); // g0 U^(z) S(p,z)
	  _propagator_add(Stmp1,Stmp1,Stmp2);// (1 + g0) U^(z) S(p,z)
	  _g5_propagator(Stmp2,Stmp1); // g5 (1 + g0) U^(z) S(p,z)
	  for (a=0;a<NF;a++) for (beta=0;beta<4;beta++){ 
	  	_propagator_assign(Sf, *_FIELD_AT(&psi_out[4*a+beta],ixmu), a,beta); //S(p,z+mu)
	  }
	  _propagator_dagger(Stmp1,Sf); //S^(p, z+mu)
	  _propagator_mul(Stmp3,Stmp1,Stmp2);	// S^(p, z+mu) g5 (1 + g0) U^(z) S(p,z)
          _g5_propagator(Stmp4, Stmp3); // g5 S^(p, z+mu) g5 (1 + g0) U^(z) S(p,z)

          _propagator_add(tr_corr[_cg0r][i],tr_corr[_cg0r][i],Stmp4); 

          // g5 S^(p, z) g5 (1 + g0) U^(z) S(p,z+mu)
	  for (a=0;a<NF;a++) for (beta=0;beta<4;beta++){ 
	  	_propagator_assign(Sf, *_FIELD_AT(&psi_in[4*a+beta],ixmu), a,beta); //S(p,z+mu)
	  }
	  _suNf_inverse_prop_multiply(Stmp1,*u1,Sf); //U^(z) S(p,z+mu)
	  Stmp3 = Stmp1; //U^(z) S(p,z+mu)
	  _g0_propagator(Stmp2, Stmp1); // g0 U^(z) S(p,z+mu)
	  _propagator_sub(Stmp1,Stmp1,Stmp2);// (1 - g0) U^(z) S(p,z+mu)
	  _g5_propagator(Stmp2,Stmp1); // g5 (1 - g0) U^(z) S(p,z+mu)
	  
	  _propagator_mul(Stmp3,Sout_dag,Stmp2); // S^(p, z) g5 (1 - g0) U^(z) S(p,z+mu)
          _g5_propagator(Stmp4, Stmp3); // g5 S^(p, z) g5 (1 - g0) U^(z) S(p,z+mu)

          _propagator_sub(tr_corr[_cg0r][i],tr_corr[_cg0r][i],Stmp4); 


//cg1
	  u1 = _4FIELD_AT(u_gauge_f,ix,1); //U_mu (x)
	  ixmu = iup(ix,1); //x + mu

          // g5 S^(p, z+mu) g5 (1 + g1) U^(z) S(p,z)
	  _suNf_inverse_prop_multiply(Stmp1,*u1,Sin); //U^(z) S(p,z)
	  Stmp3 = Stmp1; //U^(z) S(p,z)
	  _g1_propagator(Stmp2, Stmp1); // g1 U^(z) S(p,z)
	  _propagator_add(Stmp1,Stmp1,Stmp2);// (1 + g1) U^(z) S(p,z)
	  _g5_propagator(Stmp2,Stmp1); // g5 (1 + g1) U^(z) S(p,z)
	  for (a=0;a<NF;a++) for (beta=0;beta<4;beta++){ 
	  	_propagator_assign(Sf, *_FIELD_AT(&psi_out[4*a+beta],ixmu), a,beta); //S(p,z+mu)
	  }
	  _propagator_dagger(Stmp1,Sf); //S^(p, z+mu)
	  _propagator_mul(Stmp3,Stmp1,Stmp2);	// S^(p, z+mu) g5 (1 + g1) U^(z) S(p,z)
          _g5_propagator(Stmp4, Stmp3); // g5 S^(p, z+mu) g5 (1 + g1) U^(z) S(p,z)

          _propagator_add(tr_corr[_cg1r][i],tr_corr[_cg1r][i],Stmp4); 

          // g5 S^(p, z) g5 (1 + g1) U^(z) S(p,z+mu)
	  for (a=0;a<NF;a++) for (beta=0;beta<4;beta++){ 
	  	_propagator_assign(Sf, *_FIELD_AT(&psi_in[4*a+beta],ixmu), a,beta); //S(p,z+mu)
	  }
	  _suNf_inverse_prop_multiply(Stmp1,*u1,Sf); //U^(z) S(p,z+mu)
	  Stmp3 = Stmp1; //U^(z) S(p,z+mu)
	  _g1_propagator(Stmp2, Stmp1); // g1 U^(z) S(p,z+mu)
	  _propagator_sub(Stmp1,Stmp1,Stmp2);// (1 - g1) U^(z) S(p,z+mu)
	  _g5_propagator(Stmp2,Stmp1); // g5 (1 - g1) U^(z) S(p,z+mu)
	  
	  _propagator_mul(Stmp3,Sout_dag,Stmp2); // S^(p, z) g5 (1 - g1) U^(z) S(p,z+mu)
          _g5_propagator(Stmp4, Stmp3); // g5 S^(p, z) g5 (1 - g1) U^(z) S(p,z+mu)

          _propagator_sub(tr_corr[_cg1r][i],tr_corr[_cg1r][i],Stmp4); 



//cg2
	  u1 = _4FIELD_AT(u_gauge_f,ix,2); //U_mu (x)
	  ixmu = iup(ix,2); //x + mu

          // g5 S^(p, z+mu) g5 (1 + g2) U^(z) S(p,z)
	  _suNf_inverse_prop_multiply(Stmp1,*u1,Sin); //U^(z) S(p,z)
	  Stmp3 = Stmp1; //U^(z) S(p,z)
	  _g2_propagator(Stmp2, Stmp1); // g2 U^(z) S(p,z)
	  _propagator_add(Stmp1,Stmp1,Stmp2);// (1 + g2) U^(z) S(p,z)
	  _g5_propagator(Stmp2,Stmp1); // g5 (1 + g2) U^(z) S(p,z)
	  for (a=0;a<NF;a++) for (beta=0;beta<4;beta++){ 
	  	_propagator_assign(Sf, *_FIELD_AT(&psi_out[4*a+beta],ixmu), a,beta); //S(p,z+mu)
	  }
	  _propagator_dagger(Stmp1,Sf); //S^(p, z+mu)
	  _propagator_mul(Stmp3,Stmp1,Stmp2);	// S^(p, z+mu) g5 (1 + g2) U^(z) S(p,z)
          _g5_propagator(Stmp4, Stmp3); // g5 S^(p, z+mu) g5 (1 + g2) U^(z) S(p,z)

          _propagator_add(tr_corr[_cg2r][i],tr_corr[_cg2r][i],Stmp4); 

          // g5 S^(p, z) g5 (1 + g2) U^(z) S(p,z+mu)
	  for (a=0;a<NF;a++) for (beta=0;beta<4;beta++){ 
	  	_propagator_assign(Sf, *_FIELD_AT(&psi_in[4*a+beta],ixmu), a,beta); //S(p,z+mu)
	  }
	  _suNf_inverse_prop_multiply(Stmp1,*u1,Sf); //U^(z) S(p,z+mu)
	  Stmp3 = Stmp1; //U^(z) S(p,z+mu)
	  _g2_propagator(Stmp2, Stmp1); // g2 U^(z) S(p,z+mu)
	  _propagator_sub(Stmp1,Stmp1,Stmp2);// (1 - g2) U^(z) S(p,z+mu)
	  _g5_propagator(Stmp2,Stmp1); // g5 (1 - g2) U^(z) S(p,z+mu)
	  
	  _propagator_mul(Stmp3,Sout_dag,Stmp2); // S^(p, z) g5 (1 - g2) U^(z) S(p,z+mu)
          _g5_propagator(Stmp4, Stmp3); // g5 S^(p, z) g5 (1 - g2) U^(z) S(p,z+mu)

          _propagator_sub(tr_corr[_cg2r][i],tr_corr[_cg2r][i],Stmp4); 




//cg3
	  u1 = _4FIELD_AT(u_gauge_f,ix,3); //U_mu (x)
	  ixmu = iup(ix,3); //x + mu

          // g5 S^(p, z+mu) g5 (1 + g3) U^(z) S(p,z)
	  _suNf_inverse_prop_multiply(Stmp1,*u1,Sin); //U^(z) S(p,z)
	  Stmp3 = Stmp1; //U^(z) S(p,z)
	  _g3_propagator(Stmp2, Stmp1); // g3 U^(z) S(p,z)
	  _propagator_add(Stmp1,Stmp1,Stmp2);// (1 + g3) U^(z) S(p,z)
	  _g5_propagator(Stmp2,Stmp1); // g5 (1 + g3) U^(z) S(p,z)
	  for (a=0;a<NF;a++) for (beta=0;beta<4;beta++){ 
	  	_propagator_assign(Sf, *_FIELD_AT(&psi_out[4*a+beta],ixmu), a,beta); //S(p,z+mu)
	  }
	  _propagator_dagger(Stmp1,Sf); //S^(p, z+mu)
	  _propagator_mul(Stmp3,Stmp1,Stmp2);	// S^(p, z+mu) g5 (1 + g3) U^(z) S(p,z)
          _g5_propagator(Stmp4, Stmp3); // g5 S^(p, z+mu) g5 (1 + g3) U^(z) S(p,z)

          _propagator_add(tr_corr[_cg3r][i],tr_corr[_cg3r][i],Stmp4); 

          // g5 S^(p, z) g5 (1 + g3) U^(z) S(p,z+mu)
	  for (a=0;a<NF;a++) for (beta=0;beta<4;beta++){ 
	  	_propagator_assign(Sf, *_FIELD_AT(&psi_in[4*a+beta],ixmu), a,beta); //S(p,z+mu)
	  }
	  _suNf_inverse_prop_multiply(Stmp1,*u1,Sf); //U^(z) S(p,z+mu)
	  Stmp3 = Stmp1; //U^(z) S(p,z+mu)
	  _g3_propagator(Stmp2, Stmp1); // g3 U^(z) S(p,z+mu)
	  _propagator_sub(Stmp1,Stmp1,Stmp2);// (1 - g3) U^(z) S(p,z+mu)
	  _g5_propagator(Stmp2,Stmp1); // g5 (1 - g3) U^(z) S(p,z+mu)
	  
	  _propagator_mul(Stmp3,Sout_dag,Stmp2); // S^(p, z) g5 (1 - g3) U^(z) S(p,z+mu)
          _g5_propagator(Stmp4, Stmp3); // g5 S^(p, z) g5 (1 - g3)) U^(z) S(p,z+mu)

          _propagator_sub(tr_corr[_cg3r][i],tr_corr[_cg3r][i],Stmp4); 




	} //END SPACETIME
  } //END NM LOOP
	
      
}


static void init_tr_corrs(int nm){
  int k,i;
  static int first = 1;
  if(first){
	for (k=0;k<NCHANNELSR;++k){
         tr_corr[k]= malloc(sizeof(suNf_propagator)*nm);
        }
	first = 0;
  }
  if (!init){
    for(k=0; k<NCHANNELSR; k++){
	for(i=0;i<nm;i++){
	  _propagator_zero(tr_corr[k][i]);
	}
    }
    init = 1;
  }
}


void measure_renormalization(spinor_field* psi_in, spinor_field* psi_out, int nm, int pt_in, int px_in, int py_in, int pz_in, int pt_out, int px_out, int py_out, int pz_out){
  init_tr_corrs(nm);
  lprintf("measure_renormalization",50,"measure default renormalization");
  measure_renormalization_core(psi_in, psi_out, nm, pt_in, px_in, py_in, pz_in, pt_out, px_out, py_out, pz_out);
}

static void print_renormalization_core(int channel,int conf, int nm, double* mass, char* label, int pt_in, int px_in, int py_in, int pz_in, int pt_out, int px_out, int py_out, int pz_out ){
  int i,a,b,alpha,beta;

  for(i=0; i<nm; i++) { 

    lprintf("MAIN",0,"conf #%d mass=%2.6f %s %s p_in(%d,%d,%d,%d) p_out(%d,%d,%d,%d) =",conf,mass[i],label,tr_channel_names[channel], pt_in, px_in, py_in, pz_in, pt_out, px_out, py_out, pz_out );

    for(a=0;a<NF;a++) for(alpha=0;alpha<4;alpha++) for(b=0;b<NF;b++) for(beta=0;beta<4;beta++){
	lprintf("MAIN",0," ( %1.12g , %1.12g ) ", _PROP_AT(tr_corr[channel][i],a,alpha,beta,b).re, _PROP_AT(tr_corr[channel][i],a,alpha,beta,b).im );
    	}
    lprintf("MAIN",0,"\n");
     
    fflush(stdout);
	
  }
}

void print_renormalization(int conf, int nm, double* mass, char* label, int pt_in, int px_in, int py_in, int pz_in, int pt_out, int px_out, int py_out, int pz_out){
  int k,a,b,alpha,beta,i;

  if (init !=0 ){
   for(i=0; i<nm; i++) { 
    for(k=0; k<NCHANNELSR; k++) {
      for(a=0;a<NF;a++) for(alpha=0;alpha<4;alpha++) for(beta=0;beta<4;beta++) for(b=0;b<NF;b++) {
      	global_sum( &_PROP_AT(tr_corr[k][i],a,alpha,beta,b).re ,1);
      	global_sum( &_PROP_AT(tr_corr[k][i],a,alpha,beta,b).im ,1);
      }
      _propagator_mul_assign( tr_corr[k][i], (1./(double)GLB_VOLUME) );
      print_renormalization_core(k,conf,nm,mass,label,pt_in, px_in, py_in, pz_in, pt_out, px_out, py_out, pz_out);
    }
   }
  }
  init = 0;

}



