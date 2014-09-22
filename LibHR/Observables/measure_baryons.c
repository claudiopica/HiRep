/***************************************************************************\
* Copyright (c) 2014 Vincent Drach, Ari Hietanen                            *
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


// baryons 1 ; C[beta][delta]  = S^ii' StildeT^jj' S^ kk'  
void _propagator_baryon1_mul2(complex C[4][4],suNf_propagator S,suNf_propagator Stilde,int i,int ip,int j,int jp,int k,int kp){
	int  alpha, beta, gamma; 
	complex tmp[4][4]; 


	for(alpha=0;alpha<4;++alpha){ 
		for(beta=0;beta<4;++beta){	
			_complex_0( tmp[alpha][beta] );	
			for(gamma=0;gamma<4;gamma++){ 
				_complex_mul_assign( tmp[alpha][beta], S.c[i].c[alpha].c[gamma].c[ip], Stilde.c[j].c[beta].c[gamma].c[jp] ); 
			}
		}
	} 


	for(alpha=0;alpha<4;++alpha){ 
		for(beta=0;beta<4;++beta){ 
			_complex_0( C[alpha][beta] ); 
			for(gamma=0;gamma<4;gamma++){ 
				_complex_mul_assign( C[alpha][beta], tmp[alpha][gamma], S.c[k].c[gamma].c[beta].c[kp] ); 
			}
		}
	}

}

//baryons 2 ; C[beta][delta]  = S^ii'_*tr(_S^ jj' StilteT^kk'  
void _propagator_baryon2_mul2(complex C[4][4],suNf_propagator S,suNf_propagator Stilde,int i,int ip,int j,int jp,int k,int kp){
	int  alpha, beta, gamma; 
	complex tmp; 
	_complex_0(tmp);
	for(alpha=0;alpha<4;++alpha){ 
		for(gamma=0;gamma<4;gamma++){ 
			_complex_mul_assign( tmp, S.c[j].c[alpha].c[gamma].c[jp], Stilde.c[k].c[alpha].c[gamma].c[kp] ); 
		}} 
	for(alpha=0;alpha<4;++alpha){ 
		for(beta=0;beta<4;++beta){ 
			_complex_mul( C[alpha][beta], S.c[i].c[alpha].c[beta].c[ip], tmp ); 
		}} 
}



void contract_baryons(spinor_field* psi0, int tau){

	lprintf("contract_baryons",50,"Performing baryon contraction ");
	int i,j,ix,t,x,y,z,beta,tc;
	int a,b,c,ap,bp,cp;
	suNf_propagator S,Stilde,Stmp;
	int idx;
	complex static C1[4][4]; // to initialize
	complex C2[4][4]; // to initialize
	complex corr_nuc[GLB_T][4][4]; 
	double corr_nuc_re[4][4][GLB_T]; 
	double corr_nuc_im[4][4][GLB_T]; 
	double col_factor[NF][NF][NF][NF][NF][NF];

	double eps[NF][NF][NF];
//	double tmp;
	// ini corrs
	error(NG!=3,1,"contract_baryons [measure_baryons.c]" ,
	"That code does not work for that number of color !\n");
	
#if defined(REPR_FUNDAMENTAL) && NG==3 
	for (a=0;a<NF;++a){
		for (b=0;b<NF;++b){
			for (c=0;c<NF;++c){
				eps[a][b][c]=0;
			}}}

	eps[0][1][2]=1.;
	eps[1][2][0]=1.;
	eps[2][0][1]=1.;
	eps[1][0][2]=-1.;
	eps[0][2][1]=-1.;
	eps[2][1][0]=-1.;

#elif REPR_SYMMETRIC 
	// def of the color factor for the sextet 
	double eS[NF][NG][NG];
	int k,l,m;
	for ( i=0;i<NF;++i){
		for( j=0; j<NG; ++j){
			for( k=0; k<NG; ++k){
				eS[i][j][k]=0.;
			}}}
	

	eS[0][0][0]=1.;
	eS[1][0][1]=1./sqrt(2);
	eS[1][1][0]=1./sqrt(2);
	eS[2][1][1]=1.;
	eS[3][2][0]=1./sqrt(2);
	eS[3][0][2]=1./sqrt(2);
	eS[4][1][2]=1./sqrt(2);
	eS[4][2][1]=1./sqrt(2);
	eS[5][2][2]=1.;


	double tr_eS[6];
	double tmp1,tmp2,tmp3,tmp4,tmp5;
	double tmpmat1[3][3];
	double tmpmat2[3][3];

	for(k=0;k<NF;++k)
	{ 
		tr_eS[k] = eS[k][0][0] +  eS[k][1][1] +  eS[k][2][2];
	}

	for (a=0;a<NF;++a){
		for (b=0;b<NF;++b){
			for (c=0;c<NF;++c){

				tmp1 = tr_eS[a]*tr_eS[b]*tr_eS[c]			;	

				if (b==c) tmp2 = -tr_eS[a];
				else 	tmp2 =0;
				if (a==c) tmp3 =-tr_eS[b];
				else 	tmp3 =0;
				if (a==b) tmp4 = -tr_eS[c];
				else 	tmp4 =0;

				for (k=0;k<3;++k){
					for (l=0;l<3;++l){
						tmpmat1[k][l]=0.;
						for (m=0;m<3;++m){
							tmpmat1[k][l]+= eS[b][k][m]*eS[c][m][l] +  eS[c][k][m]*eS[b][m][l];
						}
					}
				}
				for (k=0;k<3;++k){
					for (l=0;l<3;++l){
						tmpmat2[k][l]=0.;
						for (m=0;m<3;++m){
							tmpmat2[k][l]+= eS[a][k][m]*tmpmat1[m][l] ;
						}
					}
				}
				tmp5 = tmpmat2[0][0] +  tmpmat2[1][1] + tmpmat2[2][2];
				eps[a][b][c] = tmp1+tmp2+tmp3+tmp4+tmp5;
				if(fabs(eps[a][b][c]) < 1e-15) eps[a][b][c]=0; 
				lprintf("COLOR_FACTOR",0," %i %i %i %3.10e \n",a,b,c,eps[a][b][c]);
			}
		}
	}

#else 
error(1,1,"contract_baryons [measure_mesons.c]" ,
	"That code does not work for that representation \n");
#endif



	for (a=0;a<NF;++a){
		for (ap=0;ap<NF;++ap){
			for (b=0;b<NF;++b){
				for (bp=0;bp<NF;++bp){
					for (c=0;c<NF;++c){
						for (cp=0;cp<NF;++cp){
							 col_factor[a][b][c][ap][bp][cp] = eps[a][b][c]*eps[ap][bp][cp];
						}}}}}}



	for (t=0; t<GLB_T; t++) {	 
		for (i=0;i<4;i++){
			for (j=0;j<4;j++){
				_complex_0(corr_nuc[t][i][j])	;
				corr_nuc_re[i][j][t]=0;
				corr_nuc_im[i][j][t]=0;
			}}}																


	// LOOP VOLUME 

	for (t=0; t<T; t++) {	 
		tc = (zerocoord[0]+t+GLB_T-tau)%GLB_T;
		for (x=0; x<X; x++) for (y=0; y<Y; y++) for (z=0; z<Z; z++) { 
			ix=ipt(t,x,y,z);					

			//get S^ab_beta_gamma
			for (a=0;a<NF;++a){
				for (beta=0;beta<4;beta++){
					idx = beta + a*4;  
					_propagator_assign(Stmp, *_FIELD_AT(&psi0[idx],ix),a,beta);
				}
			}
			_propagator_transpose(S,Stmp);
			// Stilde= Gamma S Gamma_tilde
			// C  = i g0 g2 
			// C g5 = i g1 g3
			// nuc case : Stilde = - Cg5 S Cg5 = g5g0g2 S g5g0g2 

			_g5g0g2_propagator(Stmp,S);
			_propagator_g5g0g2(Stilde,Stmp);
			// add delta ?

			for (a=0;a<NF;++a){
				for (ap=0;ap<NF;++ap){
					for (b=0;b<NF;++b){
						for (bp=0;bp<NF;++bp){
							for (c=0;c<NF;++c){
								for (cp=0;cp<NF;++cp){
									if ( fabs(col_factor[a][b][c][ap][bp][cp]) > 1e-5 )
									{							
										//products
										_propagator_baryon1_mul2(C1,S,Stilde,c,cp,b,bp,a,ap); // nuc term 1
										_propagator_baryon2_mul2(C2,S,Stilde,c,ap,a,cp,b,bp); // nuc term 2

										// multiply by color factor and accumulate in the correlator
										for (i=0;i<4;i++){
											for (j=0;j<4;j++){
												// this assumes that the col_factor is real !
												corr_nuc[tc][i][j].re -= col_factor[a][b][c][ap][bp][cp]*(C1[i][j].re - C2[i][j].re);
												corr_nuc[tc][i][j].im -= col_factor[a][b][c][ap][bp][cp]*(C1[i][j].im - C2[i][j].im);

											}
										}
									} // if col_factor != 0
								} // cp
							} //c
						}// bp
					} //b
				} //ap
			}//a
		} //END SPATIAL LOOP
	} //END TIME LOOP

	// split real and im part
	for (t=0; t<T; t++) {	 
		tc = (zerocoord[0]+t+GLB_T-tau)%GLB_T;
		for (i=0;i<4;i++){
			for (j=0;j<4;j++){
				corr_nuc_re[i][j][tc] = corr_nuc[tc][i][j].re;
				corr_nuc_im[i][j][tc] = corr_nuc[tc][i][j].im;

			}}
	} 

	// global sum
	for (i=0;i<4;i++){
		for (j=0;j<4;j++){
			global_sum(corr_nuc_re[i][j],GLB_T);
			global_sum(corr_nuc_im[i][j],GLB_T);
		}}	

	// print output 
	for (t=0;t<GLB_T;t++){ 
		for (i=0;i<4;i++){
			lprintf("CORR_N",0,"%i %3.10e %3.10e  %3.10e %3.10e  %3.10e %3.10e  %3.10e %3.10e \n",t,corr_nuc_re[i][0][t],corr_nuc_im[i][0][t],corr_nuc_re[i][1][t],corr_nuc_im[i][1][t],corr_nuc_re[i][2][t],corr_nuc_im[i][2][t],corr_nuc_re[i][3][t],corr_nuc_im[i][3][t]);
		}}

	lprintf("contract_baryons",50,"Measuring DONE! ");
}


