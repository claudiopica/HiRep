/*******************************************************************************
 *
 *Contraction routines for scattering length computation
 *VD 2014
 *******************************************************************************/


#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "global.h"
#include "io.h"
#include "random.h"
#include "error.h"
#include "geometry.h"
#include "memory.h"
#include "statistics.h"
#include "update.h"
#include "scattering.h"
#include "observables.h"
#include "suN.h"
#include "suN_types.h"
#include "dirac.h"
#include "linear_algebra.h"
#include "inverters.h"
#include "representation.h"
#include "utils.h"
#include "logger.h"
#include "communications.h"
#include "gamma_spinor.h"
#include "spin_matrix.h"

#define PI 3.141592653589793238462643383279502884197


#include "cinfo.c"

#if defined(ROTATED_SF) && defined(BASIC_SF)
#error This code does not work with the Schroedinger functional !!!
#endif

#ifdef FERMION_THETA
#error This code does not work with the fermion twisting !!!
#endif

//Apply fixed boundaryconditions by zeroing links at time slice tau to direction 0.
static void fix_T_bc(int tau){
	int index;
	int ix,iy,iz;
	suNf *u;
	if (--tau<0) tau+= GLB_T;
	lprintf("meson_measurements",15,"Setting Dirichlet boundary conidtion at global time slice %d, %d\n",tau,T_BORDER);
	if((zerocoord[0]-1<=tau && zerocoord[0]+T>tau) || (zerocoord[0]==0 && tau==GLB_T-1)) { 
		for (ix=0;ix<X_EXT;++ix) for (iy=0;iy<Y_EXT;++iy) for (iz=0;iz<Z_EXT;++iz){
			if( ( (tau==zerocoord[0]-1) || (zerocoord[0]==0 && tau==GLB_T-1)) && (NP_T>1) ){
				index=ipt_ext(0,ix,iy,iz);
			}
			else{
				index=ipt_ext(T_BORDER+tau-zerocoord[0],ix,iy,iz); 
			}
			if(index!=-1) {
				u=pu_gauge_f(index,0);
				_suNf_zero(*u);
			}
		}
	}
	lprintf("meson_measurements",50,"Boundaries set!\n");
}





void measure_pion_scattering(double* m, int nhits,int conf_num, double precision,int ts){
	spinor_field* source_ts = alloc_spinor_field_f(4,&glattice);
	spinor_field* source_tsp1 = alloc_spinor_field_f(4,&glattice);

	spinor_field* prop_ts =  alloc_spinor_field_f(4 ,&glattice);
	spinor_field* prop_tsp1 =  alloc_spinor_field_f(4 ,&glattice);

	int k;
	struct timeval start, end, etime;

	gettimeofday(&start,0);
	init_propagator_eo(1, m, precision);

	for (k=0;k<nhits;k++){
		lprintf("MAIN",0,"In the nhits loops ! \n");
		// to have dirichlet boundary conditions for testing only :
		//		fix_T_bc(ts);

	
		spinor_field_zero_f(source_ts);
		create_diluted_source_equal_atau(source_ts, ts);
		calc_propagator(prop_ts,source_ts,4);

		//Print_spinor_field(&source_ts[0],0,0,0,0);				
		//Print_spinor_field(&source_ts[1],0,0,0,0);				
		//Print_spinor_field(&source_ts[0],0,0,0,1);				
		//Print_spinor_field(&source_ts[0],1,0,0,0);				
		//Print_spinor_field(&source_ts[1],1,0,0,0);				


		spinor_field_zero_f(source_tsp1);
		create_diluted_source_equal_atau(source_tsp1, ts+1);

		//Print_spinor_field(&source_tsp1[0],0,0,0,1);				
		//Print_spinor_field(&source_tsp1[1],0,0,0,1);				
		//Print_spinor_field(&source_tsp1[0],0,0,0,1);				
		//Print_spinor_field(&source_tsp1[0],1,0,0,0);				
		//Print_spinor_field(&source_tsp1[1],1,0,0,0);				

		calc_propagator(prop_tsp1,source_tsp1,4);

		// for checks
		//measure_mesons(prop_ts, source_ts,1, 0);
		//print_mesons(1,10,1,m,"DEFAULT_SEMWALL");
		//measure_mesons(prop_tsp1, source_tsp1,1, 0);
		//print_mesons(1,10,1,m,"DEFAULT_SEMWALL");

		lprintf("MAIN",0,"Start to perform the contractions ... \n");
		contract_pion_scatt(prop_ts,prop_tsp1,k,ts);
		lprintf("MAIN",0,"Contraction done\n");
	}
	free_propagator_eo();

	free_spinor_field_f(source_ts);
	free_spinor_field_f(source_tsp1);
	free_spinor_field_f(prop_ts);
	free_spinor_field_f(prop_tsp1);

	gettimeofday(&end,0);
	timeval_subtract(&etime,&end,&start);
	lprintf("TIMING",0,"Sources generation, invert and contract for %i sources done [%ld sec %ld usec]\n",nhits,etime.tv_sec,etime.tv_usec);

	// debugging or should I code it properly for future project ?

	spinor_field* source2 = alloc_spinor_field_f(1,&glattice);
	spinor_field* prop_ts_1 =  alloc_spinor_field_f(1 ,&glattice);
	spinor_field* prop_tsp1_1 =  alloc_spinor_field_f(1 ,&glattice);
	init_propagator_eo(1, m, precision);

	create_diluted_source_equal_spinorfield1(source2,ts);

	calc_propagator(prop_ts_1,source2,1);// No dilution

	spinor_field_zero_f(source2);

	create_diluted_source_equal_spinorfield1(source2,ts+1);

	calc_propagator(prop_tsp1_1,source2,1);// No dilution

//	contract_pion_scatt_1spinorfield(prop_ts_1,prop_tsp1_1,0,ts);

	free_propagator_eo();
	free_spinor_field_f(source2);
	free_spinor_field_f(prop_tsp1_1);
	free_spinor_field_f(prop_ts_1);


}


void contract_pion_scatt(spinor_field* phi_ts,spinor_field* phi_tsp1,int k,int ts){

	int x,y,z,t,tc;
	int ix,ix_pts_pa,ix_pts;
	int beta,mu,gamma,l,m,a;
	complex tmpc;
	complex *A1,*A2;
	complex *D1,*D2;
	double norm;
	struct timeval start, end, etime;
	suNf_spin_matrix sm1,sm2, sm3,sm4;
	complex* corr[4];
	double* corr_re[4];
	double* corr_im[4];
	complex *B1[4][4],*B2[4][4], *C1[4][4], *C2[4][4];
	complex* tmp2[4];
	double* corr_2pt;
	double* corr_tmp[2];

	// debugging variables
	tmpc.re=0;tmpc.im=0;
	A1 =  malloc(sizeof(complex)*GLB_T);
	A2 =  malloc(sizeof(complex)*GLB_T);
	D1 =  malloc(sizeof(complex)*GLB_T);
	D2 =  malloc(sizeof(complex)*GLB_T);


	for (l=0;l<4;++l)
	{
		corr[l]=(complex*) malloc(sizeof(complex)*GLB_T);
		corr_re[l]=(double*) malloc(sizeof(double)*GLB_T);
		corr_im[l]=(double*) malloc(sizeof(double)*GLB_T);
		for (m=0;m<4;++m)
		{
			B1[l][m]= (complex*) malloc(sizeof(complex)*GLB_T);
			B2[l][m]= (complex*) malloc(sizeof(complex)*GLB_T);
			C1[l][m]= (complex*) malloc(sizeof(complex)*GLB_T);
			C2[l][m]= (complex*) malloc(sizeof(complex)*GLB_T);
		}

	}


	for (l=0;l<2 ;++l)  corr_tmp[l]=(double*) malloc(sizeof(double)*GLB_T);
	corr_2pt=(double*) malloc(sizeof(double)*GLB_T);

	// initilization of the complex corr //innner loop should be the last index of the array
	for (t=0; t<GLB_T; t++) 
	{
		for (l=0;l<4;l++)
		{
			_complex_0(corr[l][t]);
			corr_re[l][t] = 0.;	
			corr_im[l][t] = 0.;	
			if ( l < 2 )  corr_tmp[l][t]=0.;
		}	
		corr_2pt[t]=0.;

		_complex_0(A1[t]);
		_complex_0(A2[t]);
		_complex_0(D1[t]);
		_complex_0(D2[t]);

		for (beta=0;beta<4;beta++) for (gamma=0;gamma<4;gamma++){
			_complex_0(B1[beta][gamma][t]);
			_complex_0(B2[beta][gamma][t]);
 			_complex_0(C1[beta][gamma][t]);
			_complex_0(C2[beta][gamma][t]);
 		}



	}	

	gettimeofday(&start,0);

	for (t=0; t<T; t++)
	{

		//		if(COORD[0]==((t+ts)%GLB_T)/T) {// Check that tau is in this thread.	
		/* what is tc and offset ? MPI related ? */
		tc = (zerocoord[0]+t+GLB_T-ts)%GLB_T;

		/* loop on spatial volume */
		for (x=0; x<X; x++) for (y=0; y<Y; y++) for (z=0; z<Z; z++) {

			ix_pts = ipt(t,x,y,z);
			ix_pts_pa = iup(ix_pts,0);

			//			ix_pts_pa = ipt((t+ts+1)%GLB_T,x,y,z);
			//			ix_pts = ipt((t+ts)%GLB_T,x,y,z);
			//			printf("%i %i \n",t+ts+1,(t+ts+1)%GLB_T);

			// A term
			for (beta=0;beta<4;beta++)
			{
				_spinmatrix_assign_row(sm1, *_FIELD_AT(&phi_tsp1[beta],ix_pts_pa), beta);
				_spinmatrix_assign_row(sm2, *_FIELD_AT(&phi_ts[beta],ix_pts), beta);
			}

			_spinmatrix_mul_trace_assign(A1[tc], sm1, sm1);
			_spinmatrix_mul_trace_assign(A2[tc], sm2, sm2);

			//			if (t==0)  printf("A1.re %3.10e A2.re  %3.10e A1.im %3.10e A2.im  %3.10e\n",A1.re,A2.re,A1.im,A2.im);

			// B term



				for (mu=0;mu<4;mu++)
				{
					_spinmatrix_assign_row(sm1, *_FIELD_AT(&phi_ts[mu],ix_pts_pa), mu); //  mu : source,  beta : internal one. A_mu beta
					_spinmatrix_assign_row(sm2, *_FIELD_AT(&phi_tsp1[mu],ix_pts_pa), mu); // B_mu beta

					_spinmatrix_assign_row(sm3, *_FIELD_AT(&phi_tsp1[mu],ix_pts), mu);
					_spinmatrix_assign_row(sm4, *_FIELD_AT(&phi_ts[mu],ix_pts), mu);
				}

			// compute a new spinmatrix  : (A_mu beta * B_gamma_beta )
			// trace color index -> it becomes a non standard struct. with only spin indices
			// accumulate
			// multiply with second part and finally trace over spin indices. 
			for (mu=0;mu<4;mu++)
			{
				for (gamma=0;gamma<4;gamma++)
				{
          _spinor_prod_assign_f(B1[mu][gamma][tc],sm1.c[mu],sm2.c[gamma]);
          _spinor_prod_assign_f(B2[mu][gamma][tc],sm3.c[mu],sm4.c[gamma]);
				}
			}

			// C term   : useless C^ast = B
			for (beta=0;beta<4;beta++)
			{
				_spinmatrix_assign_row(sm1, *_FIELD_AT(&phi_tsp1[beta],ix_pts_pa), beta);
				_spinmatrix_assign_row(sm2, *_FIELD_AT(&phi_ts[beta],ix_pts_pa), beta);

				_spinmatrix_assign_row(sm3, *_FIELD_AT(&phi_ts[beta],ix_pts), beta);
				_spinmatrix_assign_row(sm4, *_FIELD_AT(&phi_tsp1[beta],ix_pts), beta);

			}

      for (mu=0;mu<4;mu++)
			{
				for (gamma=0;gamma<4;gamma++)
				{
          _spinor_prod_assign_f(C1[mu][gamma][tc],sm1.c[mu],sm2.c[gamma]);
          _spinor_prod_assign_f(C2[mu][gamma][tc],sm3.c[mu],sm4.c[gamma]);
				}
			}


			// D term
			for (beta=0;beta<4;beta++)
			{
				_spinmatrix_assign_row(sm1, *_FIELD_AT(&phi_ts[beta],ix_pts_pa), beta);
				_spinmatrix_assign_row(sm2, *_FIELD_AT(&phi_tsp1[beta],ix_pts), beta);
			}



			_spinmatrix_mul_trace_assign(D1[tc], sm1, sm1);
			_spinmatrix_mul_trace_assign(D2[tc], sm2, sm2);



		}/* end loops on space */

		corr_tmp[0][tc] = A1[tc].re;
		corr_tmp[1][tc] = A2[tc].re;
		//	}

	} /* end loop on time */

	//global sums	
	global_sum((double*)A1,2*GLB_T);
	global_sum((double*)A2,2*GLB_T);
	global_sum((double*)D1,2*GLB_T);
	global_sum((double*)D2,2*GLB_T);
	for (l=0;l<4;l++) for (m=0;m<4;m++) {
		global_sum((double*)B1[l][m],2*GLB_T);
		global_sum((double*)B2[l][m],2*GLB_T);
		global_sum((double*)C1[l][m],2*GLB_T);
		global_sum((double*)C2[l][m],2*GLB_T);
	}
	// build correlators
	for(t=0;t<GLB_T;++t) {
		_complex_mul(corr[0][t],A1[t],A2[t]); 

		for (l=0;l<4;l++) for (m=0;m<4;m++) {
			_complex_mul_assign(corr[1][t],B1[l][m][t],B2[m][l][t]);
			_complex_mul_assign(corr[2][t],C1[l][m][t],C2[m][l][t]);
		}
		_complex_mul(corr[3][t],D1[t],D2[t]);
	}

	// normalize
	for (l=0;l<4;l++) {
		if (l==1 || l ==2) norm=-GLB_VOL3*GLB_VOL3; // B and C term have a global  minus sign  in the contraction
		else	  norm= GLB_VOL3*GLB_VOL3;	 // A and B
		for(t=0;t<GLB_T;++t) {
			corr[l][t].re = corr[l][t].re/norm;
			corr[l][t].im = corr[l][t].im/norm;
		}
	}

	// build two point
	for (t=0;t<T;t++){

		tc = (zerocoord[0]+t+GLB_T-ts)%GLB_T;
		//		if(COORD[0]==((t+ts)%GLB_T)/T) {// Check that tau is in this thread.	
		//		corr_2pt[tc] = 0.5 *(corr_tmp[0][(tc+GLB_T)%GLB_T]+corr_tmp[1][tc])/(GLB_VOL3);
		corr_2pt[tc] =  (corr_tmp[1][tc])/(GLB_VOL3);
		//		}

	}

	global_sum(corr_2pt,GLB_T);

	lprintf("CORR_PP",0," Correlators C_pipi : T l nhits re im ! \n");
	for(l=0;l<4;l++ ) for(t=0;t<GLB_T;++t) lprintf("CORR_PP",0,"%i %i %i %3.10e %3.10e \n",t,l,k,corr[l][t].re,corr[l][t].im);

	lprintf("CORR_P",0," Correlators C_pi : T  nhits re ! \n");
	for(t=0;t<GLB_T;++t)  lprintf("CORR_P",0,"%i %i %3.10e  \n",t,k,corr_2pt[t]);
	fflush(stdout);
	gettimeofday(&end,0);
	timeval_subtract(&etime,&end,&start);
	lprintf("TIMING",0,"Contractions for Pion-Pion ( nhits = %i ) scattering done in  [%ld sec %ld usec]\n",k,etime.tv_sec,etime.tv_usec);

}

void contract_pion_scatt_1spinorfield(spinor_field* phi_ts,spinor_field* phi_tsp1,int k,int ts){

	int x,y,z,t,tc;
	int ix,ix_pts_pa,ix_pts;
	int beta,mu,gamma,l,a;
	complex A1,A2,B1,B2,tmpc;
	complex C1,C2,D1,D2;

	double norm;
	struct timeval start, end, etime;
	suNf_spin_matrix sm1,sm2, sm3,sm4;
	complex* corr[4];
	double* corr_re[4];
	double* corr_im[4];
	complex* tmp1[4];
	complex* tmp2[4];
	double* corr_2pt;
	double* corr_tmp[2];

	// debugging variables
	tmpc.re=0;tmpc.im=0;

	for (l=0;l<4;++l)
	{
		corr[l]=(complex*) malloc(sizeof(complex)*GLB_T);
		corr_re[l]=(double*) malloc(sizeof(double)*GLB_T);
		corr_im[l]=(double*) malloc(sizeof(double)*GLB_T);

		tmp1[l]=(complex*) malloc(sizeof(complex)*4);
		tmp2[l]=(complex*) malloc(sizeof(complex)*4);

	}

	for (l=0;l<2 ;++l)  corr_tmp[l]=(double*) malloc(sizeof(double)*GLB_T);
	corr_2pt=(double*) malloc(sizeof(double)*GLB_T);
	// initilization of the complex corr
	for (t=0; t<GLB_T; t++) 
	{
		for (l=0;l<4;l++)
		{
			_complex_0(corr[l][t]);
			corr_re[l][t] = 0.;	
			corr_im[l][t] = 0.;	
			if ( l < 2 )  corr_tmp[l][t]=0.;
		}	
		corr_2pt[t]=0.;  
	}	


	gettimeofday(&start,0);

	for (t=0; t<T; t++)
	{
		//	  if(COORD[0]==((t+ts)%GLB_T)/T) {// Check that t + ts is in this thread.	

		A1.re=0;A1.im=0.;A2.re=0;A2.im=0.;
		B1.re=0;B1.im=0.;B2.re=0;B2.im=0.;
		C1.re=0;C1.im=0.;C2.re=0;C2.im=0.;
		D1.re=0;D1.im=0.;D2.re=0;D2.im=0.;


		for (beta=0;beta<4;beta++)
		{
			for (gamma=0;gamma<4;gamma++)
			{
				_complex_0(tmp1[beta][gamma]);
				_complex_0(tmp2[beta][gamma]);
			}
		}

		/* what is tc and offset ? MPI related ? */
		tc = (zerocoord[0]+t+GLB_T-ts)%GLB_T;

		/* loop on spatial volume */
		for (x=0; x<X; x++) for (y=0; y<Y; y++) for (z=0; z<Z; z++) {

			ix_pts_pa = ipt((t+ts+1)%GLB_T,x,y,z);
			ix_pts = ipt((t+ts)%GLB_T,x,y,z);
			// NOTE THAT IT PROBABLY WONT WORK WITH MPI (at least if the number of proc in the time direction is different from 0) !!

			// A term
			_spinor_prod_assign_g(A1, *_FIELD_AT(phi_tsp1,ix_pts_pa), *_FIELD_AT(phi_tsp1,ix_pts_pa));
			_spinor_prod_assign_g(A2, *_FIELD_AT(phi_ts,ix_pts), *_FIELD_AT(phi_ts,ix_pts));



			// B term

			_spinor_prod_assign_g(B1, *_FIELD_AT(phi_ts,ix_pts_pa), *_FIELD_AT(phi_tsp1,ix_pts_pa));
			_spinor_prod_assign_g(B2, *_FIELD_AT(phi_tsp1,ix_pts), *_FIELD_AT(phi_ts,ix_pts));

			// C term   : useless C^ast = B

			_spinor_prod_assign_g(C1, *_FIELD_AT(phi_ts,ix_pts_pa), *_FIELD_AT(phi_tsp1,ix_pts_pa));
			_spinor_prod_assign_g(C2, *_FIELD_AT(phi_tsp1,ix_pts), *_FIELD_AT(phi_ts,ix_pts));


			// D term
			_spinor_prod_assign_g(D1, *_FIELD_AT(phi_ts,ix_pts_pa), *_FIELD_AT(phi_ts,ix_pts_pa));
			_spinor_prod_assign_g(D2, *_FIELD_AT(phi_tsp1,ix_pts), *_FIELD_AT(phi_tsp1,ix_pts));

		}/* end loops on space */

		_complex_mul_assign(corr[0][tc],A1,A2);
		_complex_mul_assign(corr[1][tc],B1,B2);
		_complex_mul_assign(corr[2][tc],C1,C2);
		_complex_mul_assign(corr[3][tc],D1,D2);

		corr_tmp[0][tc] = A1.re;
		corr_tmp[1][tc] = A2.re;

	} /* end loop on time */



	// split real and imaginary part of the correlators

	for (t=0;t<T;t++){

		tc = (zerocoord[0]+t+GLB_T-ts)%GLB_T;

		for (l=0;l<4;l++) {

			if (l==1 || l ==2) norm=-GLB_VOL3*GLB_VOL3; // B and C term have a global  minus sign  in the contraction
			else	  norm= GLB_VOL3*GLB_VOL3;	 // A and B

			corr_re[l][tc]=corr[l][tc].re/norm;
			corr_im[l][tc]=corr[l][tc].im/norm;
		}


		// printf("%i %i %3.10e %3.10e \n",tc,(tc+GLB_T)%GLB_T,corr_tmp[0][(tc+GLB_T)%GLB_T],corr_tmp[1][tc]);
		corr_2pt[tc] = 0.5 *(corr_tmp[0][(tc+GLB_T)%GLB_T]+corr_tmp[1][tc])/(GLB_VOL3);
		// corr_2pt[tc] =(corr_tmp[1][tc])/(GLB_X*GLB_Y*GLB_Z);

	}


	/* global sums */
	for (l=0;l<4;l++) {
		global_sum(corr_re[l],GLB_T);
		global_sum(corr_im[l],GLB_T);
	}
	global_sum(corr_2pt,GLB_T);

	lprintf("TEST_PP",0," Correlators C_pipi : T l nhits re im ! \n");
	for(l=0;l<4;l++ ) for(t=0;t<GLB_T;++t) lprintf("TEST_PP",0,"%i %i %i %3.10e %3.10e \n",t,l,k,corr_re[l][t],corr_im[l][t]);

	lprintf("TEST_P",0," Correlators C_pi : T  nhits re ! \n");
	for(t=0;t<GLB_T;++t)  lprintf("TEST_P",0,"%i %i %3.10e  \n",t,k,corr_2pt[t]);
	fflush(stdout);
	gettimeofday(&end,0);
	timeval_subtract(&etime,&end,&start);
	lprintf("TIMING",0,"Contractions for Pion-Pion ( nhits = %i ) scattering done in  [%ld sec %ld usec]\n",k,etime.tv_sec,etime.tv_usec);

	}

