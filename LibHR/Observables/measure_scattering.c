#include "spin_matrix.h"
#include "complex.h"
#include "global.h"
#include "meson_observables.h"
#include <math.h>
#include "communications.h"
#include "logger.h"
#include <stdlib.h>

#define PI 3.141592653589793238462643383279502884197

/* TJ - Extracting Vincent's meson-meson scattering measurements so that we can put them in LibHR

   Also adding momenta to the sinks
   Scattering is organised in four contractions:

   A.	src1<=====>snk1
   src2<=====>snk2

   B.   src1<----->snk1
   \   /
   X
   / \
   src2<---->snk2

   C.   src1<----->snk2
   \   /
   X
   / \
   src2<---->snk1

   D.	src1<=====>snk2
   src2<=====>snk1
   */

//Some macros for indexing
#define INDEX(px,py,pz,n_mom,tc) ((px + n_mom)*(2*n_mom+1)*(2*n_mom+1)*(GLB_T)+(py + n_mom)*(2*n_mom+1)*(GLB_T)+(pz + n_mom)*(GLB_T)+ (tc))

//Useful for splitting mesons at the sink
int tsplit(int ipt, int delta)
{
	if(delta>0)
	{
		return iup(ipt,0);
	}
	else if(delta<0)
	{
		return idn(ipt,0);
	}
	else
	{
		return ipt;
	}
}

// psi0 to psi3 - propagators for the mesons - psi0 and psi1 correspond to meson 1 and psi2 and psi3 to meson 2. Propagators psi1 and psi3 are complex conjugated in the contractions, so their momenta will be opposite to the original. Propagators psi0 and psi1 (as well as psi3 and psi4) should start on the same time slice and be generated with the same noise vector. To avoid issues with Fierz rearrangement, psi0/1 and psi2/3 should either be inverted at different time slices or with different noise vectors.
// tau - time slice of the source plane
// split - splitting of mesons at the sink, the time in the correlator is measured from tau to the EARLIER pion (pion1 if split is positive, pion2 if split is negative), this number cannot exceed the valur of T_BORDER
// n_mom - maximum magnitude of momentum for pion 1 in any direction (e.g. n_mom=1 gives 27 possibilities between (-1,-1,-1) and (1,1,1))
// p_tot_r the total momentum in the r direction. The phase factors at the sink are exp(-ip.x), so this should be equal to the total momentum at the source.
// mo - where the output is saved

void measure_scattering_AD_core(meson_observable* mo, spinor_field* psi0,spinor_field* psi1,spinor_field* psi2,spinor_field* psi3, int tau, int split, int n_mom, int p_tot_x, int p_tot_y, int p_tot_z){

	int px, py, pz, t, x, y, z, ix, ix_split, beta, tc, splittmp;
	complex  trtmp;
	double tr1re[GLB_T],tr2re[GLB_T],tr1im[GLB_T],tr2im[GLB_T];
	double pdotx, cpdotx, spdotx;
	suNf_spin_matrix sm0, sm1, sm2, sm3;
	meson_observable* motmp;

	if((split < -1) || split > 1)
	{
		lprintf("ERROR",0,"The value of split can not exceed 1");
		exit(-1);
	}
	// splittmp will be used to ensure that time is always measured to the earlier pion
	splittmp = split>0 ? 0 : split;
	for (px=-n_mom;px<=n_mom;++px) for (py=-n_mom;py<=n_mom;++py) for (pz=-n_mom;pz<=n_mom;++pz){
		for(tc=0;tc<GLB_T;++tc){
			tr1re[tc]=0.0;
			tr1im[tc]=0.0;
			tr2re[tc]=0.0;
			tr2im[tc]=0.0;
		}
		for (t=0; t<T; t++) {	 
			// Correlator time measured from tau to the earlier pion
			tc = (zerocoord[0]+t+splittmp+GLB_T-tau)%GLB_T;
			for (x=0; x<X; x++) for (y=0; y<Y; y++) for (z=0; z<Z; z++) { 
				ix=ipt(t,x,y,z);
				//ix_split=ipt(t+split,x,y,z);
				ix_split=tsplit(ix, split);

				for (beta=0;beta<4;beta++)
				{
					_spinmatrix_assign_row(sm0, *_FIELD_AT(&psi0[beta],ix), beta);
					_spinmatrix_assign_row(sm1, *_FIELD_AT(&psi1[beta],ix), beta);
					_spinmatrix_assign_row(sm2, *_FIELD_AT(&psi2[beta],ix_split), beta);
					_spinmatrix_assign_row(sm3, *_FIELD_AT(&psi3[beta],ix_split), beta);
				}
				//Contracting pion 1
				_spinmatrix_mul_trace(trtmp, sm1, sm0); // Tr(SM1^\dagger SM0)
				pdotx = 2.0*PI*( ((double) px)*(x+zerocoord[1])/GLB_X + ((double) py)*(y+zerocoord[2])/GLB_Y + ((double) pz)*(z+zerocoord[3])/GLB_Z);
				cpdotx=cos(pdotx);
				spdotx=sin(pdotx);
				tr1re[tc]+=trtmp.re*cpdotx+trtmp.im*spdotx;
				tr1im[tc]+=trtmp.im*cpdotx-trtmp.re*spdotx;

				//Contracting pion 2
				_spinmatrix_mul_trace(trtmp, sm3, sm2); // Tr(SM3^\dagger SM2)
				pdotx = 2.0*PI*( ((double) (p_tot_x - px) )*(x+zerocoord[1])/GLB_X + ((double) (p_tot_y - py))*(y+zerocoord[2])/GLB_Y + ((double) (p_tot_z - pz))*(z+zerocoord[3])/GLB_Z);
				cpdotx=cos(pdotx);
				spdotx=sin(pdotx);
				tr2re[tc]+=trtmp.re*cpdotx+trtmp.im*spdotx;
				tr2im[tc]+=trtmp.im*cpdotx-trtmp.re*spdotx;
			}
		}

		global_sum(tr1re, GLB_T);
		global_sum(tr1im, GLB_T);
		global_sum(tr2re, GLB_T);
		global_sum(tr2im, GLB_T);

		for (t=0; t<GLB_T; t++) {	 
			tc = (zerocoord[0]+t+splittmp+GLB_T-tau)%GLB_T;
			motmp=mo;
			while (motmp!=NULL){
				// I chose not to divide by the volume to conform with Ari and Rudy's conventions
				motmp->corr_re[INDEX(px,py,pz,n_mom,tc)] = (tr1re[tc]*tr2re[tc] - tr1im[tc]*tr2im[tc]);//(GLB_VOL3*GLB_VOL3);
				motmp->corr_im[INDEX(px,py,pz,n_mom,tc)] = (tr1re[tc]*tr2im[tc] + tr2re[tc]*tr1im[tc]);//(GLB_VOL3*GLB_VOL3);
				motmp=motmp->next;
			}
		}
	}
	/*  motmp=mo;
	    while (motmp!=NULL){
	    global_sum(motmp->corr_re, (2*n_mom+1)*(2*n_mom+1)*(2*n_mom+1)*GLB_T);
	    global_sum(motmp->corr_im, (2*n_mom+1)*(2*n_mom+1)*(2*n_mom+1)*GLB_T);
	    motmp=motmp->next;
	    }*/
}

//Same arguments as the AD function
void measure_scattering_BC_core(meson_observable* mo, spinor_field* psi0,spinor_field* psi1, spinor_field* psi2,spinor_field* psi3, int tau, int split, int n_mom, int p_tot_x, int p_tot_y, int p_tot_z){

	int px, py, pz, t, x, y, z, ix, ix_split, mu, gamma, tc, splittmp;
	complex trace, phase;
	double pdotx;
	suNf_spin_matrix sm00, sm11, sm21, sm30;
	complex B1[4][4], B2[4][4];
	double B1re[16*GLB_T], B1im[16*GLB_T], B2re[16*GLB_T], B2im[16*GLB_T];
	meson_observable* motmp;

	if((split < -1) || split > 1)
	{
		lprintf("ERROR",0,"The value of split=%i can not exceed 1", split);
		exit(-1);
	}
	if((n_mom!=0) || (p_tot_x != 0) || (p_tot_y != 0) || (p_tot_z != 0))
	{
		lprintf("WARNING",0,"The function measure_scattering_BC_core not tested for non-zero momentum!");
	}
	// splittmp will be used to ensure that time is always measured to the earlier pion
	splittmp = split>0 ? 0 : split;
	for (px=-n_mom;px<=n_mom;++px) for (py=-n_mom;py<=n_mom;++py) for (pz=-n_mom;pz<=n_mom;++pz){
		for(t=0;t<16*GLB_T;++t)
		{
			B1re[t]=B1im[t]=B2re[t]=B2im[t]=0.0;
		}
		for (t=0; t<T; t++) {	 
			// Correlator time measured from tau to the earlier pion
			tc = (zerocoord[0]+t+splittmp+GLB_T-tau)%GLB_T;

			//Resetting temporary variables to 0
			_complex_0(trace);
			for (mu=0;mu<4;mu++){ for (gamma=0;gamma<4;gamma++){
				_complex_0(B1[mu][gamma]);
				_complex_0(B2[mu][gamma]);
			}}

			for (x=0; x<X; x++) for (y=0; y<Y; y++) for (z=0; z<Z; z++) { 
				ix=ipt(t,x,y,z);
				//ix_split=ipt(t+split,x,y,z);
				ix_split=tsplit(ix, split);

				for (mu=0;mu<4;mu++)
				{
					// Creating spinor matrices corresponding to various propagators. 
					_spinmatrix_assign_row(sm00, *_FIELD_AT(&psi0[mu],ix), mu);
					_spinmatrix_assign_row(sm11, *_FIELD_AT(&psi1[mu],ix_split), mu);
					_spinmatrix_assign_row(sm21, *_FIELD_AT(&psi2[mu],ix_split), mu);
					_spinmatrix_assign_row(sm30, *_FIELD_AT(&psi3[mu],ix), mu);
				}

				// Multiply by phases now, rest of the code will remain unchanged
				pdotx = 2.0*PI*(((double) px)*(x+zerocoord[1])/GLB_X + ((double) py)*(y+zerocoord[2])/GLB_Y + ((double) pz)*(z+zerocoord[3])/GLB_Z);
				phase.re = cos(pdotx);
				phase.im = -sin(pdotx);
				for (mu=0;mu<4;mu++){ for (gamma=0;gamma<4;gamma++){
					_vector_mulc_f(sm00.c[mu].c[gamma], phase, sm00.c[mu].c[gamma]);
				}}

				pdotx = 2.0*PI*(((double) p_tot_x - (double) px)*(x+zerocoord[1])/GLB_X + ((double) p_tot_y - (double) py)*(y+zerocoord[2])/GLB_Y + ((double) p_tot_z - (double) pz)*(z+zerocoord[3])/GLB_Z);
				phase.re = cos(pdotx);
				phase.im = -sin(pdotx);
				for (mu=0;mu<4;mu++){ for (gamma=0;gamma<4;gamma++){
					_vector_mulc_f(sm21.c[mu].c[gamma], phase, sm21.c[mu].c[gamma]);
				}}
				// compute a new spinmatrix  : (A_mu beta * B_gamma_beta )
				// trace color index -> it becomes a non standard struct. with only spin indices
				// accumulate
				// multiply with second part and finally trace over spin indices. 
				for (mu=0;mu<4;mu++){ for (gamma=0;gamma<4;gamma++) {
					_spinor_prod_assign_f(B1[mu][gamma],sm30.c[gamma],sm00.c[mu]);
					_spinor_prod_assign_f(B2[mu][gamma],sm11.c[gamma],sm21.c[mu]);
				}}
			}
			// Communications
			for (mu=0;mu<4;mu++){ for (gamma=0;gamma<4;gamma++) {
				B1re[4*mu+gamma+16*tc] = B1[mu][gamma].re;
				B1im[4*mu+gamma+16*tc] = B1[mu][gamma].im;
				B2re[4*mu+gamma+16*tc] = B2[mu][gamma].re;
				B2im[4*mu+gamma+16*tc] = B2[mu][gamma].im;
			}}
		}

		global_sum(B1re, 16*GLB_T);
		global_sum(B1im, 16*GLB_T);
		global_sum(B2re, 16*GLB_T);
		global_sum(B2im, 16*GLB_T);
		for (t=0; t<GLB_T; t++) {	 
			tc = (zerocoord[0]+t+splittmp+GLB_T-tau)%GLB_T;
			for (mu=0;mu<4;mu++){ for (gamma=0;gamma<4;gamma++) {
				B1[mu][gamma].re = B1re[4*mu+gamma+16*tc];
				B1[mu][gamma].im = B1im[4*mu+gamma+16*tc];
				B2[mu][gamma].re = B2re[4*mu+gamma+16*tc];
				B2[mu][gamma].im = B2im[4*mu+gamma+16*tc];
			}}
			trace.re=0.0; trace.im=0.0;
			//Trace
			for (mu=0;mu<4;mu++){ for (gamma=0;gamma<4;gamma++) {
				_complex_mul_assign(trace, B1[mu][gamma], B2[gamma][mu]);
			}}
			//It will be a miracle if this works on the first try...
			motmp=mo;
			while (motmp!=NULL){
				motmp->corr_re[INDEX(px,py,pz,n_mom,tc)] = trace.re; //(GLB_VOL3*GLB_VOL3);
				motmp->corr_im[INDEX(px,py,pz,n_mom,tc)] = trace.im; //(GLB_VOL3*GLB_VOL3);
				motmp=motmp->next;
			}
		}
	}
	/*  motmp=mo;
	    while (motmp!=NULL){
	    global_sum(motmp->corr_re, (2*n_mom+1)*(2*n_mom+1)*(2*n_mom+1)*GLB_T);
	    global_sum(motmp->corr_im, (2*n_mom+1)*(2*n_mom+1)*(2*n_mom+1)*GLB_T);
	    motmp=motmp->next;
	    }*/
}

