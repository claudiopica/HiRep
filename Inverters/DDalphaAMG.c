/***************************************************************************\
* Copyright (c) 2013, Ari Hietanen & Ulrik Ishoej Soendergaard              *   
* All rights reserved.                                                      * 
\***************************************************************************/
#include <stdlib.h>

#include "geometry.h" 
#include "global.h" 
#include "error.h"
#include "logger.h"
#include "inverters.h" 
#include "linear_algebra.h"
#include "memory.h"
#include "update.h"

static spinor_field *P=NULL;


void DDalphaAMG_finalize(){
	free_spinor_field_f(P);
	P=NULL;
}

void DDalphaAMG_setup(mshift_par *par, spinor_operator M, int N, int nu, int n_inv){ 

    lprintf("DDalphaAMG",0,"Setup with N=%d, nu=%d and n_inv=%d\n",N,nu,n_inv);
    
	spinor_field *v,*tmp;				    // v is to contain the locally orthonormalized near kernel vectors
												// P is to contain the locally orthonormalized near kernel vectors
if (P==NULL){
	P=alloc_spinor_field_f(N,&glattice);		// &glattice might be changed to &even	
    atexit(DDalphaAMG_finalize);
}
	v=alloc_spinor_field_f(N+1,&glattice);
	tmp=v+N;
	
	complex *r_upper=malloc(1*sizeof(*r_upper));
	complex *r_lower=malloc(1*sizeof(*r_upper));
	
	lprintf("DDalphaAMG",0,"Removing all but near kernel with SAP.\n");
	for (int i=0;i<N;i++){
		gaussian_spinor_field(&v[i]);		
			for (int eta=1;eta<4;eta++){		// Could be optimized
				SAP_prec(eta,&cg_mshift,par, M, &v[i], tmp);
				spinor_field_copy_f(&v[i],tmp);
			}
	}	// Line 6
	
	for (int eta=1;eta<=n_inv;eta++){
	
	lprintf("DDalphaAMG",10,"Building interpolating operators (%d/%d).\n",eta,n_inv);
	// ------------------------------- Line 8: Reconstruct P from v-vectors
	// --------------- Gram-Schmidt orthogonalization
	//Define r11 := ∥x1∥2. If r11 = 0 Stop, else q1 := x1/r11.
	
	spinor_field_sqnorm_aggregate_f(&(r_upper[0].re) , &(r_lower[0].re), &v[0]);
	r_upper[0].im=0;r_lower[0].im=0;
	r_upper[0].re=sqrt(r_upper[0].re);
	r_lower[0].re=sqrt(r_lower[0].re);
	
	spinor_field_mul_aggregate_f(&P[0],1./r_upper[0].re,1./r_lower[0].re,&v[0]);
	
	// DEBUG START
			spinor_field_sqnorm_aggregate_f(&(r_upper[0].im) , &(r_lower[0].im), &P[0]);
			lprintf("DDalphaAMG",10,"Sqnorm of P[0] upper = %e, lower = %e).\n",r_upper[0].im,r_lower[0].im);
			r_upper[0].im=spinor_field_sqnorm_f(&P[0]);
			lprintf("DDalphaAMG",10,"global sqnorm of P[0] = %e).\n",r_upper[0].im);
	// DEBUG END

	for (int j=1;j<N;j++){
		spinor_field_copy_f(&P[j],&v[j]);
		for (int i=0;i<N;i++){ // Orthogonaliza against alle previous vectors
			spinor_field_prod_aggregate_f(r_upper,r_lower,&P[i],&P[j]);
			_complex_minus(r_upper[0],r_upper[0]);
			_complex_minus(r_lower[0],r_lower[0]);
			spinor_field_mulc_add_assign_aggregate_f(&P[j],r_upper[0],r_lower[0],&P[i]); 
		}
		spinor_field_sqnorm_aggregate_f(&(r_upper[0].re) , &(r_lower[0].re), &P[j]);
		r_upper[0].im=0;r_lower[0].im=0;
		r_upper[0].re=sqrt(r_upper[0].re);
		r_lower[0].re=sqrt(r_lower[0].re);
		spinor_field_mul_aggregate_f(&P[j],1./r_upper[0].re,1./r_lower[0].re,&P[j]); // normalize
	} // Later we might want to add reorthogonalization
	
	// DEBUG START
			spinor_field_sqnorm_aggregate_f(&(r_upper[0].im) , &(r_lower[0].im), &P[4]);
			lprintf("DDalphaAMG",10,"Sqnorm of P[4] upper = %e, lower = %e).\n",r_upper[0].im,r_lower[0].im);
			r_upper[0].im=spinor_field_sqnorm_f(&P[4]);
			lprintf("DDalphaAMG",10,"global sqnorm of P[4] = %e).\n",r_upper[0].im);
			spinor_field_prod_aggregate_f(r_upper,r_lower,&P[4],&P[2]);
			lprintf("DDalphaAMG",10,"Aggregate (P[4]^*,P[2]) upper = (%e , %e ).\n",r_upper[0].re,r_upper[0].im);
			lprintf("DDalphaAMG",10,"Aggregate (P[4]^*,P[2]) lower = (%e , %e ).\n",r_lower[0].re,r_lower[0].im);
	// DEBUG END
	
	} // eta
	
	free_spinor_field_f(v);
	v=NULL;
	free(r_upper);
	free(r_lower);
	lprintf("DDalphaAMG",0,"Setup done.");
}


// 1 V-cycle (C operation)
void DDalphaAMG(mshift_par *par, spinor_operator M, spinor_field *in, spinor_field *out){ 
	SAP_prec(3,&cg_mshift,par, M, in, out);
}
