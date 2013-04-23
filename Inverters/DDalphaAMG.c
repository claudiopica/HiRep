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
#include "complex.h"

static spinor_field *P=NULL;
static int NumVectors=0;
static complex *small_spinor=NULL;
static complex *small_result_spinor=NULL;
static complex *Dc=NULL;

void DDalphaAMG_finalize(){
	free_spinor_field_f(P);
	free(small_spinor);
	free(Dc);
	free(small_result_spinor);
	
	small_spinor=NULL;
	small_result_spinor=NULL;
	Dc=NULL;
	P=NULL;
	NumVectors=0;
}

/*
void print_matrix(){
lprintf("DDalphaAMG",10,"\n Real part\n");
	for(int i=0;i<2*NumVectors;i++){
		for(int j=0;j<2*NumVectors;j++){
			lprintf("DDalphaAMG",10,"%+3.3f ",Dc[i*2*NumVectors+j].re);
		}
		lprintf("DDalphaAMG",10,"\n");
	}
lprintf("DDalphaAMG",10,"\n");
lprintf("DDalphaAMG",10,"\n Imaginary part\n");
	for(int i=0;i<2*NumVectors;i++){
		for(int j=0;j<2*NumVectors;j++){
			lprintf("DDalphaAMG",10,"%+3.3f ",Dc[i*2*NumVectors+j].im);
		}
		lprintf("DDalphaAMG",10,"\n");
	}
lprintf("DDalphaAMG",10,"\n");
}
*/

//Dc = P^H D P
// Dc(i,j) -> Dc[i*NumVectors+j]
void build_coarse_operator(spinor_operator M){
	lprintf("DDalphaAMG",10,"Building coarse operator.\n");
	// Starting with the D*P part...
	spinor_field *tmp;
	tmp=alloc_spinor_field_f(4*NumVectors,&glat_nocomm);
	spinor_field *DP;
	DP=tmp+2*NumVectors;


	for(int i=0;i<NumVectors;i++){
		// Moving upper part to first N components of tmp
		// Moving lower part to last N components of tmp
		spinor_field_split_aggregate_f(&tmp[i],&tmp[i+NumVectors],&P[i]);
		 // Maybe a call to empty_buffers(tmp) is necessary
		M.dbl(&DP[i],&tmp[i]);
		M.dbl(&DP[i+NumVectors],&tmp[i+NumVectors]);
		
		// Constructing the Dc matrix of size 2Nx2N where N=NumVectors
		for (int j=0;j<NumVectors;j++){
			spinor_field_prod_aggregate_f(&Dc[j*2*NumVectors+i],&Dc[(j+NumVectors)*2*NumVectors+i],&P[j],&DP[i]);
			spinor_field_prod_aggregate_f(&Dc[j*2*NumVectors+(i+NumVectors)],&Dc[(j+NumVectors)*2*NumVectors+(i+NumVectors)],&P[j],&DP[i+NumVectors]);
		}
	}
	
	// DEBUG START
/*
 			complex *r_upper=malloc(1*sizeof(*r_upper));
 			complex *r_lower=malloc(1*sizeof(*r_lower));
 			spinor_field_sqnorm_aggregate_f(&(r_upper[0].im) , &(r_lower[0].im), &tmp[3]);
 			lprintf("DDalphaAMG",10,"Sqnorm of tmp[3] upper = %e, lower = %e).\n",r_upper[0].im,r_lower[0].im);
 			r_upper[0].im=spinor_field_sqnorm_f(&tmp[3]);
 			lprintf("DDalphaAMG",10,"global sqnorm of tmp[3] = %e).\n",r_upper[0].im);
 			spinor_field_sqnorm_aggregate_f(&(r_upper[0].im) , &(r_lower[0].im), &tmp[3+NumVectors]);
 			lprintf("DDalphaAMG",10,"Sqnorm of tmp[3+NumVectors] upper = %e, lower = %e).\n",r_upper[0].im,r_lower[0].im);
 			r_upper[0].im=spinor_field_sqnorm_f(&tmp[3+NumVectors]);
 			lprintf("DDalphaAMG",10,"global sqnorm of tmp[3+NumVectors] = %e).\n",r_upper[0].im);


 			spinor_field_sqnorm_aggregate_f(&(r_upper[0].im) , &(r_lower[0].im), &DP[3+NumVectors]);
 			lprintf("DDalphaAMG",10,"Sqnorm of DP[3+NumVectors] upper = %e, lower = %e).\n",r_upper[0].im,r_lower[0].im);
 			r_upper[0].im=spinor_field_sqnorm_f(&DP[3+NumVectors]);
 			lprintf("DDalphaAMG",10,"global sqnorm of DP[3+NumVectors] = %e).\n",r_upper[0].im);
 			free(r_upper);r_upper=NULL;
 			free(r_lower);r_lower=NULL;
 */
 	// DEBUG END
	
//	print_matrix();

	free_spinor_field_f(tmp);
}


void DDalphaAMG_setup(mshift_par *par, spinor_operator M, int N, int nu, int n_inv){ 

    lprintf("DDalphaAMG",0,"Setup with N=%d, nu=%d and n_inv=%d\n",N,nu,n_inv);
    
	spinor_field *v,*tmp;				    // v is to contain the locally orthonormalized near kernel vectors
												// P is to contain the locally orthonormalized near kernel vectors
if (P==NULL){
	P=alloc_spinor_field_f(N,&glattice);		// &glattice might be changed to &even	
	small_spinor = malloc(2*N*sizeof(*small_spinor));
	small_result_spinor = malloc(2*N*sizeof(*small_result_spinor));
	Dc = malloc(4*N*N*sizeof(*Dc));
	NumVectors=N;
    atexit(DDalphaAMG_finalize);
}
	v=alloc_spinor_field_f(N+1,&glattice);
	tmp=v+N;
	
	complex *r_upper=malloc(1*sizeof(*r_upper));
	complex *r_lower=malloc(1*sizeof(*r_lower));
	
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
	
	spinor_field_sqnorm_aggregate_f(&(r_upper[0].re) , &(r_lower[0].re), &v[0]);
	r_upper[0].im=0;r_lower[0].im=0;
	r_upper[0].re=sqrt(r_upper[0].re);
	r_lower[0].re=sqrt(r_lower[0].re);
	
	spinor_field_mul_aggregate_f(&P[0],1./r_upper[0].re,1./r_lower[0].re,&v[0]);
	
	// DEBUG START
	//		spinor_field_sqnorm_aggregate_f(&(r_upper[0].im) , &(r_lower[0].im), &P[0]);
	//		lprintf("DDalphaAMG",10,"Sqnorm of P[0] upper = %e, lower = %e).\n",r_upper[0].im,r_lower[0].im);
	//		r_upper[0].im=spinor_field_sqnorm_f(&P[0]);
	//		lprintf("DDalphaAMG",10,"global sqnorm of P[0] = %e).\n",r_upper[0].im);
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
	//		spinor_field_sqnorm_aggregate_f(&(r_upper[0].im) , &(r_lower[0].im), &P[4]);
	//		lprintf("DDalphaAMG",10,"Sqnorm of P[4] upper = %e, lower = %e).\n",r_upper[0].im,r_lower[0].im);
	//		r_upper[0].im=spinor_field_sqnorm_f(&P[4]);
	//		lprintf("DDalphaAMG",10,"global sqnorm of P[4] = %e).\n",r_upper[0].im);
	//		spinor_field_prod_aggregate_f(r_upper,r_lower,&P[4],&P[2]);
	//		lprintf("DDalphaAMG",10,"Aggregate (P[4]^*,P[2]) upper = (%e , %e ).\n",r_upper[0].re,r_upper[0].im);
	//		lprintf("DDalphaAMG",10,"Aggregate (P[4]^*,P[2]) lower = (%e , %e ).\n",r_lower[0].re,r_lower[0].im);
	// DEBUG END
	// --------------- Gram-Schmidt orthogonalization complete
	build_coarse_operator(M);
	
		// DEBUG START

	// DEBUG END
	
	
	} // eta
	
	free_spinor_field_f(v);
	v=NULL;
	free(r_upper);r_upper=NULL;
	free(r_lower);r_lower=NULL;
	lprintf("DDalphaAMG",0,"Setup done.\n");
}


// 1 V-cycle (C operation)
void DDalphaAMG(mshift_par *par, spinor_operator M, spinor_field *in, spinor_field *out){ 

	SAP_prec(3,&cg_mshift,par, M, in, out);
}

// (small) = P^H * (big)
void spinor_field_to_course_spinor_field(spinor_field *big){
	lprintf("DDalphaAMG",10,"spinor_field_to_course_spinor_field called.\n");
	for(int i=0;i<NumVectors;i++){
		spinor_field_prod_aggregate_f(&small_spinor[i],&small_spinor[i+NumVectors],&P[i],big);
	}
}
// (big) = P * (small)
void coarse_spinor_field_to_spinor_field(spinor_field *b){
	lprintf("DDalphaAMG",10,"coarse_spinor_field_to_spinor_field called.\n");
	
	spinor_field_mulc_aggregate_f(b,small_result_spinor[0],small_result_spinor[NumVectors],&P[0]);

	for(int i=0;i<NumVectors;i++){
		spinor_field_mulc_add_assign_aggregate_f(b,small_result_spinor[i],small_result_spinor[i+NumVectors],&P[i]);
	}

}

// (big) = P * (small)
/*void coarse_spinor_field_to_spinor_field(spinor_field *big){
	lprintf("DDalphaAMG",10,"coarse_spinor_field_to_spinor_field called.\n");
	spinor_field_mulc_aggregate_f(big,small_spinor[0],small_spinor[NumVectors],&P[0]);
	
	for(int i=1;i<NumVectors;i++){
		spinor_field_mulc_add_assign_aggregate_f(big,small_spinor[i],small_spinor[i+NumVectors],&P[i]);
	}

}*/
// (small)=Dc (small)
void coarse_spinor_operation(complex *a, complex *b){
	lprintf("DDalphaAMG",10,"coarse_spinor_operation called.\n");
	if(a==NULL){
	for(int i=0;i<2*NumVectors;i++){
	_complex_0(small_result_spinor[i]);
			for(int j=0;j<2*NumVectors;j++){
				_complex_mul_assign(small_result_spinor[i],Dc[i*2*NumVectors+j],small_spinor[j]);
			}
	}
	}else{
	
	for(int i=0;i<2*NumVectors;i++){
	_complex_0(a[i]);
			for(int j=0;j<2*NumVectors;j++){
				_complex_mul_assign(a[i],Dc[i*2*NumVectors+j],b[j]);
			}
	}
	}
}
void coarse_spinor_field_copy(complex *to, complex *from){
	for(int i=0;i<2*NumVectors;i++){
		to[i]=from[i];
	}
}
void coarse_mulc_add_assign(complex *to,complex c, complex *from){
	for(int i=0;i<2*NumVectors;i++){
		_complex_mul_assign(to[i],c,from[i]);
	}
}
void coarse_mulc_sub_assign(complex *to,complex c, complex *from){
	complex tmp;
	_complex_minus(tmp,c);
	for(int i=0;i<2*NumVectors;i++){
		_complex_mul_assign(to[i],tmp,from[i]);
	}
}
void coarse_zero(complex *to){

	for(int i=0;i<2*NumVectors;i++){
		_complex_0(to[i]);
	}
}

complex coarse_inner_product(complex *a,complex *b){
	complex res;
	_complex_0(res);
	
	for (int i=0;i<NumVectors;i++){
		_complex_mul_assign(res,a[i],b[i]);
	}
	
	return(res);
}
double coarse_sqnorm(complex *a){
	double res=0;
	
	for (int i=0;i<NumVectors;i++){
		res+=_complex_re(a[i])*_complex_re(a[i])+_complex_im(a[i])*_complex_im(a[i]);
	}
	
	return(res);
}

void coarse_MINRES(complex *x,complex *b){
	int iter=0;
	complex alpha;
	double tmp,residual;
	static complex *r;
	small_spinor=malloc(2*NumVectors*sizeof(r));
	static complex *p;
	small_spinor=malloc(2*NumVectors*sizeof(r));

	coarse_spinor_field_copy(r,b); // line 1 (Saad p 165)
	coarse_zero(x);
	coarse_spinor_operation(p,r);
	residual=coarse_sqnorm(r);
	
	while(residual>1e-16  ||  iter>500){	// line 2
		alpha=coarse_inner_product(p,r);
		tmp=coarse_sqnorm(p);
		_complex_mulr(alpha,1/tmp,alpha); // line 3
		coarse_mulc_add_assign(x,alpha,r); // line 4
		coarse_mulc_sub_assign(r,alpha,p); // line 5
		coarse_spinor_operation(p,r); // line 6
	residual=coarse_sqnorm(r);
	iter++;
	}	// line 7
	

	free(r);r=NULL;
	free(p);p=NULL;
}


