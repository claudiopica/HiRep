/***************************************************************************\
* Copyright (c) 2013, Ari Hietanen & Ulrik Ishoej Soendergaard              *   
* All rights reserved.                                                      * 
\***************************************************************************/
#include <stdlib.h>

#include "geometry.h" 
#include "global.h" 
#include "error.h"
#include "logger.h"
#include "io.h"
#include "moreio.h"
#include "inverters.h" 
#include "linear_algebra.h"
#include "memory.h"
#include "update.h"
#include "complex.h"

static spinor_field *P=NULL;
static int NumVectors=0;
static complex *Dc=NULL;

void DDalphaAMG_finalize(){
	free_spinor_field_f(P);
	free(Dc);
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
}*/


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
    
	spinor_field *v,*tmp,*tmp2;				    // v is to contain the locally orthonormalized near kernel vectors
												// P is to contain the locally orthonormalized near kernel vectors
if (P==NULL){
	P=alloc_spinor_field_f(N,&glattice);		// &glattice might be changed to &even	
	Dc = malloc(4*N*N*sizeof(*Dc));
	NumVectors=N;
    atexit(DDalphaAMG_finalize);
}
	v=alloc_spinor_field_f(N,&glattice);
	tmp=alloc_spinor_field_f(1,&glattice);
	tmp2=alloc_spinor_field_f(2,&glattice);
	
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
		for (int i=0;i<N;i++){ // Orthogonaliza against all previous vectors
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
	
	
	if(eta<n_inv){ // This if is not in the paper, but I think it should be
		lprintf("DDalphaAMG",10,"Improving near kernel\n");
			for(int j=0;j<N;j++){ //  ------------- Line 9
			
				M.dbl(tmp,&v[j]);
				spinor_field_sub_f(tmp2,&v[j],tmp); // Tmp now contains the residua
				DDalphaAMG(par,M,tmp,tmp2);
				spinor_field_add_assign_f(&v[j],tmp); // --Line 10
			
			}  // j //-------------   Line 12
		} // if eta 
	} // eta // ------------- Line 13
	
	free_spinor_field_f(v);v=NULL;
	free_spinor_field_f(tmp);tmp=NULL;
	free_spinor_field_f(tmp2);tmp2=NULL;
	free(r_upper);r_upper=NULL;
	free(r_lower);r_lower=NULL;
	lprintf("DDalphaAMG",0,"Setup done.\n");
}


// 1 V-cycle (C operation)
void DDalphaAMG(mshift_par *par, spinor_operator M, spinor_field *in, spinor_field *out){ 
	int nu=3;
//	SAP_prec(3,&cg_mshift,par, M, in, out);

	spinor_field *tmp,*tmp2;
	tmp=alloc_spinor_field_f(1,&glattice);
	tmp2=alloc_spinor_field_f(1,&glattice);
	complex *rhs;
	complex *sol;
	rhs=malloc(2*NumVectors*sizeof(*rhs));
	sol=malloc(2*NumVectors*sizeof(*sol));
	
	spinor_field_to_course_spinor_field_new(rhs,in);
	
	
	coarse_MINRES(sol,rhs);
	
 	coarse_spinor_field_to_spinor_field_new(out,sol);  // Now out has the course inversion result
 	
 	// At this point we should maybe syncronize
 	// Now we use the smoother on the residual
 	
 	
 	M.dbl(tmp,out);
 	spinor_field_sub_assign_f(tmp,in);
 	// Now tmp contains the residual
 	
    spinor_field_zero_f(tmp2);
 	SAP_prec(nu,&cg_mshift,par, M, tmp, tmp2);
 	spinor_field_sub_assign_f(out,tmp2);
 	
 	free(rhs);rhs=NULL;
 	free(sol);sol=NULL;
	free_spinor_field_f(tmp);
	free_spinor_field_f(tmp2);
}

// (small) = P^H * (big)
void spinor_field_to_course_spinor_field_new(complex *smalls,spinor_field *big){
	lprintf("DDalphaAMG",10,"spinor_field_to_course_spinor_field_new called.\n");
	for(int i=0;i<NumVectors;i++){
		spinor_field_prod_aggregate_f(&smalls[i],&smalls[i+NumVectors],&P[i],big);
	}
	
	lprintf("DDalphaAMG",10,"spinor_field_to_course_spinor_field_new complete.\n");
}
// (big) = P * (small)
void coarse_spinor_field_to_spinor_field_new(spinor_field *b,complex *smalls){
	lprintf("DDalphaAMG",10,"coarse_spinor_field_to_spinor_field_new called.\n");
	
	spinor_field_mulc_aggregate_f(b,smalls[0],smalls[NumVectors],&P[0]);

	for(int i=1;i<NumVectors;i++){
		spinor_field_mulc_add_assign_aggregate_f(b,smalls[i],smalls[i+NumVectors],&P[i]);
	}


}


// (small)=Dc (small)
void coarse_operation(complex *a, complex *b){
	for(int i=0;i<2*NumVectors;i++){
		_complex_0(a[i]);
		
			for(int j=0;j<2*NumVectors;j++){
				_complex_mul_assign(a[i],Dc[i*2*NumVectors+j],b[j]);
			}
	}
}



void coarse_spinor_field_copy(complex *to, complex *from){
	lprintf("DDalphaAMG",10,"coarse_spinor_field_copy called.\n");
	for(int i=0;i<2*NumVectors;i++){
		to[i]=from[i];
	}
	lprintf("DDalphaAMG",10,"coarse_spinor_field_copy complete.\n");
}
void coarse_mulc_add_assign(complex *to,complex c, complex *from){
	for(int i=0;i<2*NumVectors;i++){
		_complex_mul_assign(to[i],c,from[i]);
	}
}
void coarse_mulc_sub_assign(complex *to, complex c, complex *from){
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
	
	for (int i=0;i<2*NumVectors;i++){
		_complex_prod_assign(res,a[i],b[i]);
	}
	
	return(res);
}
double coarse_sqnorm(complex *a){
	double res=0;
	
	for (int i=0;i<2*NumVectors;i++){
		res+=_complex_re(a[i])*_complex_re(a[i])+_complex_im(a[i])*_complex_im(a[i]);
	}
	
	return(res);
}

void coarse_MINRES(complex *x,complex *b){

	lprintf("DDalphaAMG",10,"coarse_MINRES called.\n");
	int iter=0;
	complex alpha;
	double tmp,residual,bnorm;
	complex *r;
	r=malloc(2*NumVectors*sizeof(*r));
	complex *p;
	p=malloc(2*NumVectors*sizeof(*p));

	coarse_spinor_field_copy(r,b); // line 1 (Saad p 165)	
	coarse_zero(x);	
	coarse_operation(p,r); 

	residual=coarse_sqnorm(r);
	bnorm=coarse_sqnorm(b);
	
	
//	lprintf("coarse_MINRES",10,"Residual(0) = %e.\n", residual);
	
	while(residual/bnorm>1e-5  &&  iter<100){	// line 2
		alpha=coarse_inner_product(p,r); // I think this should be the right one

		tmp=coarse_sqnorm(p);
		_complex_mulr(alpha,1/tmp,alpha); // line 3
		coarse_mulc_add_assign(x,alpha,r); // line 4
		
		coarse_mulc_sub_assign(r,alpha,p); // line 5		
		coarse_operation(p,r); // line 6

	residual=coarse_sqnorm(r);
	iter++;
	if(iter%10==0){lprintf("coarse_MINRES",10,"r(%d) = %e.\n", iter,residual/bnorm);}
	}	// line 7
	
	
	lprintf("DDalphaAMG",10,"coarse_MINRES Terminated r(%d) = %e .\n",iter,residual/bnorm);

	free(r);r=NULL;
	free(p);p=NULL;
}


