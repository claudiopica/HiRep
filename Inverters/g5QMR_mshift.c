#include "inverters.h"
#include "linear_algebra.h"
#include "complex.h"
#include "malloc.h"
#include "update.h"
#include "utils.h"
#include <stdio.h>
#include <math.h>
#include <assert.h>

#define _print_par(c) \
  printf("[%d] " #c " = %e\n",cgiter,c)

/* Calculate the Givens rotation (c,s) that maps:
 * |  c   s | | a | = | 1/r |
 * | -s^+ c | | b |   | 0 |
 */
#define _Givens_rot(r,c,s,a,b) \
  { if (fabs(a)<1.e-10) { (r)=(b); (c)=0.; (s)=1.; }		\
    else { (s)=(b)/(a); (c)=1./sqrt(1.+(s)*(s)); (s)*=(c); (r)=1./((c)*(a)+(s)*(b)); } \
  } 

/*
 * performs the multi-shifted QMR inversion for g5-hermitean matrices:
 * out[i] = (M-(par->shift[i]))^-1 in
 * returns the number of cg iterations done.
 */
static int g5QMR_mshift_core(short *valid, mshift_par *par, spinor_operator M, suNf_spinor *in, suNf_spinor_dble **out){

	suNf_spinor_dble **q1,**q2;
	suNf_spinor_dble *p1, *p2, *Mp, *sd;
	suNf_spinor_dble *sptmp, *memall;

	double alpha, beta, delta, rho, innorm2; 
	double *r, *s1, *s2, *c1, *c2;

	int i;
	int cgiter;
	unsigned int notconverged;
	unsigned int spinorlen;

	unsigned short *flags;

	/* fare qualche check sugli input */
	/* par->n deve essere almeno 1! */
	assert(par->n>0);
	/*
		 printf("numero vettori n=%d\n",par->n);
		 for (i=0; i<(par->n); ++i) {
		 printf("shift[%d]=%f\n",i,par->shift[i]);
		 printf("out[%d]=%p\n",i,out[i]);      
		 }
		 */

	/* allocate spinors fields and aux real variables */
	/* implementation note: to minimize the number of malloc calls
	 * objects of the same type are allocated together
	 */
	get_spinor_len(&spinorlen);
	q1 = (suNf_spinor_dble **)malloc(sizeof(suNf_spinor_dble*)*2*(par->n));
	q2 = q1+(par->n);
	q1[0] = (suNf_spinor_dble *)malloc(sizeof(suNf_spinor_dble)*(2*(par->n)+4)*(spinorlen));
	q2[0] = q1[0]+(par->n)*(spinorlen);
	memall=q1[0];
	for (i=1; i<(par->n); ++i) {
		q1[i] = q1[i-1]+(spinorlen);
		q2[i] = q2[i-1]+(spinorlen);
	}
	p1 = q2[par->n-1]+(spinorlen);
	p2 = p1+(spinorlen);
	Mp = p2+(spinorlen);
	sd = Mp+(spinorlen);

	r = (double *)malloc(sizeof(double)*5*(par->n));
	s1 = r+(par->n);
	s2 = s1+(par->n);
	c1 = s2+(par->n);
	c2 = c1+(par->n);

	flags=(unsigned short *)malloc(sizeof(unsigned short)*(par->n));

	/* init recursion */
	cgiter = 0;
	notconverged=par->n;

	/* spinor_field_copy_dble_f(p2, in); */
	assign_s2sd(spinorlen,p2, in); /* trial solution = 0 */
	innorm2=spinor_field_sqnorm_dble_f(p2);
	if(par->n==1) { /* in this case is not a multishift and we use as starting vector out[0] */
		assign_sd2s(spinorlen,(suNf_spinor*)sd,out[0]);
		M((suNf_spinor*)Mp,(suNf_spinor*)sd);
		assign_s2sd(spinorlen,Mp,(suNf_spinor*)Mp);
		spinor_field_mul_add_assign_dble_f(Mp,-par->shift[0],out[0]);
		spinor_field_sub_dble_f(p2,p2,Mp);
	}
	rho=sqrt(spinor_field_sqnorm_dble_f(p2));
	spinor_field_mul_dble_f(p2,1./rho,p2);
	spinor_field_zero_dble_f(p1);  
	for (i=0; i<(par->n); ++i) {
		r[i]=rho;
		c2[i]=c1[i]=1.;
		s1[i]=s2[i]=0.;
		if (par->n!=1) spinor_field_zero_dble_f(out[i]); /* if no multishift we start with the trial solution */
		spinor_field_zero_dble_f(q1[i]);
		spinor_field_zero_dble_f(q2[i]);
		flags[i]=1;
	}
	delta = spinor_field_g5_prod_re_dble_f(p2,p2);
	beta = delta;

	/* cg recursion */
	do {
		++cgiter;

		assign_sd2s(spinorlen,(suNf_spinor*)sd,p2);
		M((suNf_spinor*)Mp,(suNf_spinor*)sd);
		assign_s2sd(spinorlen,Mp,(suNf_spinor*)Mp);
		spinor_field_mul_add_assign_dble_f(Mp,-par->shift[0],p2);

		/* compute alpha */
		alpha = spinor_field_g5_prod_re_dble_f(p2,Mp)/delta;

		/* update p1, p2 */
		spinor_field_mul_add_assign_dble_f(Mp,-beta,p1);
		spinor_field_mul_add_assign_dble_f(Mp,-alpha,p2);
		sptmp=p1;
		p1=p2;
		p2=Mp;
		Mp=sptmp;

		/* update rho */
		rho=sqrt(spinor_field_sqnorm_dble_f(p2));

		for (i=0; i<(par->n); ++i) { /* update solutions */
			if (flags[i]) {
				double a, t, e, d, m;
				a=(alpha-par->shift[i]+par->shift[0]);
				t=s1[i]*beta;
				e=c2[i]*c1[i]*beta+s2[i]*a;
				m=c2[i]*a-s2[i]*c1[i]*beta;
				s1[i]=s2[i];
				c1[i]=c2[i];
				_Givens_rot(d,c2[i],s2[i],m,rho);

				/* update q */
				spinor_field_lc_dble_f(Mp,-t*d,q1[i],-e*d,q2[i]);
				spinor_field_mul_add_assign_dble_f(Mp,d,p1);
				sptmp=q1[i];
				q1[i]=q2[i];
				q2[i]=Mp;
				Mp=sptmp;/* swap q1[i]<-q2[i]<-Mp and Mp point to q1 */

				/* update solution */
				spinor_field_mul_add_assign_dble_f(out[i],c2[i]*r[i],q2[i]);

				/* update residuum */
				r[i]*=-s2[i];	

				if(fabs(r[i]*r[i])*(double)(1+cgiter)<par->err2*innorm2){
					flags[i]=0;
					--notconverged;
#ifndef NDEBUG
					printf("[%d] converged at iter: %d\n",i,cgiter);
#endif
				}

			}    
		}

		if (fabs(rho)> 1.e-10) {
			double olddelta;

			/*normalize p2 */
			spinor_field_mul_dble_f(p2,1./rho,p2);

			/* update delta and beta */
			olddelta=delta;
			delta = spinor_field_g5_prod_re_dble_f(p2,p2);

			if(fabs(delta)<1.e-6) { /* the method has failed ! */
				printf("g5QMR: Method failed! (delta=%e) (cgiter=%d)\n",delta,cgiter);
				notconverged=0; /* exit loop */
			}

			/* update beta */
			beta=rho*delta/olddelta;

		} else { /* system has converged */
			notconverged=0;
			printf("g5QMR: rho = 0 ! (system converged)\n");
		}


#ifndef NDEBUG
		if(cgiter%100==0){
			printf("g5QMR: [%d] res[0]=%e notconverged=%d\n",cgiter,fabs(r[0]*r[0])*(double)(1+cgiter)/innorm2,notconverged);   
		}
#endif

	} while ((par->max_iter==0 || cgiter<par->max_iter) && notconverged);

	/* test results */
	for(i=0;i<par->n;++i){
		float norm;
		assign_sd2s(spinorlen,(suNf_spinor*)sd,out[i]);
		M((suNf_spinor*)Mp,(suNf_spinor*)sd);
		assign_s2sd(spinorlen,Mp,(suNf_spinor*)Mp);
		if(par->shift[i]!=0.) {
			spinor_field_mul_add_assign_dble_f(Mp,-par->shift[i],out[i]);
		}
		assign_s2sd(spinorlen,sd,in);
		spinor_field_mul_add_assign_dble_f(Mp,-1.0,sd);
		norm=spinor_field_sqnorm_dble_f(Mp)/((double)spinor_field_sqnorm_f(in));
		if (fabs(norm)>par->err2){
			valid[i]=0;
#ifndef NDEBUG
			printf("g5QMR Failed: err2[%d] = %e  shift = %e\n",i,norm,par->shift[i]);
#endif
		} else {
			valid[i]=1;
		}
	}
#ifndef NDEBUG
	printf("g5QMR: %d matrix multiplications\n",cgiter);
#endif

	/* free memory */
	free(memall);
	free(q1);
	free(r);

	free(flags);

	/* return number of cg iter */
	if(par->n==1) cgiter++;
	return cgiter;
}

static float sh;
static spinor_operator g5Herm;
static void Herm(suNf_spinor *out, suNf_spinor *in){
  g5Herm(out,in);
	if(sh!=0.) {
		spinor_field_mul_add_assign_f(out,-sh,in);
	}
	spinor_field_g5_f(out,out);
}

int g5QMR_mshift(mshift_par *par, spinor_operator M, suNf_spinor *in, suNf_spinor_dble **out){
	int cgiter=0;
	int n;
	mshift_par orig;
	short *valid;
	int loccg;
	int debugint;

	orig=*par; /* save par */
	valid=malloc(sizeof(short)*orig.n);
	cgiter+=g5QMR_mshift_core(valid,par,M,in,out);
	debugint=cgiter;

	/* if some vector has not converged try non-multishift */
	par->n=1;
	for(n=0;n<orig.n;++n){
		if (valid[n]==0) {
#ifndef NDEBUG
			printf("Doing non multishift on vector %d\n",n);
#endif
			par->shift=orig.shift+n;
			/* par->max_iter=100; */
			loccg=g5QMR_mshift_core(valid+n,par,M,in,out+n);
			cgiter+=loccg;

			if (valid[n]==0) {
				/* revert to MINRES */
				unsigned int spinorlen;
				MINRES_par Mpar;
				
				get_spinor_len(&spinorlen);
				Mpar.err2=par->err2;
				Mpar.max_iter=0;

				g5Herm=M;
				sh=par->shift[0];

				spinor_field_g5_f(in,in);
				assign_sd2s(spinorlen,(suNf_spinor*)out[n],out[n]);
				loccg=MINRES(&Mpar,&Herm,in,(suNf_spinor*)out[n],(suNf_spinor*)out[n]);
				spinor_field_g5_f(in,in);
				assign_s2sd(spinorlen,out[n],(suNf_spinor*)out[n]);
				/* give up */
				printf("Used MINRES on vector %d [%d]\n",n,loccg);
				cgiter+=loccg;
			}
		}
	}

	/* this is for debug purposes */
	printf("g5QMR matrix mul: %d/%d\n",debugint,cgiter);
	
	*par=orig; /* restore par */

	free(valid);
	return cgiter;

}
