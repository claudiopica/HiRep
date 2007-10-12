#include "global.h"
#include "linear_algebra.h"
#include "inverters.h"
#include "suN.h"
#include "dirac.h"
#include "utils.h"
#include "update.h"
#include "error.h"
#include "logger.h"
#include "memory.h"
#include <malloc.h>
#include <math.h>
#include <assert.h>
#include <string.h>


/* #define TESTINGMODE */


static double hmass;
static suNf_spinor *h2tmp;

static void H(suNf_spinor *out, suNf_spinor *in);


static void H(suNf_spinor *out, suNf_spinor *in){
	g5Dphi(hmass,out,in);
}


void H2(suNf_spinor *out, suNf_spinor *in){
  g5Dphi(hmass, h2tmp, in);
  g5Dphi(hmass, out, h2tmp);
}



/* void dirac_eva_onemass(int nev,int nevt,int kmax, */
/*         int imax,double omega1,double omega2,double mass, */
/*         suNf_spinor *ev[],double *d,int *status) { */
/* 	int i, j, x, maxh2iter, ie; */
/* 	complex *ptmp1, *ptmp2, *alpha; */
/* 	double *h2values; */
/* 	unsigned int volume; */
/* 	suNf_spinor **h2vectors; */
/* 	suNf_spinor *ws[2]; */
/* 	complex *smallH; */
/* 	complex *smallvectors; */
/* #ifdef TESTINGMODE	 */
/*  	int k; */
/*  	complex smalltest, prod; */
/* #endif */

/* 	int init; */
/* Specifies whether all eigenvectors should be initialized */
/* (init=0), or only the last nevt-nev eigenvectors (init=1) */
/* or none of them (init!=0 or 1) */

/* 	double ubnd; */
/* Upper bound on the eigenvalues of the operator */

/*  	double norm; */
/* 	suNf_spinor *test=0; */
/* 	 */
/* 	get_spinor_len(&volume); */
/* 	error(volume!=VOLUME,1,"dirac_eva_onemass","The spinors must be defined in each site"); */
/* 	 */
/* 	test = alloc_spinor_field_f(); */
/* 	h2values = (double*)malloc(sizeof(double)*nevt); */
/* 	h2vectors = (suNf_spinor**)malloc(sizeof(suNf_spinor*)*nevt); */
/* 	for(i = 0; i < nevt; i++) */
/* 		h2vectors[i] = alloc_spinor_field_f(); */
/* 	ws[0] = alloc_spinor_field_f(); */
/* 	ws[1] = alloc_spinor_field_f(); */
/* 	h2tmp = alloc_spinor_field_f(); */
/* 	smallH = (complex*)malloc(sizeof(complex)*nevt*nevt); */
/* 	smallvectors = (complex*)malloc(sizeof(complex)*nevt*nevt); */
/* 	 */
/* 	maxh2iter = max_H2(&ubnd,mass); */
/* 	init = 0; */

/* 	hmass = mass; */
/* 	 */
/* 	ie = eva(VOLUME,nev,nevt,init,kmax,imax,ubnd,omega1,omega2,&H2,ws,h2vectors,h2values,status); */
/* 	error(ie!=0,1,"dirac_eva_onemass","Failed to compute Dirac eigenvectors"); */

/* 	lprintf("DIRAC_EVA_ONEMASS",10,"Computed %d eigenvetors of H2.\n",nevt); */

/* #ifdef TESTINGMODE	 */
/* 	for(i = 0; i < nevt; i++) { */
/* 		H2(test,h2vectors[i]); */
/* 		spinor_field_mul_add_assign_f(test,-h2values[i],h2vectors[i]); */
/* 		norm=spinor_field_sqnorm_f(test); */
/* 		lprintf("DIRAC_EVA_ONEMASS",40,"H2 eigenvectors test [%d,%e] = %e\n",i,h2values[i],norm); */
/* 	} */

/* 	norm = 0.; */
/* 	for(i = 0; i < nevt; i++) */
/* 	for(j = 0; j < nevt; j++) { */
/* 	   prod.re = spinor_field_prod_re_f(h2vectors[i],h2vectors[j]); */
/* 	   prod.im = spinor_field_prod_im_f(h2vectors[i],h2vectors[j]); */
/* 	   if(i == j) norm += (prod.re-1.)*(prod.re-1.); */
/* 	   else norm += prod.re*prod.re; */
/* 	   norm += prod.im*prod.im; */
/* 	} */
/* 	norm = sqrt(norm); */
/* 	lprintf("DIRAC_EVA_ONEMASS",40,"H2 eigenvectors norm test = %e\n",norm); */

/* 	for(i = 0; i < nevt; i++) { */
/* 		double sign = 1.0; */
/* 		double norm1; */
/* 		H(test,h2vectors[i]); */
/* 		spinor_field_mul_add_assign_f(test,-sqrt(h2values[i]),h2vectors[i]); */
/* 		norm=spinor_field_sqnorm_f(test); */
/* 		H(test,h2vectors[i]); */
/* 		spinor_field_mul_add_assign_f(test,sqrt(h2values[i]),h2vectors[i]); */
/* 		norm1=spinor_field_sqnorm_f(test); */
/* 		if(norm1 < norm) { */
/* 			sign = -1.0; */
/* 			norm = norm1; */
/* 		} */
/* 		lprintf("DIRAC_EVA_ONEMASS",40,"H eigenvectors test [%d,%e] = %e\n",i,sign*sqrt(h2values[i]),norm); */
/* 	} */
/* #endif */

/* 	for(i = 0; i < nevt; i++) { */
/* 		H(test,h2vectors[i]); */
/* 		for(j = 0; j < i; j++) { */
/* 			smallH[i*nevt+j].re = smallH[j*nevt+i].re = spinor_field_prod_re_f(h2vectors[j],test); */
/* 			smallH[i*nevt+j].im = -spinor_field_prod_im_f(h2vectors[j],test); */
/* 			smallH[j*nevt+i].im = -smallH[i*nevt+j].im; */
/* 		} */
/* 		smallH[i*nevt+i].re = spinor_field_prod_re_f(h2vectors[j],test); */
/* 		smallH[i*nevt+i].im = 0.; */
/* 	} */


/* 	jacobi2(nevt, smallH, d, smallvectors); */


/* #ifdef TESTINGMODE	 */
/* 	for(i = 0; i < nevt; i++) { */
/* 		norm = 0.; */
/* 		for(j = 0; j < nevt; j++) { */
/* 			smalltest.re = -d[i]*smallvectors[j*nevt+i].re; */
/* 			smalltest.im = -d[i]*smallvectors[j*nevt+i].im; */
/* 			for(k = 0; k < nevt; k++) { */
/* 				smalltest.re += smallH[j*nevt+k].re*smallvectors[k*nevt+i].re */
/* 												-smallH[j*nevt+k].im*smallvectors[k*nevt+i].im; */
/* 				smalltest.im += smallH[j*nevt+k].re*smallvectors[k*nevt+i].im */
/* 												+smallH[j*nevt+k].im*smallvectors[k*nevt+i].re; */
/* 			} */
/* 			norm += smalltest.re*smalltest.re + smalltest.im*smalltest.im; */
/* 		} */
/* 		norm = sqrt(norm); */
/* 		lprintf("DIRAC_EVA_ONEMASS",50,"smallH eigenvectors test [%d,%e] = %e\n",i,d[i],norm); */
/* 	} */
/* 	 */
/* 	norm=0.; */
/* 	for(i = 0; i < nevt; i++) */
/* 	for(k = 0; k < nevt; k++) { */
/* 		prod.re = prod.im = 0.; */
/* 		for(j = 0; j < nevt; j++) { */
/* 			prod.re += smallvectors[j*nevt+i].re*smallvectors[j*nevt+k].re */
/* 						+smallvectors[j*nevt+i].im*smallvectors[j*nevt+k].im; */
/* 			prod.im += smallvectors[j*nevt+i].re*smallvectors[j*nevt+k].im */
/* 						-smallvectors[j*nevt+i].im*smallvectors[j*nevt+k].re; */
/* 		} */
/*     if(i == k) norm += (prod.re-1.)*(prod.re-1.); */
/*     else norm += prod.re*prod.re; */
/*     norm += prod.im*prod.im; */
/* 	} */
/* 	norm = sqrt(norm); */
/* 	lprintf("DIRAC_EVA_ONEMASS",50,"smallH eigenvectors norm test = %e\n",norm); */
/* #endif */


/* 	for(x = 0; x < VOLUME*sizeof(suNf_spinor)/sizeof(complex); x++) */
/* 		for(i = 0; i < nevt; i++) { */
/* 			ptmp1 = (complex*)(ev[i])+x; */
/* 			ptmp1->re = ptmp1->im = 0.0; */
/* 			for(j = 0; j < nevt; j++) { */
/* 				ptmp2 = (complex*)(h2vectors[j])+x; */
/* 				alpha = smallvectors + (j*nevt+i); */
/* 				ptmp1->re += alpha->re * ptmp2->re - alpha->im * ptmp2->im; */
/* 				ptmp1->im += alpha->re * ptmp2->im + alpha->im * ptmp2->re; */
/* 			} */
/* 		} */


/* 	for(i = 0; i < nevt; i++) { */
/* 		H(test,ev[i]); */
/* 		spinor_field_mul_add_assign_f(test,-d[i],ev[i]); */
/* 		norm=spinor_field_sqnorm_f(test); */
/* 		lprintf("DIRAC_EVA_ONEMASS",10,"H eigenvectors test [%d,%e] = %e\n",i,d[i],norm); */
/* 	} */
/* 	 */
/* #ifdef TESTINGMODE	 */
/*   norm = 0.; */
/*   for(i = 0; i < nevt; i++) */
/*   for(j = 0; j < nevt; j++) { */
/*     prod.re = spinor_field_prod_re_f(ev[i],ev[j]); */
/*     prod.im = spinor_field_prod_im_f(ev[i],ev[j]); */
/*     if(i == j) norm += (prod.re-1.)*(prod.re-1.); */
/*     else norm += prod.re*prod.re; */
/*     norm += prod.im*prod.im; */
/*   } */
/*   norm = sqrt(norm); */
/* 	lprintf("DIRAC_EVA_ONEMASS",40,"H eigenvectors norm test = %e\n",norm); */


/* 	for(i = 0; i < nevt; i++) { */
/* 		lprintf("DIRAC_EVA_ONEMASS",40,"H eigenvalues [%d] = %f\n",i,d[i]); */
/* 	} */
/* #endif */

/* 	free_field(test); */
/* 	free(h2values); */
/* 	for(i = 0; i < nevt; i++) */
/* 		free_field(h2vectors[i]); */
/* 	free_field(ws[0]); */
/* 	free_field(ws[1]); */
/* 	free_field(h2tmp); */
/* 	free(smallH); */
/* } */




void dirac_eva_onemass(int nev,int nevt,int kmax,
        int imax,double omega1,double omega2,double mass,
        suNf_spinor *ev[],double *d,int *status) {
	int i, j, maxh2iter, ie;
	unsigned int volume;
	suNf_spinor *ws[2];
	suNf_spinor *plus, *minus;
 	complex prod;

	int init;
/* Specifies whether all eigenvectors should be initialized */
/* (init=0), or only the last nevt-nev eigenvectors (init=1) */
/* or none of them (init!=0 or 1) */

	double ubnd;
/* Upper bound on the eigenvalues of the operator */

 	double norm;
	suNf_spinor *test=0;
	
	get_spinor_len(&volume);
	error(volume!=VOLUME,1,"dirac_eva_onemass","The spinors must be defined in each site");
	
	test = alloc_spinor_field_f(1);
	ws[0] = alloc_spinor_field_f(1);
	ws[1] = alloc_spinor_field_f(1);
	plus = alloc_spinor_field_f(1);
	minus = alloc_spinor_field_f(1);
	h2tmp = alloc_spinor_field_f(1);
	
	maxh2iter = max_H2(&ubnd,mass);
	init = 0;

	hmass = mass;
	
	ie = eva(VOLUME,nev,nevt,init,kmax,imax,ubnd,omega1,omega2,&H2,ws,ev,d,status);

	lprintf("DIRAC_EVA_ONEMASS",10,"Computed %d eigenvetors of H2.\n",nev);

#ifdef TESTINGMODE	
	for(i = 0; i < nevt; i++) {
		H2(test,ev[i]);
		spinor_field_mul_add_assign_f(test,-d[i],ev[i]);
		norm=spinor_field_sqnorm_f(test);
		lprintf("DIRAC_EVA_ONEMASS",40,"H2 eigenvectors test [%d,%e] = %e\n",i,d[i],norm);
	}

	norm = 0.;
	for(i = 0; i < nev; i++)
	for(j = 0; j < nev; j++) {
	   prod.re = spinor_field_prod_re_f(ev[i],ev[j]);
	   prod.im = spinor_field_prod_im_f(ev[i],ev[j]);
	   if(i == j) norm += (prod.re-1.)*(prod.re-1.);
	   else norm += prod.re*prod.re;
	   norm += prod.im*prod.im;
	}
	norm = sqrt(norm);
	lprintf("DIRAC_EVA_ONEMASS",40,"H2 eigenvectors norm test = %e\n",norm);

#endif

	for(i = 0; i < nev; i++) {
		double nplus, nminus;
		H(test,ev[i]);
		spinor_field_lc_f(minus, 1., ev[i], -1./sqrt(d[i]), test);
		spinor_field_lc_f(plus, 1., ev[i], 1./sqrt(d[i]), test);
		nplus = sqrt(spinor_field_sqnorm_f(plus));
		nminus = sqrt(spinor_field_sqnorm_f(minus));

		if(nplus > nminus) {
			spinor_field_mul_f(ev[i],1./nplus,plus);
			d[i] = sqrt(d[i]);
		} else {
			spinor_field_mul_f(ev[i],1./nminus,minus);
			d[i] = -sqrt(d[i]);
		}
		for(j = i+1; j < nev; j++) {
	 		prod.re = -spinor_field_prod_re_f(ev[i],ev[j]);
	 		prod.im = -spinor_field_prod_im_f(ev[i],ev[j]);
	 		spinor_field_mulc_add_assign_f(ev[j],prod,ev[i]);
	 		norm = sqrt(spinor_field_sqnorm_f(ev[j]));
			spinor_field_mul_f(ev[j],1./norm,ev[j]);
		}
	}


	for(i = 0; i < nev; i++) {
		H(test,ev[i]);
		spinor_field_mul_add_assign_f(test,-d[i],ev[i]);
		norm=spinor_field_sqnorm_f(test);
		lprintf("DIRAC_EVA_ONEMASS",10,"H eigenvectors test [%d,%e] = %e\n",i,d[i],norm);
	}
	
#ifdef TESTINGMODE
  double dnorm = 0.;
  norm = 0.;
  for(i = 0; i < nevt; i++)
  for(j = 0; j < nevt; j++) {
    prod.re = spinor_field_prod_re_f(ev[i],ev[j]);
    prod.im = spinor_field_prod_im_f(ev[i],ev[j]);
    if(i == j) {
    	norm += (prod.re-1.)*(prod.re-1.);
    	dnorm += (prod.re-1.)*(prod.re-1.) + prod.im*prod.im;
    } else norm += prod.re*prod.re;
    norm += prod.im*prod.im;
  }
  norm = sqrt(norm);
  dnorm = sqrt(dnorm);
	lprintf("DIRAC_EVA_ONEMASS",40,"H eigenvectors norm test = %e\n",norm);
	lprintf("DIRAC_EVA_ONEMASS",40,"H eigenvectors diagonal norm test = %e\n",dnorm);
	lprintf("DIRAC_EVA_ONEMASS",40,"H eigenvectors off-diagonal norm test = %e\n",sqrt(norm*norm-dnorm*dnorm));


	for(i = 0; i < nevt; i++) {
		lprintf("DIRAC_EVA_ONEMASS",40,"H eigenvalues [%d] = %f\n",i,d[i]);
	}
#endif

	free_field(test);
	free_field(plus);
	free_field(minus);
	free_field(ws[0]);
	free_field(ws[1]);
	free_field(h2tmp);
}



void dirac_eva(int nev,int nevt,int kmax,
        int imax,double omega1,double omega2,int n_masses,double *mass,
        suNf_spinor *ev[],double *d,int *status) {
/* ev[i*nevt+j], i<n_masses, j<nevt */
/* d[i*nevt+j], i<n_masses, j<nevt */
	int i;
	for(i = 0; i < n_masses; i++) {
		dirac_eva_onemass(nev,nevt,kmax,imax,omega1,omega2,mass[i],ev+i*nevt,d+i*nevt,status);
	}
}

