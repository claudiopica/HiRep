#include "inverters.h"
#include "linear_algebra.h"
#include "complex.h"
#include "malloc.h"
#include "update.h"
#include <stdio.h>
#include <math.h>
#include <assert.h>

#define _print_par(c) \
  printf("[%d] " #c " = %e\n",cgiter,c)

/*
 * performs the multi-shifted QMR inversion for g5-hermitean matrices:
 * out[i] = (M-(par->shift[i]))^-1 in
 * returns the number of cg iterations done.
 */
int MINRES_mshift(mshift_par *par, spinor_operator M, suNf_spinor *in, suNf_spinor **out){

  suNf_spinor **q1,**q2;
  suNf_spinor *p1, *p2, *Mp;
  suNf_spinor *sptmp, *memall;

  double alpha, beta, oldbeta,innorm2; 
  double *r, *s1, *s2, *c1, *c2, *rho1, *rho2, *rp;

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
  q1 = (suNf_spinor **)malloc(sizeof(suNf_spinor*)*2*(par->n));
  q2 = q1+(par->n);
  q1[0] = (suNf_spinor *)malloc(sizeof(suNf_spinor)*(2*(par->n)+3)*(spinorlen));
  q2[0] = q1[0]+(par->n)*(spinorlen);
  memall=q1[0];
  for (i=1; i<(par->n); ++i) {
    q1[i] = q1[i-1]+(spinorlen);
    q2[i] = q2[i-1]+(spinorlen);
  }
  p1 = q2[par->n-1]+(spinorlen);
  p2 = p1+(spinorlen);
  Mp = p2+(spinorlen);

  r = (double *)malloc(sizeof(double)*8*(par->n));
  s1 = r+(par->n);
  s2 = s1+(par->n);
  c1 = s2+(par->n);
  c2 = c1+(par->n);
  rho1 = c2+(par->n);
  rho2 = rho1+(par->n);
  rp = rho2+(par->n);

  flags=(unsigned short *)malloc(sizeof(unsigned short)*(par->n));

  /* init recursion */
  cgiter = 0;
  notconverged=par->n;

  spinor_field_copy_f(p2, in); /* trial solution = 0 */
  innorm2=spinor_field_sqnorm_f(p2);
  beta=sqrt(innorm2);
  spinor_field_mul_f(p2,1./beta,p2);
  spinor_field_zero_f(p1);  
  for (i=0; i<(par->n); ++i) {
    r[i]=rho2[i]=beta;
    rho1[i]=1.;
    c2[i]=-1.;
    rp[i]=s1[i]=s2[i]=c1[i]=0.;
    spinor_field_zero_f(out[i]);
    spinor_field_zero_f(q1[i]);
    spinor_field_zero_f(q2[i]);
    flags[i]=1;
  }

  /* cg recursion */
  do {
    ++cgiter;

    M(Mp,p2);
    
    /* compute alpha */
    alpha = spinor_field_prod_re_f(Mp,p2);

    /* update p1, p2 */
    spinor_field_mul_add_assign_f(Mp,-beta,p1);
    spinor_field_mul_add_assign_f(Mp,-alpha,p2);
    sptmp=p1;
    p1=p2;
    p2=Mp;
    Mp=sptmp;

    /* update beta */
    oldbeta=beta;
    beta=sqrt(spinor_field_sqnorm_f(p2));
    
    /* normalize p2 */
    spinor_field_mul_f(p2,1./beta,p2);

    for (i=0; i<(par->n); ++i) { /* update solutions */
      if(flags[i]) {
	double d, h, a, k;
	a=(i==0)?alpha:(alpha-par->shift[i-1]);
	d=(a-rp[i]*c1[i])*s2[i];
	h=oldbeta*s1[i];
	rp[i]=-oldbeta*c1[i]*s2[i]-a*c2[i];
	k=sqrt(rp[i]*rp[i]+beta*beta);
	c1[i]=c2[i];
	c2[i]=rp[i]/k;
	s1[i]=s2[i];
	s2[i]=beta/k;

	spinor_field_lc_f(Mp,-h/rho1[i],q1[i],-d/rho2[i],q2[i]);
	sptmp=q1[i];
	q1[i]=q2[i];
	q2[i]=sptmp; /* swap q1[i]<->q2[i] */
	spinor_field_add_f(q2[i],p1,Mp);
	
	/* update rho */
	rho1[i]=rho2[i];
	rho2[i]=k;
	
	/* update solution */
	spinor_field_mul_add_assign_f(out[i],r[i]*c2[i]/k,q2[i]);
	
	/* update residuum */
	r[i]*=s2[i];

	if((r[i]*r[i])<par->err2*innorm2){
	  flags[i]=0;
	  --notconverged;
	  /* printf("[%d] converged at iter: %d\n",i,cgiter); */
	}
	
	/* _print_par(r[i]); */
      }
	
    }    


  } while ((par->max_iter==0 || cgiter<par->max_iter) && notconverged);

  /* test results */
#ifndef NDEBUG
  for(i=0;i<par->n;++i){
    float norm;
    M(Mp,out[i]);
    if(i!=0) {
      spinor_field_mul_add_assign_f(Mp,-par->shift[i-1],out[i]);
    }
    spinor_field_mul_add_assign_f(Mp,-1.0,in);
    norm=spinor_field_sqnorm_f(Mp)/spinor_field_sqnorm_f(in);
    if (fabs(norm)>par->err2)
      printf("MINRES Failed: err2[%d] = %e\n",i,norm);
  }
  printf("MINRES: %d matrix multiplications\n",cgiter);
#endif
  
   
  /* free memory */
  free(memall);
  free(q1);
  free(r);

  free(flags);

  /* return number of cg iter */
  return cgiter;
}
