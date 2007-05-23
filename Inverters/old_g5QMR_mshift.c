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
int g5QMR_mshift(QMR_mshift_par *par, spinor_operator M, suNf_spinor *in, suNf_spinor **out){

  suNf_spinor **q1,**q2;
  suNf_spinor *p1, *p2, *Mp;
  suNf_spinor *sptmp, *memall;

  double alpha, beta, oldbeta, delta, olddelta, rho, norm2in; 
  double *r, *s1, *s2, *c1, *c2;

  int i;
  int cgiter;
  unsigned int notconverged;
  
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
  q1 = (suNf_spinor **)malloc(sizeof(suNf_spinor*)*2*(par->n));
  q2 = q1+(par->n);
  q1[0] = (suNf_spinor *)malloc(sizeof(suNf_spinor)*(2*(par->n)+3)*(par->spinorlen));
  q2[0] = q1[0]+(par->n)*(par->spinorlen);
  memall=q1[0];
  for (i=1; i<(par->n); ++i) {
    q1[i] = q1[i-1]+(par->spinorlen);
    q2[i] = q2[i-1]+(par->spinorlen);
  }
  p1 = q2[par->n-1]+(par->spinorlen);
  p2 = p1+(par->spinorlen);
  Mp = p2+(par->spinorlen);

  r = (double *)malloc(sizeof(double)*5*(par->n));
  s1 = r+(par->n);
  s2 = s1+(par->n);
  c1 = s2+(par->n);
  c2 = c1+(par->n);

  flags=(unsigned short *)malloc(sizeof(unsigned short)*(par->n));

  /* init recursion */
  cgiter = 0;
  notconverged=par->n;

  spinor_field_copy_f(p2, in); /* trial solution = 0 */
  norm2in=rho=sqrt(spinor_field_sqnorm_f(p2));
  spinor_field_mul_f(p2,1./rho,p2);
  spinor_field_zero_f(p1);  
  for (i=0; i<(par->n); ++i) {
    r[i]=rho;
    olddelta=c2[i]=c1[i]=1.;
    s1[i]=s2[i]=0.;
    spinor_field_zero_f(out[i]);
    spinor_field_zero_f(q1[i]);
    spinor_field_zero_f(q2[i]);
    flags[i]=1;
  }

  delta = spinor_field_g5_prod_re_f(p2,p2);
  beta = rho*delta;

  /* cg recursion */
  do {
    ++cgiter;

    M(Mp,p2);

    /* compute alpha */
    alpha = spinor_field_g5_prod_re_f(p2,Mp)/delta;

    /* update p1, p2 */
    spinor_field_mul_add_assign_f(Mp,-beta,p1);
    spinor_field_mul_add_assign_f(Mp,-alpha,p2);
    sptmp=p1;
    p1=p2;
    p2=Mp;
    Mp=sptmp;

    /* normalize p2 */
    rho=sqrt(spinor_field_sqnorm_f(p2));
    spinor_field_mul_f(p2,1./rho,p2);

    /* update delta */
    olddelta=delta;
    delta = spinor_field_g5_prod_re_f(p2,p2);

    /* update beta */
    oldbeta=beta;
    beta=rho*delta/olddelta;
    
    for (i=0; i<(par->n); ++i) { /* update solutions */
      if (flags[i]) {
	double a, t, e, et, dt, n, g;
	a=(i==0)?alpha:(alpha-par->shift[i-1]);
	t=s1[i]*oldbeta;
	et=c1[i]*oldbeta;
	e=c2[i]*et+s2[i]*a;
	dt=c2[i]*a-s2[i]*et;
	_print_par(dt);
	n=sqrt(dt*dt+beta*beta);
	c1[i]=c2[i];
	c2[i]=fabs(dt)/n;
	s1[i]=s2[i];
	s2[i]=(fabs(dt)<1.e-8)?1.:c2[i]*beta/dt; /* test it ! */
	g=1./(c2[i]*dt+s2[i]*beta);

	/* update q */
	spinor_field_lc_f(Mp,-t*g,q1[i],-e*g,q2[i]);
	spinor_field_mul_add_assign_f(Mp,g,p1);
	sptmp=q1[i];
	q1[i]=q2[i];
	q2[i]=Mp;
	Mp=sptmp;/* swap q1[i]<-q2[i]<-Mp and Mp point to q1 */
      
	/* update solution */
	spinor_field_mul_add_assign_f(out[i],c2[i]*r[i],q2[i]);
	
	/* update residuum */
	r[i]*=s2[i];

	if((r[i]*r[i])<par->err2*norm2in){
	  flags[i]=0;
	  --notconverged;
	  printf("[%d] converged at iter: %d\n",i,cgiter);
	}
	
	_print_par(r[i]);	
	
      }    
    }
    
    printf("notconverged=%d\n",notconverged);
    

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
    norm=spinor_field_sqnorm_f(Mp);
    if (fabs(norm)>5.*par->err2)
      printf("g5QMR Failed: err2[%d] = %e\n",i,norm);
  }
  printf("g5QMR: %d matrix multiplications\n",cgiter);
#endif
   
  /* free memory */
  free(memall);
  free(q1);
  free(r);

  free(flags);

  /* return number of cg iter */
  return cgiter;
}
