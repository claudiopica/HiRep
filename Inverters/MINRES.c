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
int MINRES(MINRES_par *par, spinor_operator M, suNf_spinor *in, suNf_spinor *out, suNf_spinor *trial){

  suNf_spinor *q1,*q2;
  suNf_spinor *p1, *p2, *Mp;
  suNf_spinor *sptmp, *memall;

  double alpha, beta, oldbeta, innorm2; 
  double r, s1, s2, c1, c2, rho1, rho2, rp;
  double d, h, k;

  int cgiter;
  unsigned short notconverged;
	unsigned int spinorlen;

  /* fare qualche check sugli input */
  /* par->n deve essere almeno 1! */
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
  memall = (suNf_spinor *)malloc(sizeof(suNf_spinor)*(5)*(spinorlen));
  q1=memall;
  q2= q1+(spinorlen);
  p1 = q2+(spinorlen);
  p2 = p1+(spinorlen);
  Mp = p2+(spinorlen);

  /* init recursion */
  cgiter = 0;
  notconverged=1;

  spinor_field_copy_f(p2, in);
  if(trial!=0) {
    M(p1,trial);
    spinor_field_sub_assign_f(p2,p1);
    if(out!=trial){
      spinor_field_copy_f(out,trial);
    }
    
  } else {
    spinor_field_zero_f(out);
  }

  innorm2=spinor_field_sqnorm_f(in);
  beta=sqrt(spinor_field_sqnorm_f(p2));
  spinor_field_mul_f(p2,1./beta,p2);
  spinor_field_zero_f(p1);  
  r=rho2=beta;
  rho1=1.;
  c2=-1.;
  rp=s1=s2=c1=0.;
  spinor_field_zero_f(q1);
  spinor_field_zero_f(q2);

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


    d=(alpha-rp*c1)*s2;
    h=oldbeta*s1;
    rp=-oldbeta*c1*s2-alpha*c2;
    k=sqrt(rp*rp+beta*beta);
    c1=c2;
    c2=rp/k;
    s1=s2;
    s2=beta/k;

    spinor_field_lc_f(Mp,-h/rho1,q1,-d/rho2,q2);
    sptmp=q1;
    q1=q2;
    q2=sptmp; /* swap q1[i]<->q2[i] */
    spinor_field_add_f(q2,p1,Mp);
	
    /* update rho */
    rho1=rho2;
    rho2=k;
	
    /* update solution */
    spinor_field_mul_add_assign_f(out,r*c2/k,q2);
	
    /* update residuum */
    r*=s2;

    if((r*r)<par->err2*innorm2){
      notconverged=0;
      /* printf("[%d] converged at iter: %d\n",i,cgiter); */
    }
	
	
  } while ((par->max_iter==0 || cgiter<par->max_iter) && notconverged);

  /* test results */
  {
    float norm;
    M(Mp,out);
    spinor_field_mul_add_assign_f(Mp,-1.0,in);
    norm=spinor_field_sqnorm_f(Mp)/spinor_field_sqnorm_f(in);
    if (fabs(norm)>par->err2)
      printf("MINRES Failed: err2 = %e\n",norm);
  }
#ifndef NDEBUG
  printf("MINRES: %d matrix multiplications\n",cgiter);
#endif
  
   
  /* free memory */
  free(memall);

  /* return number of cg iter */
  return cgiter;
}
