/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *   
* All rights reserved.                                                      * 
\***************************************************************************/

#include "inverters.h"
#include "linear_algebra.h"
#include "complex.h"
#include "memory.h"
#include "update.h"
#include "logger.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>

#define _print_par(c) \
  printf("[%d] " #c " = %e\n",cgiter,c)

/*
 * performs the multi-shifted MINRES inversion for hermitean matrices:
 * out[i] = (M-(par->shift[i]))^-1 in
 * returns the number of cg iterations done.
 */
static int MINRES_mshift_core(short int *flags,mshift_par *par, spinor_operator M, spinor_field *in, spinor_field *out){

  spinor_field **q1,**q2;
  spinor_field *p1, *p2, *Mp;
  spinor_field *sptmp, *memall;

  double alpha, beta, oldbeta,innorm2; 
  double *r, *s1, *s2, *c1, *c2, *rho1, *rho2, *rp;

  int i;
  int cgiter;
  unsigned int notconverged;

  /* fare qualche check sugli input */
  /* par->n deve essere almeno 1! */
  assert(par->n>0);
#ifndef CHECK_SPINOR_MATCHING
  for(i=0;i<par->n;++i)
   _TWO_SPINORS_MATCHING(in,&out[i]);
#endif

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

  memall = alloc_spinor_field_f(2*(par->n)+3,in->type);
  q1 = (spinor_field**)malloc(sizeof(spinor_field*)*par->n);
  q2 = (spinor_field**)malloc(sizeof(spinor_field*)*par->n);
  for(i=0; i<par->n; i++) {
    q1[i] = memall+i;
    q2[i] = memall+par->n+i;
  }
  p1 = memall+2*par->n;
  p2 = p1+1;
  Mp = p2+1;

  r = (double *)malloc(sizeof(double)*8*(par->n));
  s1 = r+(par->n);
  s2 = s1+(par->n);
  c1 = s2+(par->n);
  c2 = c1+(par->n);
  rho1 = c2+(par->n);
  rho2 = rho1+(par->n);
  rp = rho2+(par->n);

  /* init recursion */
  cgiter = 0;
  notconverged=par->n;

  if (par->n==1){ /* non multishift case */
    spinor_field_copy_f(p2, in);
    M(p1,&out[0]);
    spinor_field_mul_add_assign_f(p1,-par->shift[0],&out[0]);
    ++cgiter;
    spinor_field_sub_assign_f(p2,p1);
    /* use out[0] as trial solution */
  } else {
    spinor_field_copy_f(p2, in); /* trial solution = 0 */
  }
  innorm2=spinor_field_sqnorm_f(in);
  beta=sqrt(spinor_field_sqnorm_f(p2));
  spinor_field_mul_f(p2,1./beta,p2);
  spinor_field_zero_f(p1);  
  for (i=0; i<(par->n); ++i) {
    r[i]=rho2[i]=beta;
    rho1[i]=1.;
    c2[i]=-1.;
    rp[i]=s1[i]=s2[i]=c1[i]=0.;
    if(par->n!=1) spinor_field_zero_f(&out[i]);
    spinor_field_zero_f(q1[i]);
    spinor_field_zero_f(q2[i]);
    flags[i]=1;
  }

  /* cg recursion */
  do {
    ++cgiter;

    M(Mp,p2);
    spinor_field_mul_add_assign_f(Mp,-par->shift[0],p2);

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
	a=alpha-par->shift[i]+par->shift[0];
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
	spinor_field_mul_add_assign_f(&out[i],r[i]*c2[i]/k,q2[i]);

	/* update residuum */
	r[i]*=s2[i];

	if((r[i]*r[i])<par->err2*innorm2){
	  flags[i]=0;
	  --notconverged;
	  /* printf("[%d] converged at iter: %d\n",i,cgiter); */
	}

      }

    }    

  } while ((par->max_iter==0 || cgiter<par->max_iter) && notconverged);

  /* test results */
  for(i=0;i<par->n;++i){
    double norm;
    M(Mp,&out[i]);
    ++cgiter;
    spinor_field_mul_add_assign_f(Mp,-par->shift[i],&out[i]);
    spinor_field_sub_f(Mp,Mp,in);
    norm=spinor_field_sqnorm_f(Mp)/innorm2;
    flags[i]=1;
    if (fabs(norm)>par->err2){
      flags[i]=0;
      lprintf("INVERTER",30,"MINRES failed on vect %d: err2 = %1.8e > %1.8e\n",i,norm,par->err2);
    } else {
      lprintf("INVERTER",20,"MINRES inversion: err2 = %1.8e < %1.8e\n",norm,par->err2);
    } 
  }


  /* free memory */
  free_spinor_field_f(memall);
  free(q1);
  free(q2);
  free(r);

  /* return number of cg iter */
  return cgiter;
}

int MINRES_mshift(mshift_par *par, spinor_operator M, spinor_field *in, spinor_field *out){
  int iter,msiter;
  int i;
  mshift_par par_save=*par;
  short int *valid = malloc(sizeof(*valid)*(par->n));
	
  iter=MINRES_mshift_core(valid, par, M, in, out);
  msiter=iter;

  par->n=1;
  for(i=0;i<par_save.n;++i) {
    int rep=0;
    while (valid[i]==0) {
      par->shift=par_save.shift+i;
      iter+=MINRES_mshift_core(valid+i, par, M, in, out+i);
      if((++rep)%5==0)
	lprintf("INVERTER",-10,"MINRES_mshift recursion = %d (precision too high?)\n",rep);
    }
  }

  *par=par_save;

  free(valid);
	
  lprintf("INVERTER",10,"MINRES_mshift: MVM = %d/%d\n",msiter,iter);

  return iter;

}
