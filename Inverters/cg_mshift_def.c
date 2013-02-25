/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *   
* All rights reserved.                                                      * 
\***************************************************************************/

#include "inverters.h"
#include "linear_algebra.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "global.h"
#include "memory.h"
#include "logger.h"
#include <assert.h>

/*
 * performs the multi-shifted CG inversion:
 * out[i] = (M-(par->shift[i]))^-1 in
 * returns the number of cg iterations done.
 */
static int cg_mshift_core(short int *sflags, mshift_par *par, spinor_operator M, spinor_field *in, spinor_field *out) {

  spinor_field *k,*r,*Mk;
  spinor_field *p;
  double omega, oldomega, gamma;
  double alpha, lambda, delta;
  double innorm2;
  double *z1, *z2, *z3;

  int i;
  int cgiter;
  unsigned short notconverged;

  /* fare qualche check sugli input */
  assert(par->n>0);
#ifndef CHECK_SPINOR_MATCHING
  for(i=0;i<par->n;++i)
    _TWO_SPINORS_MATCHING(in,&out[i]);
#endif

  /* allocate spinors fields and aux real variables */
  p = alloc_spinor_field_f(3+par->n,in->type);
  k=p+par->n;
  r=k+1;
  Mk=r+1;

  z1 = malloc(sizeof(*z1)*(par->n));
  z2 = malloc(sizeof(*z2)*(par->n));
  z3 = malloc(sizeof(*z3)*(par->n));

  /* init recursion */
  cgiter = 0;
  omega = 1.;
  gamma = 0.;
  innorm2=spinor_field_sqnorm_f(in);
  if(par->n==1) { /* non multishift case */
    /* use out[0] as initial guess */
    M(Mk,&out[0]);
    ++cgiter;
    spinor_field_mul_add_assign_f(Mk,-par->shift[0],&out[0]);
    spinor_field_sub_f(r,in,Mk);

  } else { /* initial guess = 0 for multishift */
    spinor_field_copy_f(r, in);
  }

  spinor_field_copy_f(k, r);
  delta=spinor_field_sqnorm_f(r);
  for (i=0; i<(par->n); ++i) {
    z1[i]=z2[i]=1.;
    spinor_field_copy_f(&p[i], r);
    if(par->n!=1) spinor_field_zero_f(&out[i]);
    sflags[i]=1;
  }

  /* cg recursion */
  do {
    M(Mk,k);
    alpha = spinor_field_prod_re_f(k,Mk);
    oldomega = omega;
    omega = - delta/alpha;
    for (i=0; i<(par->n); ++i) {
      if(sflags[i]) {
        z3[i] = oldomega*z1[i]*z2[i]/(omega*gamma*(z1[i]-z2[i])+z1[i]*oldomega*(1.+par->shift[i]*omega));
        spinor_field_mul_add_assign_f(&out[i],-omega*z3[i]/z2[i],&p[i]);
      }
    }
    spinor_field_mul_add_assign_f(r,omega,Mk);
    lambda=spinor_field_sqnorm_f(r);
    gamma=lambda/delta;
    delta=lambda;

    spinor_field_mul_f(k,gamma,k);
    spinor_field_add_assign_f(k,r);
    notconverged=0; /* assume that all vectors have converged */
    for (i=0; i<(par->n); ++i) {
      /* check convergence of vectors */
      if(delta*z3[i]*z3[i]<par->err2*innorm2) sflags[i]=0; 
      if(sflags[i]){
        notconverged++;
        spinor_field_mul_f(&p[i],gamma*z3[i]*z3[i]/(z2[i]*z2[i]),&p[i]);
        spinor_field_mul_add_assign_f(&p[i],z3[i],r);
        z1[i]=z2[i];
        z2[i]=z3[i];
      }
    }

    /* Uncomment this to print cg recursion parameters 
       lprintf("CGTEST",0,"[ %d ] alpha=%e\n",cgiter,alpha);
       lprintf("CGTEST",0,"[ %d ] omega=%e\n",cgiter,omega);
       lprintf("CGTEST",0,"[ %d ] still runnning=%d\n",cgiter,notconverged);
       for (i=0;i<par->n;++i) lprintf("CGTEST",0,"z3[%d]=%e; ",i,z3[i]);
       lprintf("CGTEST",0,"\n[ %d ] gamma=%e\n",cgiter,gamma);
       lprintf("CGTEST",0,"[ %d ] delta=%e\n",cgiter,delta);
       */

    ++cgiter;
  } while ((par->max_iter==0 || cgiter<par->max_iter) && notconverged);

  /* test results */
  for(i=0;i<par->n;++i){
    double norm;
    M(Mk,&out[i]);
    ++cgiter;
    spinor_field_mul_add_assign_f(Mk,-par->shift[i],&out[i]);
    spinor_field_sub_f(Mk,Mk,in);
    norm=spinor_field_sqnorm_f(Mk)/spinor_field_sqnorm_f(in);
    sflags[i]=1;
    if (fabs(norm)>par->err2){
      sflags[i]=0;
      lprintf("INVERTER",30,"CG failed on vect %d: err2 = %1.8e > %1.8e\n",i,norm,par->err2);
    } else {
      lprintf("INVERTER",20,"CG inversion: err2 = %1.8e < %1.8e\n",norm,par->err2);
    }
  }

  /* free memory */
  free_spinor_field_f(p);
  free(z1); free(z2); free(z3);

  /* return number of cg iter */
  return cgiter;
}

/* We assume to have an exact eigenspace and use the following strategy:
 * 1) split the input vector into 2 pieces: one orthogonal and one parallel to the eigenspace
 * defined by the projector P;
 * 2) compute the inverse in the eigenspace defined by P using Pinv;
 * 3) compute the inverse in the orthogonal eigenspace using cg;
 * add together the two pieces and again run cg to be sure of convergence.
 */
int cg_mshift_def(mshift_par *par, spinor_operator M, spinor_operator P, spinor_operator_m Pinv, spinor_field *in, spinor_field *out){ 
  int cgiter,msiter;
  int i;
  mshift_par par_save=*par;
  short int sflags[par->n];
  spinor_field *p1, *p2;

  /* allocate spinors fields and aux real variables */
  p1 = alloc_spinor_field_f(2,in->type);
  p2 = p1 + 1;

  /* split input vector */
  P(p1,in); /* orthogonal part */
  spinor_field_sub_f(p2,in,p1); /* parallel part */

  /* invert orthogonal part */
  cgiter=cg_mshift_core(sflags, par, M, p1, out);
  msiter=cgiter; /* save multishift iterations for logging */

  /* add back parallel part with Pinv */
  for(i=0;i<par_save.n;++i) {
    Pinv(p1,p2,par->shift[i]); /* store temporary res in p1 */
    spinor_field_add_assign_f(out+i,p1);
  }


  /* check result */
  par->n=1;
  for(i=0;i<par_save.n;++i) {
    int rep=0;
    sflags[i]=0;
    while (sflags[i]==0) {
      par->shift=par_save.shift+i;
      cgiter+=cg_mshift_core(sflags+i, par, M, in, out+i);
      if((++rep)%5==0)
        lprintf("INVERTER",-10,"CG_mshift recursion = %d (precision too high?)\n",rep);
    }
  }

  *par=par_save;

  lprintf("INVERTER",10,"CG_mshift: MVM = %d/%d\n",msiter,cgiter);

  /* free memory */
  free_spinor_field_f(p1);

  return cgiter;
}
