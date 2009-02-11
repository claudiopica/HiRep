/***************************************************************************\
 * Copyright (c) 2008, Claudio Pica                                         *   
 * All rights reserved.                                                     * 
 \***************************************************************************/

#include "inverters.h"
#include "linear_algebra.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "global.h"
#include "memory.h"
#include "utils.h"
#include "logger.h"
#include <assert.h>

/*
 * performs the multi-shifted CG inversion:
 * out[i] = (M-(par->shift[i]))^-1 in
 * returns the number of cg iterations done.
 */
static int cg_mshift_flt_core(short int *sflags, mshift_par *par, spinor_operator_flt M, spinor_field *in, spinor_field *out){

  spinor_field *k,*r,*Mk;
  spinor_field *p;
  spinor_field_flt *Mk_flt, *k_flt;
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
  Mk_flt = alloc_spinor_field_f_flt(2,in->type);
  k_flt = Mk_flt + 1;

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
    assign_sd2s(k_flt,&out[0]);
    M(Mk_flt,k_flt);
    assign_s2sd(Mk,Mk_flt);
    //M(Mk,&out[0]);
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
/*    sflags[i]=1; */
  }

  /* cg recursion */
  do {
    assign_sd2s(k_flt,k);
    M(Mk_flt,k_flt);
    assign_s2sd(Mk,Mk_flt);
    //M(Mk,k);
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
      if(delta*z3[i]*z3[i]>par->err2*innorm2) ++notconverged;
      if(sflags[i]){
        spinor_field_mul_f(&p[i],gamma*z3[i]*z3[i]/(z2[i]*z2[i]),&p[i]);
        spinor_field_mul_add_assign_f(&p[i],z3[i],r);
        z1[i]=z2[i];
        z2[i]=z3[i];
      }

      /* old code
      if(delta*z3[i]*z3[i]<par->err2*innorm2) sflags[i]=0;
      if(sflags[i]){
        notconverged++;
        spinor_field_mul_f(&p[i],gamma*z3[i]*z3[i]/(z2[i]*z2[i]),&p[i]);
        spinor_field_mul_add_assign_f(&p[i],z3[i],r);
        z1[i]=z2[i];
        z2[i]=z3[i];
      }
      */
    }

    /* Uncomment this to print cg recursion parameters 
       lprintf("CGTEST",0,"[ %d ] alpha=%e\n",cgiter,alpha);
       lprintf("CGTEST",0,"[ %d ] omega=%e\n",cgiter,omega);
       lprintf("CGTEST",0,"[ %d ] still runnning=%d\n",cgiter,notconverged);
       for (i=0;i<par->n;++i) lprintf("CGTEST",0,"z3[%d]=%e; ",i,z3[i]);
       lprintf("CGTEST",0,"\n[ %d ] gamma=%e\n",cgiter,gamma);
       lprintf("CGTEST",0,"[ %d ] delta=%e\n",cgiter,delta);
       */
       //lprintf("CGTEST",0,"[ %d ] still runnning=%d\n",cgiter,notconverged);

    ++cgiter;
  } while ((par->max_iter==0 || cgiter<par->max_iter) && notconverged);

  /* free memory */
  free_spinor_field(p);
  free_spinor_field_flt(Mk_flt);
  free(z1); free(z2); free(z3);

  /* return number of cg iter */
  return cgiter;
}


int cg_mshift_flt(mshift_par *par, spinor_operator M, spinor_operator_flt F, spinor_field *in, spinor_field *out){ 
  int siter=0,diter=0;
  int i;
  mshift_par local_par=*par;
  short int sflags[par->n];
  short int loc_flags[par->n];
  double norm[par->n];

  double innorm2;
  int notconverged;
  
  spinor_field *out_core;
  spinor_field *res, *res2, *tmp;

  /* check types */
   _TWO_SPINORS_MATCHING(in,out); 
  _ARRAY_SPINOR_MATCHING(out,i,par->n);

  /* allocate memory for single-precision solutions and residual vectors */
//  res_flt = alloc_spinor_field_f_flt(1+par->n,in->type);
//  out_flt = res_flt + 1;
  res = alloc_spinor_field_f(3+par->n,in->type);
  res2 = res + 1;
  tmp = res2 + 1;
  out_core = tmp + 1;

  /* set all flags to 1
   * set all out to zero execpt if par->n==1
   */
  for (i=0; i<(par->n); ++i) {
    sflags[i]=1; 
    if(par->n!=1) spinor_field_zero_f(&out[i]);
  }

  /* compute input norm2 */
  innorm2=spinor_field_sqnorm_f(in);

  /* begin external loop */
  do {
    int first=1;
    notconverged=0;
    local_par.err2=1.e-100; /* small number */

    for (i=0; i<par->n; ++i) {
      if (sflags[i]!=0) {
        /* compute residual vector */
        M(tmp,&out[i]); ++diter;
        spinor_field_sub_f(res2,in,tmp);
        spinor_field_mul_add_assign_f(res2,par->shift[i],&out[i]);
        /* test for convergence */
        norm[i]=spinor_field_sqnorm_f(res2);
        lprintf("CGDEBUG",20,"norm %d = %e relerr=%e\n",i,norm[i],norm[i]/innorm2);
        if (norm[i]<innorm2*par->err2) {
          /* this shift has reached convergence */
          sflags[i]=0;
          norm[i]/=innorm2;
          lprintf("CGDEBUG",20,"shift %d converged. (%d)\n",i,siter);
        } else {
          /* this shift has yet converged */
          ++notconverged;
          /* prepare input residual vector */
          if(first){
            /* we take the residual vector from the first non-converged shift */
            first=0;
            if(local_par.err2*innorm2<norm[i]) local_par.err2=norm[i]/innorm2;
            norm[i]=sqrt(norm[i]);
            spinor_field_mul_f(res,1./norm[i],res2); /* normalize input res vector */

          } else {
            /* take the scalar product of current residual vector with the 
             * normalized first one as residual norm */
            norm[i]=spinor_field_prod_im_f(res, res2);
            lprintf("CGDEBUG",20,"Im=%e check norm=%e\n",norm[i],spinor_field_sqnorm_f(res));
            //norm[i]=spinor_field_prod_re_f(res, res2);
            norm[i]=sqrt(spinor_field_sqnorm_f(res2));

          }
        }
      }
    }
            lprintf("CGDEBUG",20,"==============================================\n");


    if (notconverged) {
      /* Do single precision inversion */
      //assign_sd2s(res_flt,res); 
      //local_par.err2=1.e-14; /* can do something smarter ? */
      //if (par->err2>local_par.err2) local_par.err2=par->err2;
      local_par.err2=par->err2/local_par.err2/2.;
      if(local_par.err2<1.e-8) local_par.err2=1.e-8;
      lprintf("CGDEBUG",20,"err2 = %e\n",local_par.err2);
      for(i=0;i<par->n;++i) loc_flags[i]=sflags[i]; /* set core flags */
      siter+=cg_mshift_flt_core(loc_flags, &local_par, F, res, out_core); /* save single precision iterations */

      /* accumulate solution in double precision */
      for(i=0;i<par->n; ++i) {
        if(sflags[i]!=0){
          //assign_s2sd(res2,&out_flt[i]);
          spinor_field_mul_add_assign_f(&out[i],norm[i],&out_core[i]);
        }
      }
    }

  } while(notconverged);

  
  for (i=0;i<par->n;++i){
    lprintf("INVERTER",20,"CG inversion: err2 = %1.8e < %1.8e\n",norm[i],par->err2);
  }

  //free_spinor_field_flt(res_flt);
  free_spinor_field(res);

  lprintf("INVERTER",10,"CG_mshift: MVM = %d (single) - %d (double)\n",siter,diter);

  return siter+diter;
}

