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
static int cg_mshift_flt_core(short int *sflags, mshift_par *par, spinor_operator_flt M, spinor_field_flt *in, spinor_field_flt *out){

  spinor_field_flt *k,*r,*Mk;
  spinor_field_flt *p;
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
  p = alloc_spinor_field_f_flt(3+par->n,in->type);
#ifdef WITH_GPU
  alloc_spinor_field_f_flt_gpu(3+par->n,p);
#endif //WITH_GPU
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
  innorm2=spinor_field_sqnorm_f_flt(in);
  if(par->n==1) { /* non multishift case */
    /* use out[0] as initial guess */
    M(Mk,&out[0]);
    ++cgiter;
    spinor_field_mul_add_assign_f_flt(Mk,-par->shift[0],&out[0]);
    spinor_field_sub_f_flt(r,in,Mk);

  } else { /* initial guess = 0 for multishift */
    spinor_field_copy_f_flt(r, in);
  }
  spinor_field_copy_f_flt(k, r);
  delta=spinor_field_sqnorm_f_flt(r);
  for (i=0; i<(par->n); ++i) {
    z1[i]=z2[i]=1.;
    spinor_field_copy_f_flt(&p[i], r);
    if(par->n!=1) spinor_field_zero_f_flt(&out[i]);
/*    sflags[i]=1; */
  }

  /* cg recursion */
  do {
    M(Mk,k);
    alpha = spinor_field_prod_re_f_flt(k,Mk);
    oldomega = omega;
    omega = - delta/alpha;
    for (i=0; i<(par->n); ++i) {
      if(sflags[i]) {
        z3[i] = oldomega*z1[i]*z2[i]/(omega*gamma*(z1[i]-z2[i])+z1[i]*oldomega*(1.+par->shift[i]*omega));
        spinor_field_mul_add_assign_f_flt(&out[i],-omega*z3[i]/z2[i],&p[i]);
      }
    }
    spinor_field_mul_add_assign_f_flt(r,omega,Mk);
    lambda=spinor_field_sqnorm_f_flt(r);
    gamma=lambda/delta;
    delta=lambda;

    spinor_field_mul_f_flt(k,gamma,k);
    spinor_field_add_assign_f_flt(k,r);
    notconverged=0; /* assume that all vectors have converged */
    for (i=0; i<(par->n); ++i) {
      /* check convergence of vectors */
      if(delta*z3[i]*z3[i]>par->err2*innorm2) ++notconverged;
      if(sflags[i]){
        spinor_field_mul_f_flt(&p[i],gamma*z3[i]*z3[i]/(z2[i]*z2[i]),&p[i]);
        spinor_field_mul_add_assign_f_flt(&p[i],z3[i],r);
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

  /* free memory */
#ifdef WITH_GPU
  free_spinor_field_flt_gpu(p);
#endif //WITH_GPU
  free_spinor_field_flt(p);

  free(z1); free(z2); free(z3);

  /* return number of cg iter */
  return cgiter;
}

static spinor_field_flt *p_glb;

/*
 * performs shifted CG inversion:
 * out = (M-par)^-1 in
 * returns the number of cg iterations done.
 */
static int cg_flt_core(double shift,  mshift_par *par, spinor_operator_flt M, spinor_field_flt *in, spinor_field_flt *out){
  spinor_field_flt *k,*r,*Mk;
  spinor_field_flt *p;
  double omega, oldomega, gamma;
  double alpha, lambda, delta;
  double innorm2;
  double z1, z2, z3;

  int i;
  int cgiter;
  unsigned short notconverged;
  static int first=1;

  /* fare qualche check sugli input */
#ifndef CHECK_SPINOR_MATCHING
  _TWO_SPINORS_MATCHING(in,out);
#endif

  if (first){
    /* allocate spinors fields and aux real variables */
    p_glb = alloc_spinor_field_f_flt(4,in->type);
#ifdef WITH_GPU
    alloc_spinor_field_f_flt_gpu(4,p_glb);
#endif //WITH_GPU
    first = 0;
  }

  p=p_glb;
  k=p+1;
  r=k+1;
  Mk=r+1;
  /* init recursion */
  cgiter = 0;
  omega = 1.;
  gamma = 0.;
  innorm2=spinor_field_sqnorm_f_flt(in);
  /* use out[0] as initial guess */
  M(Mk,out);
  ++cgiter;
  spinor_field_mul_add_assign_f_flt(Mk,-shift,out);
  spinor_field_sub_f_flt(r,in,Mk);

  spinor_field_copy_f_flt(k, r);
  delta=spinor_field_sqnorm_f_flt(r);

  z1=z2=1.;
  spinor_field_copy_f_flt(p, r);

  /* cg recursion */
  do {
    M(Mk,k);
    alpha = spinor_field_prod_re_f_flt(k,Mk);
    oldomega = omega;
    omega = - delta/alpha;

    z3 = oldomega*z1*z2/(omega*gamma*(z1-z2)+z1*oldomega*(1.+shift*omega));
    spinor_field_mul_add_assign_f_flt(out,-omega*z3/z2,p);

    spinor_field_mul_add_assign_f_flt(r,omega,Mk);
    lambda=spinor_field_sqnorm_f_flt(r);
    gamma=lambda/delta;
    delta=lambda;

    spinor_field_mul_f_flt(k,gamma,k);
    spinor_field_add_assign_f_flt(k,r);
    
    /* check convergence of vectors */
    notconverged = (delta*z3*z3 > par->err2*innorm2);

    spinor_field_mul_f_flt(p,gamma*z3*z3/(z2*z2),p);
    spinor_field_mul_add_assign_f_flt(p,z3,r);
    z1=z2;
    z2=z3;

    //    Uncomment this to print cg recursion parameters 
    /*       lprintf("CGTEST",0,"[ %d ] alpha=%e\n",cgiter,alpha);
       lprintf("CGTEST",0,"[ %d ] omega=%e\n",cgiter,omega);
       lprintf("CGTEST",0,"[ %d ] still runnning=%d\n",cgiter,notconverged);
       lprintf("CGTEST",0,"z3=%e; ",i,z3);
       lprintf("CGTEST",0,"\n[ %d ] gamma=%e\n",cgiter,gamma);
       lprintf("CGTEST",0,"[ %d ] delta=%e\n",cgiter,delta);*/
       
    ++cgiter;
  } while ((par->max_iter==0 || cgiter<par->max_iter) && notconverged);

  /* return number of cg iter */
  return cgiter;
}


/* cg mshift with single precision acceleration */
/* results are still accurate to double precision! */
int cg_mshift_flt(mshift_par *par, spinor_operator M, spinor_operator_flt F, spinor_field *in, spinor_field *out){ 
  int siter=0,diter=0;
  int i;
  mshift_par local_par=*par;
  short int sflags[par->n]; /* this is used in the core routine to select which shifts to update */
  short int loc_flags[par->n];
  double norm[par->n];

  double innorm2;
  int notconverged;
  
  spinor_field_flt *out_flt, *res_flt,*tmp_flt;
  spinor_field *res, *res2, *tmp;

  /* check types */
  _TWO_SPINORS_MATCHING(in,out); 
  _ARRAY_SPINOR_MATCHING(out,i,par->n);

  /* allocate memory for single-precision solutions and residual vectors */
  res_flt = alloc_spinor_field_f_flt(2+par->n,in->type);
#ifdef WITH_GPU
  alloc_spinor_field_f_flt_gpu(2+par->n,res_flt);
#endif //WITH_GPU

  out_flt = res_flt + 1;
  tmp_flt = out_flt + par->n;
  res = alloc_spinor_field_f(3,in->type);
#ifdef WITH_GPU
  alloc_spinor_field_f_gpu(3,res);
#endif //WITH_GPU

  res2 = res + 1;
  tmp = res2 + 1;

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
          /* this shift has not yet converged */
          ++notconverged;
          /* prepare input residual vector */
          if(first){
            /* we take the residual vector from the first non-converged shift */
            first=0;
            local_par.err2=norm[i]/innorm2;
            norm[i]=sqrt(norm[i]);
            spinor_field_mul_f(res,1./norm[i],res2); /* normalize input res vector */
          } else {
            /* take the scalar product of current residual vector with the 
             * normalized first one as residual norm */
            norm[i]=spinor_field_prod_im_f(res, res2);
            lprintf("CGDEBUG",20,"Im=%e check norm=%e\n",norm[i],spinor_field_sqnorm_f(res));
            norm[i]=spinor_field_prod_re_f(res, res2);
            /* norm[i]=sqrt(spinor_field_sqnorm_f(res2)); */

          }
        }
      }
    }
    lprintf("CGDEBUG",20,"==============================================\n");


    if (notconverged) {
      /* Do single precision inversion */
#ifdef WITH_GPU
      assign_sd2s_gpu(res_flt,res); 
#else
      assign_sd2s(res_flt,res); 
#endif
      local_par.err2=par->err2/local_par.err2*0.9;
#define MAX_PREC 1.e-10
      if(local_par.err2<MAX_PREC) local_par.err2=MAX_PREC;
#undef MAX_PREC
      /* local_par.err2=1.e-10; */
      lprintf("CGDEBUG",20,"err2 = %e\n",local_par.err2);
      /* should not be needed since core do not modify flags */
      for(i=0;i<par->n;++i) loc_flags[i]=sflags[i]; /* set core flags */ 
      siter+=cg_mshift_flt_core(loc_flags, &local_par, F, res_flt, out_flt); /* save single precision iterations */
      
      /* accumulate solution in double precision */
      for(i=0;i<par->n; ++i) {
        if(sflags[i]!=0){
#ifdef WITH_GPU
          assign_s2sd_gpu(res2,&out_flt[i]);
#else
          assign_s2sd(res2,&out_flt[i]);
#endif 
          spinor_field_mul_add_assign_f(&out[i],norm[i],res2);
        }
      }
    }
  } while(notconverged);

  
  for (i=0;i<par->n;++i){
    lprintf("INVERTER",20,"CG inversion: err2 = %1.8e < %1.8e\n",norm[i],par->err2);
  }

#ifdef WITH_GPU
  free_spinor_field_flt_gpu(res_flt);
#endif //WITH_GPU
  free_spinor_field_flt(res_flt);

#ifdef WITH_GPU
  free_spinor_field_gpu(res);
#endif //WITH_GPU
  free_spinor_field(res);

  lprintf("INVERTER",10,"CG_mshift: MVM = %d (single) - %d (double)\n",siter,diter);

  return siter+diter;
}

/* cg mshift with single precision acceleration */
/* results are still accurate to double precision! */
int cg_mshift_flt2(mshift_par *par, spinor_operator M, spinor_operator_flt F, spinor_field *in, spinor_field *out){ 
  int siter=0,diter=0;
  int i,j;
  mshift_par local_par=*par;
  short int sflags[par->n]; /* this is used in the core routine to select which shifts to update */
  //  short int loc_flags[par->n];
  double norm[par->n];

  double innorm2;
  int notconverged;
  
  spinor_field_flt *out_flt, *res_flt,*tmp_flt;
  spinor_field *res, *tmp;

  /* check types */
  _TWO_SPINORS_MATCHING(in,out); 
  _ARRAY_SPINOR_MATCHING(out,i,par->n);

  /* allocate memory for single-precision solutions and residual vectors */
  res_flt = alloc_spinor_field_f_flt(2+par->n,in->type);
#ifdef WITH_GPU
  alloc_spinor_field_f_flt_gpu(2+par->n,res_flt);
#endif //WITH_GPU

  out_flt = res_flt + 1;
  tmp_flt = out_flt + par->n;
  res = alloc_spinor_field_f(2,in->type);
#ifdef WITH_GPU
  alloc_spinor_field_f_gpu(2,res);
#endif //WITH_GPU

  tmp = res + 1;

  /* compute input norm2 */
  innorm2=spinor_field_sqnorm_f(in);
  
  /* set all flags to 1
   * set all out to zero execpt if par->n==1
   */
  for (i=0; i<(par->n); ++i) {
    sflags[i]=1;
    if(par->n!=1) spinor_field_zero_f(&out[i]);
  }

#ifdef WITH_GPU
  assign_sd2s_gpu(res_flt,in); 
#else
  assign_sd2s(res_flt,in); 
#endif

#define MAX_PREC 1.e-10
  if(local_par.err2<MAX_PREC) local_par.err2=MAX_PREC;
#undef MAX_PREC
  siter=cg_mshift_flt_core(sflags, &local_par, F, res_flt, out_flt); /* save single precision iterations */

  for (i=0;i<(par->n); ++i){
#ifdef WITH_GPU
    assign_s2sd_gpu(&out[i],&out_flt[i]);
#else
    assign_s2sd(&out[i],&out_flt[i]); 
#endif
  }

  for (i=0;i<(par->n); ++i){
    for (;;){
      M(tmp,&out[i]); ++diter;
      spinor_field_sub_f(res,in,tmp);
      spinor_field_mul_add_assign_f(res,par->shift[i],&out[i]);
      /* test for convergence */
      norm[i]=spinor_field_sqnorm_f(res);
      if (norm[i]<innorm2*par->err2) {
	/* this shift has reached convergence */
	norm[i]/=innorm2;
	lprintf("CGDEBUG",20,"shift %d converged. err: %1.10g (%d)\n",i,norm[i],siter);
	break;
      }
      lprintf("CGDEBUG",20,"shift %d, err: %1.10g, acc: %1.10g (%d)\n",i,norm[i]/innorm2,par->err2/local_par.err2,siter);

#ifdef WITH_GPU
      assign_sd2s_gpu(res_flt,res); 
#else
      assign_sd2s(res_flt,res); 
#endif
      local_par.err2=par->err2*innorm2/norm[i]*0.9;
#define MAX_PREC 1.e-10
      if(local_par.err2<MAX_PREC) local_par.err2=MAX_PREC;
#undef MAX_PREC
      lprintf("CGDEBUG",20,"\nerr2: %1.10g\n",local_par.err2);      
      siter+=cg_flt_core(par->shift[i], &local_par, F, res_flt, &out_flt[i]); /* save single precision iterations */
      /*      for (j=0; j<(par->n); ++j) {
	sflags[j]=1;
      }
      sflags[i]=1;
      siter+=cg_mshift_flt_core(sflags, &local_par, F, res_flt, out_flt); /* save single precision iterations */
#ifdef WITH_GPU
      assign_s2sd_gpu(res,&out_flt[i]);
#else
      assign_s2sd(res,&out_flt[i]);
#endif 
      spinor_field_add_assign_f(&out[i],res);
    } 
  }  

  for (i=0;i<par->n;++i){
    lprintf("INVERTER",20,"CG inversion: err2 = %1.8e < %1.8e\n",norm[i],par->err2);
  }

#ifdef WITH_GPU
  free_spinor_field_flt_gpu(res_flt);
#endif //WITH_GPU
  free_spinor_field_flt(res_flt);

#ifdef WITH_GPU
  free_spinor_field_gpu(res);
#endif //WITH_GPU
  free_spinor_field(res);

  lprintf("INVERTER",10,"CG_mshift: MVM = %d (single) - %d (double)\n",siter,diter);

  return siter+diter;
}

