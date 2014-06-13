/***************************************************************************\
 * Copyright (c) 2008, Agostino Patella, Claudio Pica                        *   
 * All rights reserved.                                                      * 
\***************************************************************************/

#include "inverters.h"
#include "linear_algebra.h"
#include "complex.h"
#include "memory.h"
#include "update.h"
#include "utils.h"
#include "logger.h"
#include "communications.h"
#include "global.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>


#undef NDEBUG


#define _print_par(c)				\
  printf("[%d] " #c " = %e\n",cgiter,c)

/* Calculate the Givens rotation (c,s) that maps:
 * |  c   s | | a | = | 1/r |
 * | -s^+ c | | b |   | 0 |
 */
#define _Givens_rot(r,c,s,a,b)						\
  { if (fabs(a)<1.e-15) { (r)=1./(b); (c)=0.; (s)=1.; }			\
    else { (s)=(b)/(a); (c)=1./sqrt(1.+(s)*(s)); (s)*=(c); (r)=1./((c)*(a)+(s)*(b)); } \
  } 


/*
 * performs the multi-shifted QMR inversion for g5-hermitean matrices:
 * out[i] = (M-(par->shift[i]))^-1 in
 * returns the number of cg iterations done.
 */
static int g5QMR_mshift_core(short *valid, mshift_par *par, spinor_operator M, spinor_field *in, spinor_field *out){
  
  spinor_field **q1,**q2;
  spinor_field *p1, *p2, *Mp;
  spinor_field *sptmp, *memall;
  
  double alpha, beta, delta, rho, innorm2; 
  double *r, *s1, *s2, *c1, *c2;
  double maxm;
  
  int i;
  int cgiter;
  unsigned int notconverged;
  
  unsigned short *flags;
  
  /* fare qualche check sugli input */
  /* par->n deve essere almeno 1! */
  assert(par->n>0);
  for(i=0;i<par->n;++i)
    _TWO_SPINORS_MATCHING(in,&out[i]);
  
  /* allocate spinors fields and aux real variables */
  /* implementation note: to minimize the number of malloc calls
   * objects of the same type are allocated together
   */
#ifdef WITH_GPU
  alloc_mem_t=GPU_MEM; /* allocate only on GPU */
#endif
  memall = alloc_spinor_field_f(2*(par->n)+3,in->type);
  alloc_mem_t=std_mem_t;
  q1 = (spinor_field**)malloc(sizeof(spinor_field*)*par->n);
  q2 = (spinor_field**)malloc(sizeof(spinor_field*)*par->n);
  for(i=0; i<par->n; i++) {
    q1[i] = memall+i;
    q2[i] = memall+par->n+i;
  }
  p1 = memall+2*par->n;
  p2 = p1+1;
  Mp = p2+1;
  //  sd = Mp+1;
  
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
  innorm2=spinor_field_sqnorm_f(in);
  if(par->n==1) { /* in this case is not a multishift and we use as starting vector out[0] */
    M.dbl(Mp,&out[0]);
    spinor_field_mul_add_assign_f(Mp,-par->shift[0],&out[0]);
    spinor_field_sub_f(p2,p2,Mp);
  }
  rho=sqrt(spinor_field_sqnorm_f(p2));
  lprintf("INVERTER",50,"g5QMR: rho init: %1.8e\n",rho*rho/innorm2);
  
  spinor_field_mul_f(p2,1./rho,p2);
  spinor_field_zero_f(p1);
  for (i=0; i<(par->n); ++i) {
    r[i]=rho;
    c2[i]=c1[i]=1.;
    s1[i]=s2[i]=0.;
    if (par->n!=1) spinor_field_zero_f(&out[i]); /* if no multishift we start with the trial solution */
    spinor_field_zero_f(q1[i]);
    spinor_field_zero_f(q2[i]);
    flags[i]=1;
  }
  delta = spinor_field_g5_prod_re_f(p2,p2);
  beta = delta;
  
  /* cg recursion */
  do {
    ++cgiter;
    
    M.dbl(Mp,p2);
    spinor_field_mul_add_assign_f(Mp,-par->shift[0],p2);
    
    /* compute alpha */
    alpha = spinor_field_g5_prod_re_f(p2,Mp)/delta;
    
    /* update p1, p2 */
    spinor_field_mul_add_assign_f(Mp,-beta,p1);
    spinor_field_mul_add_assign_f(Mp,-alpha,p2);
    sptmp=p1;
    p1=p2;
    p2=Mp;
    Mp=sptmp;
    
    /* update rho */
    rho=sqrt(spinor_field_sqnorm_f(p2));
    
    maxm=1.e-10; /* to check if the error is going down */
    
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
        
        maxm=(maxm>fabs(m))?maxm:fabs(m);
        
        /* update q */
        spinor_field_lc_f(Mp,-t*d,q1[i],-e*d,q2[i]);
        spinor_field_mul_add_assign_f(Mp,d,p1);
        sptmp=q1[i];
        q1[i]=q2[i];
        q2[i]=Mp;
        Mp=sptmp;/* swap q1[i]<-q2[i]<-Mp and Mp point to q1 */
        
        /* update solution */
        spinor_field_mul_add_assign_f(&out[i],c2[i]*r[i],q2[i]);
        
        /* update residuum */
        r[i]*=-s2[i];	
        
        if(fabs(r[i]*r[i])*(double)(1+cgiter)<par->err2*innorm2){
          flags[i]=0;
          --notconverged;
          lprintf("INVERTER",50,"g5QMR: vect %d converged at iter: %d\n",i,cgiter);
        } else {
          lprintf("INVERTER",50,"g5QMR: vect %d iter: %d res: %1.8e s2: %1.8e m=%1.8e\n",i,cgiter,fabs(r[i]*r[i])*(double)(1+cgiter)/innorm2,s2[i],m);
        }
        
      }    
    }
    
    if (fabs(rho)> 1.e-13) {
      double olddelta;
      
      /*normalize p2 */
      spinor_field_mul_f(p2,1./rho,p2);
      
      /* update delta and beta */
      olddelta=delta;
      delta = spinor_field_g5_prod_re_f(p2,p2);
      
      lprintf("INVERTER",40,"g5QMR: delta=%1.8e rho=%1.8e [iter=%d]\n",delta,rho,cgiter);
      if(fabs(delta)<1.e-6 || maxm<1.e-3) { /* the method has failed ! */
        lprintf("INVERTER",40,"g5QMR: method failed! delta=%e [iter=%d]\n",delta,cgiter);
        notconverged=0; /* exit loop */
      }
      
      /* update beta */
      beta=rho*delta/olddelta;
      
    } else { /* system has converged */
      notconverged=0;
      lprintf("INVERTER",40,"g5QMR: rho < 1.e-13 ! (system converged)\n");
    }
    
    
#ifndef NDEBUG
    if(cgiter%100==0){
      lprintf("INVERTER",40,"g5QMR: [%d] res[0]=%e notconverged=%d\n",cgiter,fabs(r[0]*r[0])*(double)(1+cgiter)/innorm2,notconverged);   
    }
#endif
    
  } while ((par->max_iter==0 || cgiter<par->max_iter) && notconverged);
  
  /* test results */
  for(i=0;i<par->n;++i){
    double norm;
    M.dbl(Mp,&out[i]);
    ++cgiter;
    if(par->shift[i]!=0.) {
      spinor_field_mul_add_assign_f(Mp,-par->shift[i],&out[i]);
    }
    spinor_field_sub_f(Mp,Mp,in);
    norm=spinor_field_sqnorm_f(Mp)/innorm2;
    valid[i]=1;
    if (fabs(norm)>par->err2){
      valid[i]=0;
      lprintf("INVERTER",30,"g5QMR failed on vect %d: err2 = %1.8e > %1.8e\n",i,norm,par->err2);
    } else {
      lprintf("INVERTER",20,"g5QMR inversion: err2 = %1.8e < %1.8e\n",norm,par->err2);
    } 
  }
  
  /* free memory */
  free_spinor_field_f(memall);
  free(q1);
  free(q2);
  free(r);
  
  free(flags);
  
  /* return number of cg iter */
  if(par->n==1) ++cgiter;
  return cgiter;
}





static int g5QMR_core_flt(short *valid, double err2, int max_iter, spinor_operator M, spinor_field_flt *in, spinor_field_flt *out){
  
  spinor_field_flt *q1,*q2;
  spinor_field_flt *p1, *p2, *Mp;
  spinor_field_flt *sptmp, *memall;
#ifndef NDEBUG
  spinor_field_flt *sdbg;
#endif
  
  double alpha, beta, delta, rho, innorm2; 
  double r, s1, s2, c1, c2;
  double maxm;
  
  int cgiter;
  unsigned int notconverged;
  
  unsigned short flag;
  
  /* fare qualche check sugli input */
  _TWO_SPINORS_MATCHING(in,out);
  
  /* allocate spinors fields and aux real variables */
  /* implementation note: to minimize the number of malloc calls
   * objects of the same type are allocated together
   */
#ifdef WITH_GPU
  alloc_mem_t=GPU_MEM; /* allocate only on GPU */
#endif
  memall = alloc_spinor_field_f_flt(2+3,in->type);
  alloc_mem_t=std_mem_t;
  q1 = memall;
  q2 = q1+1;
  p1 = q2+1;
  p2 = p1+1;
  Mp = p2+1;
  //  sd = Mp+1;

#ifndef NDEBUG
#ifdef WITH_GPU
  alloc_mem_t=GPU_MEM; /* allocate only on GPU */
#endif
  sdbg = alloc_spinor_field_f_flt(1,in->type);
  alloc_mem_t=std_mem_t;
#endif

  /* init recursion */
  cgiter = 0;
  notconverged=1;
  
  spinor_field_copy_f_flt(p2, in); /* trial solution = 0 */
  innorm2=spinor_field_sqnorm_f_flt(in);
  M.flt(Mp,out);
  spinor_field_sub_f_flt(p2,p2,Mp);
  rho=sqrt(spinor_field_sqnorm_f_flt(p2));
  lprintf("INVERTER",50,"g5QMR_core_flt: rho init: %1.8e\n",rho*rho/innorm2);

#ifndef NDEBUG
  M.flt(sdbg,out);
  spinor_field_sub_assign_f_flt(sdbg,in);
  lprintf("INVERTER",20,"g5QMR_core_flt: [%d] res=%e notconverged=%d\n",cgiter,spinor_field_sqnorm_f_flt(sdbg)/innorm2,notconverged);   
#endif
  
  spinor_field_mul_f_flt(p2,1./rho,p2);
  spinor_field_zero_f_flt(p1);
  r=rho;
  c2=c1=1.;
  s1=s2=0.;
  spinor_field_zero_f_flt(q1);
  spinor_field_zero_f_flt(q2);
  flag=1;
  delta = spinor_field_g5_prod_re_f_flt(p2,p2);
  beta = delta;
  
  /* cg recursion */
  do {
    ++cgiter;
    
    M.flt(Mp,p2);
    
    /* compute alpha */
    alpha = spinor_field_g5_prod_re_f_flt(p2,Mp)/delta;
    
    /* update p1, p2 */
    spinor_field_mul_add_assign_f_flt(Mp,-beta,p1);
    spinor_field_mul_add_assign_f_flt(Mp,-alpha,p2);
    sptmp=p1;
    p1=p2;
    p2=Mp;
    Mp=sptmp;
    
    /* update rho */
    rho=sqrt(spinor_field_sqnorm_f_flt(p2));
    
    maxm=1.e-10; /* to check if the error is going down */
    
    if (flag) {
      float a, t, e, d, m;
      a=alpha;
      t=s1*beta;
      e=c2*c1*beta+s2*a;
      m=c2*a-s2*c1*beta;
      s1=s2;
      c1=c2;
      if (fabs(m)<1.e-15) {
        d=1./rho; c2=0.; s2=1.;
      } else {
        s2=rho/m; c2=1./sqrt(1.+s2*s2); s2*=c2; d=1./(c2*m+s2*rho);
      }
      
      maxm=(maxm>fabs(m))?maxm:fabs(m);
      
      /* update q */
      spinor_field_lc_f_flt(Mp,-t*d,q1,-e*d,q2);
      spinor_field_mul_add_assign_f_flt(Mp,d,p1);
      sptmp=q1;
      q1=q2;
      q2=Mp;
      Mp=sptmp;/* swap q1<-q2<-Mp and Mp point to q1 */
      
      /* update solution */
      spinor_field_mul_add_assign_f_flt(out,c2*r,q2);
      
      /* update residuum */
      r*=-s2;
      
      if(fabs(r*r)*(float)(1+cgiter)<err2*innorm2){
        flag=0;
        --notconverged;
        lprintf("INVERTER",50,"g5QMR_core_flt: vect converged at iter: %d\n",cgiter);
      } else {
        lprintf("INVERTER",50,"g5QMR_core_flt: iter: %d res: %1.8e s2: %1.8e m=%1.8e\n",cgiter,fabs(r*r)*(float)(1+cgiter)/innorm2,s2,m);
      }
      
    }
    
    if (fabs(rho)> 1.e-5) {
      float olddelta;
      
      /*normalize p2 */
      spinor_field_mul_f_flt(p2,1./rho,p2);
      
      /* update delta and beta */
      olddelta=delta;
      delta = spinor_field_g5_prod_re_f_flt(p2,p2);
      
      lprintf("INVERTER",40,"g5QMR_core_flt: delta=%1.8e rho=%1.8e [iter=%d]\n",delta,rho,cgiter);
      if(fabs(delta)<1.e-6 || maxm<1.e-3) { /* the method has failed ! */
        lprintf("INVERTER",40,"g5QMR_core_flt: method failed! delta=%e [iter=%d]\n",delta,cgiter);
        notconverged=0; /* exit loop */
      }
      
      /* update beta */
      beta=rho*delta/olddelta;
      
    } else { /* system has converged */
      notconverged=0;
      lprintf("INVERTER",40,"g5QMR_core_flt: rho < 1.e-5 ! (system converged)\n");
    }
    
    
#ifndef NDEBUG
    if(cgiter%1==0){
      M.flt(sdbg,out);
      spinor_field_sub_assign_f_flt(sdbg,in);
      lprintf("INVERTER",20,"g5QMR_core_flt: [%d] res=%e notconverged=%d\n",cgiter,spinor_field_sqnorm_f_flt(sdbg)/innorm2,notconverged);   
    }
#endif
    
  } while ((max_iter==0 || cgiter<max_iter) && notconverged);
  
  /* test results */
  float norm;
  M.flt(Mp,out);
  ++cgiter;
  spinor_field_sub_f_flt(Mp,Mp,in);
  norm=spinor_field_sqnorm_f_flt(Mp)/innorm2;
  *valid=1;
  if (fabs(norm)>err2){
    *valid=0;
    lprintf("INVERTER",30,"g5QMR_core_flt failed: err2 = %1.8e > %1.8e\n",norm,err2);
  } else {
    lprintf("INVERTER",20,"g5QMR_core_flt inversion: err2 = %1.8e < %1.8e\n",norm,err2);
  } 
  
  /* free memory */
  free_spinor_field_f_flt(memall);

#ifndef NDEBUG
  free_spinor_field_f_flt(sdbg);
#endif
  
  /* return number of cg iter */
  ++cgiter;
  return cgiter;
}




static int g5QMR_core(short *valid, double err2, int max_iter, spinor_operator M, spinor_field *in, spinor_field *out){
  
  spinor_field *q1,*q2;
  spinor_field *p1, *p2, *Mp;
  spinor_field *sptmp, *memall;
#ifndef NDEBUG
  spinor_field *sdbg;
#endif
  
  double alpha, beta, delta, rho, innorm2; 
  double r, s1, s2, c1, c2;
  double maxm;
  
  int cgiter;
  unsigned int notconverged;
  
  unsigned short flag;
  
  /* fare qualche check sugli input */
  _TWO_SPINORS_MATCHING(in,out);
  
  /* allocate spinors fields and aux real variables */
  /* implementation note: to minimize the number of malloc calls
   * objects of the same type are allocated together
   */
#ifdef WITH_GPU
  alloc_mem_t=GPU_MEM; /* allocate only on GPU */
#endif
  memall = alloc_spinor_field_f(2+3,in->type);
  alloc_mem_t=std_mem_t;
  q1 = memall;
  q2 = q1+1;
  p1 = q2+1;
  p2 = p1+1;
  Mp = p2+1;
  //  sd = Mp+1;

#ifndef NDEBUG
#ifdef WITH_GPU
  alloc_mem_t=GPU_MEM; /* allocate only on GPU */
#endif  
  sdbg = alloc_spinor_field_f(1,in->type);
  alloc_mem_t=std_mem_t;
#endif

  /* init recursion */
  cgiter = 0;
  notconverged=1;
  
  spinor_field_copy_f(p2, in); /* trial solution = 0 */
  innorm2=spinor_field_sqnorm_f(in);
  M.dbl(Mp,out);
  spinor_field_sub_f(p2,p2,Mp);
  rho=sqrt(spinor_field_sqnorm_f(p2));
  lprintf("INVERTER",50,"g5QMR_core: rho init: %1.8e\n",rho*rho/innorm2);

#ifndef NDEBUG
  M.dbl(sdbg,out);
  spinor_field_sub_assign_f(sdbg,in);
  lprintf("INVERTER",20,"g5QMR_core: [%d] res=%e notconverged=%d\n",cgiter,spinor_field_sqnorm_f(sdbg)/innorm2,notconverged);   
#endif
  
  spinor_field_mul_f(p2,1./rho,p2);
  spinor_field_zero_f(p1);
  r=rho;
  c2=c1=1.;
  s1=s2=0.;
  spinor_field_zero_f(q1);
  spinor_field_zero_f(q2);
  flag=1;
  delta = spinor_field_g5_prod_re_f(p2,p2);
  beta = delta;
  
  /* cg recursion */
  do {
    ++cgiter;
    
    M.dbl(Mp,p2);
    
    /* compute alpha */
    alpha = spinor_field_g5_prod_re_f(p2,Mp)/delta;
    
    /* update p1, p2 */
    spinor_field_mul_add_assign_f(Mp,-beta,p1);
    spinor_field_mul_add_assign_f(Mp,-alpha,p2);
    sptmp=p1;
    p1=p2;
    p2=Mp;
    Mp=sptmp;
    
    /* update rho */
    rho=sqrt(spinor_field_sqnorm_f(p2));
    
    maxm=1.e-10; /* to check if the error is going down */
    
    if (flag) {
      double a, t, e, d, m;
      a=alpha;
      t=s1*beta;
      e=c2*c1*beta+s2*a;
      m=c2*a-s2*c1*beta;
      s1=s2;
      c1=c2;
      if (fabs(m)<1.e-15) {
        d=1./rho; c2=0.; s2=1.;
      } else {
        s2=rho/m; c2=1./sqrt(1.+s2*s2); s2*=c2; d=1./(c2*m+s2*rho);
      }
      
      maxm=(maxm>fabs(m))?maxm:fabs(m);
      
      /* update q */
      spinor_field_lc_f(Mp,-t*d,q1,-e*d,q2);
      spinor_field_mul_add_assign_f(Mp,d,p1);
      sptmp=q1;
      q1=q2;
      q2=Mp;
      Mp=sptmp;/* swap q1<-q2<-Mp and Mp point to q1 */
      
      /* update solution */
      spinor_field_mul_add_assign_f(out,c2*r,q2);
      
      /* update residuum */
      r*=-s2;
      
      if(fabs(r*r)*(double)(1+cgiter)<err2*innorm2){
        flag=0;
        --notconverged;
        lprintf("INVERTER",50,"g5QMR_core: vect converged at iter: %d\n",cgiter);
      } else {
        lprintf("INVERTER",50,"g5QMR_core: iter: %d res: %1.8e s2: %1.8e m=%1.8e\n",cgiter,fabs(r*r)*(double)(1+cgiter)/innorm2,s2,m);
      }
      
    }
    
    if (fabs(rho)> 1.e-5) {
      double olddelta;
      
      /*normalize p2 */
      spinor_field_mul_f(p2,1./rho,p2);
      
      /* update delta and beta */
      olddelta=delta;
      delta = spinor_field_g5_prod_re_f(p2,p2);
      
      lprintf("INVERTER",40,"g5QMR_core: delta=%1.8e rho=%1.8e [iter=%d]\n",delta,rho,cgiter);
      if(fabs(delta)<1.e-6 || maxm<1.e-3) { /* the method has failed ! */
        lprintf("INVERTER",40,"g5QMR_core: method failed! delta=%e [iter=%d]\n",delta,cgiter);
        notconverged=0; /* exit loop */
      }
      
      /* update beta */
      beta=rho*delta/olddelta;
      
    } else { /* system has converged */
      notconverged=0;
      lprintf("INVERTER",40,"g5QMR_core: rho < 1.e-5 ! (system converged)\n");
    }
    
#ifndef NDEBUG
    if(cgiter%1==0){
      M.dbl(sdbg,out);
      spinor_field_sub_assign_f(sdbg,in);
      lprintf("INVERTER",20,"g5QMR_core: [%d] res=%e notconverged=%d\n",cgiter,spinor_field_sqnorm_f(sdbg)/innorm2,notconverged);   
    }
#endif
    
  } while ((max_iter==0 || cgiter<max_iter) && notconverged);
  
  /* test results */
  double norm;
  M.dbl(Mp,out);
  ++cgiter;
  spinor_field_sub_f(Mp,Mp,in);
  norm=spinor_field_sqnorm_f(Mp)/innorm2;
  *valid=1;
  if (fabs(norm)>err2){
    *valid=0;
    lprintf("INVERTER",30,"g5QMR_core failed: err2 = %1.8e > %1.8e\n",norm,err2);
  } else {
    lprintf("INVERTER",20,"g5QMR_core inversion: err2 = %1.8e < %1.8e\n",norm,err2);
  } 
  
  /* free memory */
  free_spinor_field_f(memall);

#ifndef NDEBUG
  free_spinor_field_f(sdbg);
#endif

  
  /* return number of cg iter */
  ++cgiter;
  return cgiter;
}

 


static double sh;
static spinor_operator g5Herm;
static void Herm_dbl(spinor_field *out, spinor_field *in){
  g5Herm.dbl(out,in);
  if(sh!=0.) {
    spinor_field_mul_add_assign_f(out,-sh,in);
  }
  spinor_field_g5_f(out,out);
}

static void Herm_flt(spinor_field_flt *out, spinor_field_flt *in){
  g5Herm.flt(out,in);
  if(sh!=0.) {
    spinor_field_mul_add_assign_f_flt(out,(float)(-sh),in);
  }
  spinor_field_g5_f_flt(out,out);
}

static spinor_operator Herm={&Herm_dbl,&Herm_flt};

int g5QMR_mshift(mshift_par *par, spinor_operator M, spinor_field *in, spinor_field *out){
  int cgiter;
  int n;
  mshift_par orig;
  short *valid;
  int loccg;
  int msiter;

  orig=*par; /* save par */
  valid=malloc(sizeof(short)*orig.n);
  cgiter=g5QMR_mshift_core(valid,par,M,in,out);
  msiter=cgiter;

  /* if some vector has not converged try non-multishift */
  par->n=1;
  for(n=0;n<orig.n;++n){
    if (valid[n]==0) {
#ifndef NDEBUG
      lprintf("INVERTER",20,"Doing non multishift on vector %d\n",n);
#endif
      par->shift=orig.shift+n;
      par->max_iter=0;
      loccg=g5QMR_mshift_core(valid+n,par,M,in,out+n);
      cgiter+=loccg;

      if (valid[n]==0) {
	/* revert to MINRES */
	MINRES_par Mpar;
				
	lprintf("INVERTER",20,"g5QMR_mshift: using MINRES on vect %d\n",n);

	Mpar.err2=par->err2;
	Mpar.max_iter=0;

	g5Herm=M;
	sh=par->shift[0];

	spinor_field_g5_f(in,in); /* multiply input by g5 for MINRES */
	loccg=MINRES(&Mpar,Herm,in,&out[n],(n==0)?&out[n]:&out[n-1]);
	spinor_field_g5_f(in,in); /* restore input vector */

	cgiter+=loccg;

      }

    }
  }

  /* this is for debug purposes */
	
  *par=orig; /* restore par */
  
  free(valid);
  
  lprintf("INVERTER",10,"g5QMR_mshift: MVM = %d/%d\n",msiter,cgiter);
  
  return cgiter;
  
}


int g5QMR_fltacc( g5QMR_fltacc_par* par, spinor_operator M, spinor_field *in, spinor_field *out)
{
  int cgiter=0,cgiter_flt=0,cgiter_minres=0;
  short valid;
  spinor_field_flt *in_flt, *out_flt, *res_flt;
  spinor_field *res, *Mp;
  double err2, innorm2;

#ifdef WITH_GPU
  alloc_mem_t=GPU_MEM; /* allocate only on GPU */
#endif
  res = alloc_spinor_field_f(2,in->type);
  Mp = res+1;
  in_flt = alloc_spinor_field_f_flt(3,in->type);
  alloc_mem_t=std_mem_t; /* set the allocation memory type back */
  out_flt = in_flt+1;
  res_flt = out_flt+1;

  innorm2 = spinor_field_sqnorm_f(in);

  /* 0. res = in - Mp*out */
  spinor_field_copy_f(res,in);
  M.dbl(Mp, out);
  spinor_field_sub_assign_f(res,Mp);

  err2 = spinor_field_sqnorm_f(res);
  lprintf("INVERTER",50,"g5QMR_fltacc: res=%e\n", err2/innorm2);

   while(err2/innorm2>par->err2) {
    /* 1. M^{-1} res = out + M^{-1} in */
    /* 2. out + M^{-1} res -> out */
    /* 3. res = in - M.out */

    assign_sd2s(res_flt,res);
    spinor_field_zero_f_flt(out_flt);
    cgiter_flt += g5QMR_core_flt(&valid,par->err2_flt,par->max_iter_flt,M,res_flt,out_flt);
    if (valid==0) {
#if 0
      lprintf("INVERTER",20,"g5QMR_fltacc: using GMRES\n");
      inverter_par Mpar;
      Mpar.err2=par->err2_flt;
      Mpar.max_iter=0;
      Mpar.kry_dim=10;
      g5Herm=M;
      sh=0.;
      //      spinor_field_g5_f_flt(res_flt,res_flt); /* multiply input by g5 for MINRES */
      cgiter_minres+=GMRES_flt(&Mpar,M,res_flt,out_flt,out_flt);
      //      spinor_field_g5_f_flt(res_flt,res_flt); /* restore input vector */
#endif
      lprintf("INVERTER",20,"g5QMR_fltacc: using MINRES\n");
      MINRES_par Mpar;
      Mpar.err2=par->err2_flt;
      Mpar.max_iter=0;
      g5Herm=M;
      sh=0.;
      spinor_field_g5_f_flt(res_flt,res_flt); /* multiply input by g5 for MINRES */
      cgiter_minres+=MINRES_flt(&Mpar,Herm,res_flt,out_flt,out_flt);
      spinor_field_g5_f_flt(res_flt,res_flt); /* restore input vector */
    }

    assign_s2sd(res,out_flt);
    spinor_field_add_assign_f(out, res);

    M.dbl(res,out);
    spinor_field_sub_f(res,in,res);

    err2 = spinor_field_sqnorm_f(res);
    cgiter++;

    lprintf("INVERTER",50,"g5QMR_fltacc: res=%e\n", err2/innorm2);
  }

  lprintf("INVERTER",10,"g5QMR_fltacc: MVM (g5QMR_flt,g5QMR,MINRES) = %d ; %d ; %d\n",cgiter_flt,cgiter,cgiter_minres);

  free_spinor_field_f(res);
  free_spinor_field_f_flt(in_flt);

  return cgiter+cgiter_flt+cgiter_minres;
}


int g5QMR_flt(inverter_par *par, spinor_operator M, spinor_field_flt *in, spinor_field_flt *out){
  int iter, rep=0;
  short int valid=0;

  iter = g5QMR_core_flt(&valid,par->err2_flt,par->max_iter_flt,M,in,out);
  while (!valid && (par->max_iter==0 || iter < par->max_iter)){
    iter += g5QMR_core_flt(&valid,par->err2_flt,par->max_iter_flt,M,in,out);
    if ((++rep %10 == 0 )){
      lprintf("INVERTER",-10,"g5QMR recursion = %d (precision too high?)\n",rep);
    }
  }
  lprintf("INVERTER",10,"g5QMR_flt: MVM_flt = %d\n",iter);
  return iter;
}
