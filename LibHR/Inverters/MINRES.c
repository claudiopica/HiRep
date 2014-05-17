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
#include "communications.h"
#include <stdlib.h>
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
static int MINRES_core(short int *valid, MINRES_par *par, spinor_operator M, spinor_field *in, spinor_field *out, spinor_field *trial){

  spinor_field *q1,*q2;
  spinor_field *p1, *p2, *Mp;
  spinor_field *sptmp, *memall;

  double alpha, beta, oldbeta, innorm2; 
  double r, s1, s2, c1, c2, rho1, rho2, rp;
  double d, h, k;

  int cgiter;
  unsigned short notconverged;

  /* fare qualche check sugli input */
  _TWO_SPINORS_MATCHING(in,out);
  if(trial!=NULL) {_TWO_SPINORS_MATCHING(in,trial);}

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

  memall = alloc_spinor_field_f(5,in->type);
  q1=memall;
  q2= q1+1;
  p1 = q2+1;
  p2 = p1+1;
  Mp = p2+1;

  /* init recursion */
  cgiter = 0;
  notconverged=1;

  spinor_field_copy_f(p2, in);
  if(trial!=NULL) {
    M(p1,trial);
    ++cgiter;
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
    }
		/* just for debug 
		else {
      lprintf("INVERTER",30,"MINRES iter %d res: %1.8e\n",cgiter,(r*r)/innorm2);
    }
		*/
	
	
  } while ((par->max_iter==0 || cgiter<par->max_iter) && notconverged);

  /* test results */
  M(Mp,out);
  ++cgiter;
  spinor_field_sub_f(Mp,Mp,in);
  innorm2=spinor_field_sqnorm_f(Mp)/innorm2;
  *valid=1;
  if (fabs(innorm2)>par->err2) {
    *valid=0;
    lprintf("INVERTER",30,"MINRES failed: err2 = %1.8e > %1.8e\n",innorm2,par->err2);
  } else {
    lprintf("INVERTER",20,"MINRES inversion: err2 = %1.8e < %1.8e\n",innorm2,par->err2);
  } 
   
  /* free memory */
  free_spinor_field_f(memall);

  /* return number of cg iter */
  return cgiter;
}

int MINRES(MINRES_par *par, spinor_operator M, spinor_field *in, spinor_field *out, spinor_field *trial){
  int iter,rep=0;
  short int valid;

  iter=MINRES_core(&valid, par, M, in, out, trial);
  while(!valid && (par->max_iter==0 || iter<par->max_iter)) {
    iter+=MINRES_core(&valid, par, M, in, out, out);
    if((++rep)%5==0)
      lprintf("INVERTER",-10,"MINRES recursion = %d (precision too high?)\n",rep);
  }

  lprintf("INVERTER",10,"MINRES: MVM = %d\n",iter);

  return iter;
}



static double spinor_field_prod_re_f_f2d(spinor_field_flt *s1, spinor_field_flt *s2)
{
  double res=0.;
  
  _TWO_SPINORS_FOR_SUM(s1,s2,res) {
    float prod;
    _spinor_prod_re_f(prod,*_SPINOR_PTR(s1),*_SPINOR_PTR(s2));
    res+=(double)prod;
  }
#ifdef WITH_MPI
  global_sum(&res,1);
#endif
  return res;
}

static double spinor_field_sqnorm_f_f2d(spinor_field_flt *s1)
{
  double res=0.; 
   
  _ONE_SPINOR_FOR_SUM(s1,res) {
    float prod;
    _spinor_prod_re_f(prod,*_SPINOR_PTR(s1),*_SPINOR_PTR(s1));
    res+=(double)prod;
  }
#ifdef WITH_MPI
  global_sum(&res,1);
#endif
  return res;
}


static int MINRES_core_flt(short int *valid, MINRES_par *par, spinor_operator_flt M, spinor_field_flt *in, spinor_field_flt *out, spinor_field_flt *trial){

  spinor_field_flt *q1,*q2;
  spinor_field_flt *p1, *p2, *Mp;
  spinor_field_flt *sptmp, *memall;

  double alpha, beta, oldbeta, innorm2; 
  double r, s1, s2, c1, c2, rho1, rho2, rp;
  double d, h, k;

  int cgiter;
  unsigned short notconverged;

  /* fare qualche check sugli input */
  _TWO_SPINORS_MATCHING(in,out);
  if(trial!=NULL) {_TWO_SPINORS_MATCHING(in,trial);}
     
  /* allocate spinors fields and aux real variables */
  /* implementation note: to minimize the number of malloc calls
   * objects of the same type are allocated together
   */

  memall = alloc_spinor_field_f_flt(5,in->type);
  q1=memall;
  q2= q1+1;
  p1 = q2+1;
  p2 = p1+1;
  Mp = p2+1;

  /* init recursion */
  cgiter = 0;
  notconverged=1;

  spinor_field_copy_f_flt(p2, in);
  if(trial!=NULL) {
    M(p1,trial);
    ++cgiter;
    spinor_field_sub_assign_f_flt(p2,p1);
    if(out!=trial){
      spinor_field_copy_f_flt(out,trial);
    }
    
  } else {
    spinor_field_zero_f_flt(out);
  }

  innorm2=spinor_field_sqnorm_f_f2d(in);
  beta=sqrt(spinor_field_sqnorm_f_f2d(p2));
  spinor_field_mul_f_flt(p2,(float)(1./beta),p2);
  spinor_field_zero_f_flt(p1);  
  r=rho2=beta;
  rho1=1.;
  c2=-1.;
  rp=s1=s2=c1=0.;
  spinor_field_zero_f_flt(q1);
  spinor_field_zero_f_flt(q2);

  /* cg recursion */
  do {
    ++cgiter;

    M(Mp,p2);
    
    /* compute alpha */
    alpha = spinor_field_prod_re_f_f2d(Mp,p2);

    /* update p1, p2 */
    spinor_field_mul_add_assign_f_flt(Mp,-(float)beta,p1);
    spinor_field_mul_add_assign_f_flt(Mp,-(float)alpha,p2);
    sptmp=p1;
    p1=p2;
    p2=Mp;
    Mp=sptmp;

    /* update beta */
    oldbeta=beta;
    beta=sqrt(spinor_field_sqnorm_f_f2d(p2));
    
    /* normalize p2 */
    spinor_field_mul_f_flt(p2,(float)(1./beta),p2);


    d=(alpha-rp*c1)*s2;
    h=oldbeta*s1;
    rp=-oldbeta*c1*s2-alpha*c2;
    k=sqrt(rp*rp+beta*beta);
    c1=c2;
    c2=rp/k;
    s1=s2;
    s2=beta/k;

    spinor_field_lc_f_flt(Mp,(float)(-h/rho1),q1,(float)(-d/rho2),q2);
    sptmp=q1;
    q1=q2;
    q2=sptmp; /* swap q1[i]<->q2[i] */
    spinor_field_add_f_flt(q2,p1,Mp);
	
    /* update rho */
    rho1=rho2;
    rho2=k;
	
    /* update solution */
    spinor_field_mul_add_assign_f_flt(out,(float)(r*c2/k),q2);
	
    /* update residuum */
    r*=s2;

    if((r*r)<par->err2*innorm2){
      notconverged=0;
    }
		/* just for debug 
		else {
      lprintf("INVERTER",30,"MINRES iter %d res: %1.8e\n",cgiter,(r*r)/innorm2);
    }
		*/
	
	
  } while ((par->max_iter==0 || cgiter<par->max_iter) && notconverged);

  /* test results */
  M(Mp,out);
  ++cgiter;
  spinor_field_sub_f_flt(Mp,Mp,in);
  innorm2=spinor_field_sqnorm_f_f2d(Mp)/innorm2;
  *valid=1;
  if (fabs(innorm2)>par->err2) {
    *valid=0;
    lprintf("INVERTER",30,"MINRES failed: err2 = %1.8e > %1.8e\n",innorm2,par->err2);
  } else {
    lprintf("INVERTER",20,"MINRES inversion: err2 = %1.8e < %1.8e\n",innorm2,par->err2);
  } 
   
  /* free memory */
  free_spinor_field_f_flt(memall);

  /* return number of cg iter */
  return cgiter;
}

int MINRES_flt(MINRES_par *par, spinor_operator_flt M, spinor_field_flt *in, spinor_field_flt *out, spinor_field_flt *trial){
  int iter,rep=0;
  short int valid;

  iter=MINRES_core_flt(&valid, par, M, in, out, trial);
  while(!valid && (par->max_iter==0 || iter<par->max_iter)) {
    iter+=MINRES_core_flt(&valid, par, M, in, out, out);
    if((++rep)%5==0)
      lprintf("INVERTER",-10,"MINRES recursion = %d (precision too high?)\n",rep);
  }

  lprintf("INVERTER",10,"MINRES: MVM = %d\n",iter);

  return iter;
}


