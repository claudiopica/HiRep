/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *   
* All rights reserved.                                                      * 
\***************************************************************************/

#include "global.h"
#include "inverters.h"
#include "linear_algebra.h"
#include "rational_functions.h"
#include "error.h"
#include "memory.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "logger.h"

/* include rational approximations database */
#include "approx_data.db"

/* Functions for manipulation of rational_app structures */

/* r_app_alloc:
 * read the specified app->order and allocate enough
 * space for the coefficients
 */
void r_app_alloc(rational_app *app) {
  app->_asize=app->order;
  app->a=malloc(sizeof(*(app->a))*(app->order*2+1));
  app->b=app->a+(app->order+1);
}

/* r_app_free:
 * free memory space used by the approximation
 */
void r_app_free(rational_app *app) {
  app->_asize=0;
  free(app->a);
  app->a=0;
  app->b=0;
}

/* r_app_realloc:
 * read the new app->order and make enough space for 
 * the coefficients
 */
void r_app_realloc(rational_app *app){
  if (app->order>app->_asize) {
    r_app_free(app);
    r_app_alloc(app);
  }
}

/* Functions for coefficients computation */

/* converts the coef from the root/poles form to the partial fraction exp
 * before : r(x)=a[0]*(x-a[1])/(x-b[0])*...*(x-a[n+1])/(x-b[n])
 * after  : r(x)=a[0]+a[1]/(x-b[0])+a[2]/(x-b[1])+...+a[n+1]/(x-b[n])
 */
static void r_app_rp2pfe(rational_app *app) {
  unsigned int i,j;
  double *ctmp;

  ctmp=malloc(sizeof(*ctmp)*(app->order));
  for(i=0;i<app->order;++i)
    ctmp[i]=app->a[i+1];

  /* only a[1]...a[n+1] changes */
  for(i=0;i<app->order;++i) {
    app->a[i+1]=app->a[0];
    for(j=0;j<app->order;++j) {
      app->a[i+1]*=(i==j)?(app->b[i]-ctmp[j]):((app->b[i]-ctmp[j])/(app->b[i]-app->b[j]));
    }
  }

  free(ctmp);
  
}

/* converts the coef from the root/poles form to the partial fraction exp
 * as r_app_rp2pfe but for the inverse function
 */
static void r_app_rp2pfe_inv(rational_app *app) {
  unsigned int i;
  double ctmp;
  
  app->a[0]=1./app->a[0];
  /* swap poles and roots: a[i+1] <-> b[i] */
  for(i=0;i<app->order;++i) {
    ctmp=app->a[i+1];
    app->a[i+1]=app->b[i];
    app->b[i]=ctmp;
  }
  
  r_app_rp2pfe(app);

}

/* rescale the approximation of a factor k 
 * used by r_app_set to adjust approximation intervals
 */
void r_app_rescale(rational_app *app, double k) {
  unsigned int i;
  double fexp=-fabs(((double)app->n)/((double)app->d));
  app->a[0]*=pow(k,fexp);
  for(i=0;i<app->order;){
    app->b[i]*=k;
    ++i;
    app->a[i]*=k;
  }
}

/* DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG */
/* THE FOLLOWING ARE DEBUG FUNCTIONS TO TEST RATIONAL APPROXIMATIONS */
/*
  static void print_delta(rational_app *app, double min, double max) {
  const int n=1000;
  double x=min;

  while(x<max){
  int i;
  double check=pow(x,-fabs(((double)app->n)/((double)app->d)));
  double fx=app->a[0];
  for (i=0;i<app->order;++i)
  fx*=(x-app->a[i+1])/(x-app->b[i]);
  check=fabs(fx-check)/check;
  if (app->rel_error<check)
  printf("x=%1.5e relerr=%1.5e\n",x,check);
  x+=(max-min)/(double)n;
  }
  }

  static void print_delta2(rational_app *app, double min, double max) {
  const int n=1000;
  double x=min;

  while(x<max){
  int i;
  double check=pow(x,((double)app->n)/((double)app->d));
  double fx=app->a[0];
  for (i=0;i<app->order;++i)
  fx+=app->a[i+1]/(x-app->b[i]);
  check=fabs(fx-check)/check;
  if (app->rel_error<check)
  printf("x=%1.5e relerr=%1.5e\n",x,check);
  x+=(max-min)/(double)n;
  }
  }
*/
/* DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG */
/* END END END END END END END END END END END END END END END END */ 

/* fill up rational_app structure with a r_app
 * of degree app->order which is an approximant for the function
 * f(x)=x^(app->n/app->d) in the interval [min,max]
 * with relative error app->rel_error
 */
void r_app_set(rational_app *app, double min, double max) {
  /* scan approc_data database for the best approximation */
  int bn,bd,bo;
  double berr,bmin,bmax;
  double req_e=min/max;
  const double *best=0;
  const double *cur_app=&approx_data[0];
  while(*cur_app!=0.) {
    bn=(int) cur_app[0];
    bd=(int) cur_app[1];
    bo=(int) cur_app[2];
    berr= cur_app[3];
    bmin=cur_app[4];
    bmax=cur_app[5];
    if (bn==abs(app->n) &&
        bd==abs(app->d) &&
        berr<(app->rel_error) &&
        (bmin/bmax)<req_e
        ){
      /* we have found a candidate approximation */
      /* the selection rules are:
       * 1) the approximation function must match ;
       * 2) the relative error must be less than the required rel_error;
       * 3) the approx interval must contain the required one.
       */
      if (best==0 || ((int)best[2])>bo) {
        /* this is a better approx: either the first approx found or
         * one with greater error (i.e. smaller order) but still
         * the required precision
         */
        best=cur_app;
      }
    }
    /* go to next approx */
    cur_app += 2*bo+7;
  }
  /* if we cannot find a rational approx give an error */
  error((best==0),1,"r_app_set",
        "Failed to find a suitable rational approximation.\nPlease adjust database and try again.\n");

  /* now set up coefficients in app */
  bo=(int) best[2];
  berr= best[3];
  bmin= best[4];
  bmax= best[5];
  app->order=bo;

  lprintf("RAPPROX",0,"Best approx found: o=%d err=%1.5e bmin=%1.5e bmax=%1.5e\n",bo,berr,bmin,bmax);

  r_app_realloc(app);
  app->a[0]=best[6];
  for(bn=0;bn<bo;) {
    app->b[bn]=best[7+bo+bn];
    ++bn;
    app->a[bn]=best[6+bn];
  }
  /* rescale to requested interval */
  r_app_rescale(app,max/bmax);
  app->min=min;
  app->max=max;
	
  /* change to pfe representation */
  if((app->n*app->d)<0) { /* check sign for inverse function */
    /* this is indeed correct: the database contains coef for f(x)=x^-k */
    r_app_rp2pfe(app);
  } else {
    r_app_rp2pfe_inv(app);
  }

}

/*
 * computes: out = (a[0]*I + \sum_i a[i]*(Q-b[i])^-1 ) in
 * this MUST work in the case: out==in 
 * where Q is a linear operator acting on spinors
 * This implementation uses CG_mshift => Q must be hermitean and positive definite!
 */
void rational_func(rational_app *coef, spinor_operator Q, spinor_field *out, spinor_field *in) {
  static mshift_par inv_par;
  spinor_field *inv_out;
  unsigned int i;

  /* check input types */
#ifndef CHECK_SPINOR_MATCHING
  _TWO_SPINORS_MATCHING(in,out);
#endif

  /* allocate cg output vectors */
  inv_out = alloc_spinor_field_f(coef->order,in->type);

  /* set up cg parameters */
  inv_par.n = coef->order;
  inv_par.shift = coef->b;
  inv_par.err2=coef->rel_error/coef->order;    /* CAMBIARE: METTERE PARAMETRI COMUNI ALL'UPDATE */
  inv_par.err2*=inv_par.err2;
#define EPSILON 1.e-25
  if(inv_par.err2<EPSILON) inv_par.err2=EPSILON;
#undef EPSILON
  inv_par.max_iter=0; /* no limit */
   
  /* compute inverse vectors */
  if (inv_par.n==1) { spinor_field_zero_f(inv_out); }
  cg_mshift(&inv_par, Q, in, inv_out);

  /* sum all the contributions */
  spinor_field_mul_f(out,coef->a[0],in);
  for (i=1; i<(coef->order)+1; ++i) {
    spinor_field_mul_add_assign_f(out,coef->a[i],&inv_out[i-1]);
  }

  /* free memory */
  free_spinor_field_f(inv_out);

}
