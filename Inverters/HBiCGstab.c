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

/*
 * performs the BiCGstab inversion for the hermitean matrix M:
 * returns the number of cg iterations done.
 */
int HBiCGstab(MINRES_par *par, spinor_operator M, spinor_field *in, spinor_field *out){

  spinor_field *s;
  spinor_field *r, *r1, *o, *Ms, *Mo, *o0;
  spinor_field *sptmp;

  double delta, phi; 
  double alpha, beta, chi;
  double rtmp1,rtmp2,rtmp3;

  int cgiter;
  char sflags;
  unsigned short notconverged;
   
  /* allocate spinors fields and aux real variables */
  /* implementation note: to minimize the number of malloc calls
   * objects of the same type are allocated together
   */
#ifndef CHECK_SPINOR_MATCHING
    _TWO_SPINORS_MATCHING(in,out);
#endif

  s = alloc_spinor_field_f(7,in->type);
  r = s+1;
  r1 = r+1;
  o = r1+1;
  Ms = o+1;
  Mo = Ms+1;
  o0 = Mo+1;

  /* init recursion */
  cgiter = 0;
  notconverged = 1;
  beta=1.;
  spinor_field_copy_f(r, in);

  rho=z1=z2=beta;
  alpha=0.;
  spinor_field_copy_f(s, in);
  spinor_field_zero_f(out);
  sflags=1;

  chi=0;
  /* choose omega so that delta and phi are not zero */
  /* spinor_field_copy_f(o, in);  omega = in ; this may be changed */
  /* gaussian_spinor_field(o0); */
  spinor_field_copy_f(o0, in); 
  spinor_field_copy_f(o, o0);
  delta = spinor_field_prod_re_f(o, r);

  M(Ms,&s[0]);
  phi=spinor_field_prod_re_f(o0, Ms)/delta; /* o = in (see above) */


  /* cg recursion */
  do {
    ++cgiter;

    /* compute beta[0] and store in ctmp1 the old value (needed in the shifted update */
    rtmp1=beta;
    beta=-1./phi; /* b=1/phi */

    /* compute omega and chi[0] */
    spinor_field_mul_f(o,beta,Ms);
    spinor_field_add_assign_f(o,r);
    
    M(Mo, o);

    oldchi=chi;
    chi=spinor_field_prod_re_f(Mo,o)/spinor_field_sqnorm_f(Mo);
    
    /* compute r1 */
    spinor_field_mul_f(r1,chi,Mo);
    spinor_field_sub_f(r1,o,r1); 

    /* update delta and alpha[0] */
    oldalpha=alpha;
    alpha=-beta/(chi*delta);
    delta = spinor_field_prod_re_f(o0,r1); /* in = omega al passo zero */
    alpha*=delta;

    /* compute new out[0] */
    rtmp2=-beta;
    spinor_field_lc_add_assign_f(out,rtmp2,&s[0],chi,o);

    /* compute new s[0] */
    rtmp3=-alpha*chi;
    spinor_field_lc_f(Mo,alpha,&s[0],rtmp3,Ms); /* use Mo as temporary storage */
    spinor_field_add_f(&s[0],r1,Mo);

    /* assign r<-r1 */
    sptmp=r;
    r=r1;
    r1=sptmp; /*exchange pointers */

    /* update phi */
    M(Ms,&s[0]);
    phi = spinor_field_prod_re_f(o0, Ms)/delta;

    rtmp1=spinor_field_sqnorm_f(r);

    if(rtmp1<par->err2)
      notconverged=0;
    
#ifndef NDEBUG
    lprintf("INVERTER",100,"HBiCGstab iter %d res: %1.8e\n",cgiter,rtmp1);
    fflush(stdout);
#endif

  } while ((par->max_iter==0 || cgiter<par->max_iter) && notconverged);


  if (notconverged){
    lprintf("INVERTER",20,"HBiCGstab failed: err2 = %1.8e > %1.8e\n",rtmp1,par->err2);
  } else {
    lprintf("INVERTER",20,"HBiCGstab succeded: err2 = %1.8e < %1.8e\n",rtmp1,par->err2);
  }


  /* test results */
#ifndef NDEBUG
    double norm;
    M(Ms,out);
    spinor_field_mul_add_assign_f(Ms,-1.0,in);
    norm=spinor_field_sqnorm_f(Ms)/spinor_field_sqnorm_f(in);
    if (fabs(norm)>par->err2)
      lprintf("INVERTER",30,"HBiCGstab failed: err2 = %1.8e > %1.8e\n",norm,par->err2);
#endif
   
  /* free memory */
  free_spinor_field_f(s);

  /* return number of cg iter */
  return cgiter;
}





int HBiCGstab_flt(MINRES_par *par, spinor_operator_flt M, spinor_field_flt *in, spinor_field_flt *out){

  spinor_field_flt *s;
  spinor_field_flt *r, *r1, *o, *Ms, *Mo, *o0;
  spinor_field_flt *sptmp;

  float delta, phi; 
  float z1, z2, alpha, beta, chi, rho;
  float rtmp1,rtmp2,rtmp3, oldalpha,oldchi;

  int cgiter;
  char sflags;
  unsigned short notconverged;
   
  /* allocate spinors fields and aux real variables */
  /* implementation note: to minimize the number of malloc calls
   * objects of the same type are allocated together
   */
#ifndef CHECK_SPINOR_MATCHING
    _TWO_SPINORS_MATCHING(in,out);
#endif

  s = alloc_spinor_field_f_flt(7,in->type);
  r = s+1;
  r1 = r+1;
  o = r1+1;
  Ms = o+1;
  Mo = Ms+1;
  o0 = Mo+1;

  /* init recursion */
  cgiter = 0;
  notconverged = 1;
  beta=1.;
  spinor_field_copy_f_flt(r, in);

  rho=z1=z2=beta;
  alpha=0.;
  spinor_field_copy_f_flt(s, in);
  spinor_field_zero_f_flt(out);
  sflags=1;

  chi=0;
  /* choose omega so that delta and phi are not zero */
  /* spinor_field_copy_f(o, in);  omega = in ; this may be changed */
  /* gaussian_spinor_field(o0); */
  spinor_field_copy_f_flt(o0, in); 
  spinor_field_copy_f_flt(o, o0);
  delta = spinor_field_prod_re_f_flt(o, r);

  M(Ms,&s[0]);
  phi=spinor_field_prod_re_f_flt(o0, Ms)/delta; /* o = in (see above) */


  /* cg recursion */
  do {
    ++cgiter;

    /* compute beta[0] and store in ctmp1 the old value (needed in the shifted update */
    rtmp1=beta;
    beta=-1./phi; /* b=1/phi */

    /* compute omega and chi[0] */
    spinor_field_mul_f_flt(o,beta,Ms);
    spinor_field_add_assign_f_flt(o,r);
    
    M(Mo, o);

    oldchi=chi;
    chi=spinor_field_prod_re_f_flt(Mo,o)/spinor_field_sqnorm_f_flt(Mo);
    
    /* compute r1 */
    spinor_field_mul_f_flt(r1,chi,Mo);
    spinor_field_sub_f_flt(r1,o,r1); 

    /* update delta and alpha[0] */
    oldalpha=alpha;
    alpha=-beta/(chi*delta);
    delta = spinor_field_prod_re_f_flt(o0,r1); /* in = omega al passo zero */
    alpha*=delta;

    /* compute new out[0] */
    rtmp2=-beta;
    spinor_field_lc_add_assign_f_flt(out,rtmp2,&s[0],chi,o);

    /* compute new s[0] */
    rtmp3=-alpha*chi;
    spinor_field_lc_f_flt(Mo,alpha,&s[0],rtmp3,Ms); /* use Mo as temporary storage */
    spinor_field_add_f_flt(&s[0],r1,Mo);

    /* assign r<-r1 */
    sptmp=r;
    r=r1;
    r1=sptmp; /*exchange pointers */

    /* update phi */
    M(Ms,&s[0]);
    phi = spinor_field_prod_re_f_flt(o0, Ms)/delta;

    rtmp1=spinor_field_sqnorm_f_flt(r);

    if(rtmp1<(float)(par->err2))
      notconverged=0;

#ifndef NDEBUG
    lprintf("INVERTER",100,"HBiCGstab_flt iter %d res: %1.8e\n",cgiter,rtmp1);
    fflush(stdout);
#endif
  } while ((par->max_iter==0 || cgiter<par->max_iter) && notconverged);


  if (notconverged){
    lprintf("INVERTER",30,"HBiCGstab_flt failed: err2 = %1.8e > %1.8e\n",rtmp1,par->err2);
  } else {
    lprintf("INVERTER",20,"HBiCGstab_flt succeded: err2 = %1.8e < %1.8e\n",rtmp1,par->err2);
  }


  /* test results */
#ifndef NDEBUG
    double norm;
    M(Ms,out);
    spinor_field_mul_add_assign_f_flt(Ms,-1.0,in);
    norm=spinor_field_sqnorm_f_flt(Ms)/spinor_field_sqnorm_f_flt(in);
    if (fabsf(norm)>(float)(par->err2))
      lprintf("INVERTER",30,"HBiCGstab_flt failed: err2 = %1.8e > %1.8e\n",norm,par->err2);
#endif
   
  /* free memory */
  free_spinor_field_f_flt(s);

  /* return number of cg iter */
  return cgiter;
}


