/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *   
* All rights reserved.                                                      * 
\***************************************************************************/

#include "inverters.h"
#include "linear_algebra.h"
#include "complex.h"
#include "update.h"
#include "memory.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

/* _compute_z(z3+i,z1+i,z2+i,&ctmp1,beta,alpha,&(par->shift[i-1])); */
/* res = (z1*betam1)/(beta*alpha*(z1-z2)+z1*betam1*(1+sigma*beta)) (abbiamo diviso per z2) */
__inline static void _compute_z(complex *res, complex *z1, complex *z2, 
				complex *betam1, complex *beta, complex *alpha,
				double *sigma) 
{
  complex ctmp1, ctmp2, ctmp3, ctmp4;
  
  _complex_mul(ctmp1,*z1,*betam1);
  _complex_mul(ctmp2,*beta,*alpha); 
  _complex_sub(*res,*z1,*z2);
  _complex_mulr(ctmp4,*sigma,*beta);
  _complex_mul(ctmp3,ctmp2,*res);
  ctmp2=ctmp1;
  _complex_add_1(ctmp4);
  _complex_mul_assign(ctmp3,ctmp1,ctmp4);
  /*  _complex_mul(ctmp2,ctmp1,*z2); */
  _complex_div(*res,ctmp2,ctmp3);

}

#define _print_complex(c) \
  printf("[%d] " #c " = (%e,%e)\n",cgiter,c.re,c.im)

/*
 * performs the multi-shifted CG inversion:
 * out[i] = (M-(par->shift[i]))^-1 in
 * returns the number of cg iterations done.
 */
int BiCGstab_mshift(mshift_par *par, spinor_operator M, spinor_field *in, spinor_field *out){

  spinor_field *s;
  spinor_field *r, *r1, *o, *Ms, *Mo, *o0;
  spinor_field *sptmp;

  complex delta, phi; 
  complex *z1, *z2, *z3, *alpha, *beta, *chi, *rho; /* alpha is unnecessary */
  complex ctmp1, ctmp2, ctmp3,ctmp4, ctmp5, oldalpha;
  double rtmp1;

  int i;
  int cgiter;
  char *sflags;
  unsigned short notconverged;
	unsigned int spinorlen;
   
  /* fare qualche check sugli input */
  /* par->n deve essere almeno 2! */
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
  s = alloc_spinor_field_f((par->n)+6);
  r = s+par->n;
  r1 = r+1;
  o = r1+1;
  Ms = o+1;
  Mo = Ms+1;
  o0 = Mo+1;

  z1 = (complex *)malloc(sizeof(complex)*7*(par->n));
  z2 = z1+(par->n);
  z3 = z2+(par->n);
  alpha = z3+(par->n);
  beta = alpha+(par->n);
  chi = beta+(par->n);
  rho = chi+(par->n);

  sflags = (char *)malloc(sizeof(char)*(par->n));
   
  /* init recursion */
  cgiter = 0;
  notconverged = 1;
  _complex_1(beta[0]);
  spinor_field_copy_f(r, in);
  for (i=0; i<(par->n); ++i) {
    rho[i]=z1[i]=z2[i]=beta[0];
    _complex_0(alpha[i]);
    spinor_field_copy_f(&s[i], in);
    spinor_field_zero_f(&out[i]);
    sflags[i]=1;
  }
  /* choose omega so that delta and phi are not zero */
  /* spinor_field_copy_f(o, in);  omega = in ; this may be changed */
  gaussian_spinor_field(o0); 
  /* spinor_field_copy_f(o0, in); */
  spinor_field_copy_f(o, o0);
  delta = spinor_field_prod_f(o, r);

  M(Ms,&s[0]);
  ctmp1 = spinor_field_prod_f(o0, Ms); /* o = in (see above) */
  _complex_div(phi,ctmp1,delta);

  _print_complex(delta);
  _print_complex(phi);
  _print_complex(alpha[0]);
  _print_complex(beta[0]);


  /* cg recursion */
  do {
    ++cgiter;

    /* compute beta[0] and store in ctmp1 the old value (needed in the shifted update */
    ctmp1=beta[0];
    _complex_inv(beta[0], phi); /* b=1/phi */
    _complex_minus(beta[0],beta[0]); /* b=-b */

    /* compute omega and chi[0] */
    spinor_field_mulc_f(o,beta[0],Ms);
    spinor_field_add_assign_f(o,r);
    
    M(Mo,o);

    ctmp2=spinor_field_prod_f(Mo,o);
    rtmp1=1./spinor_field_sqnorm_f(Mo);
    _complex_mulr(chi[0],rtmp1,ctmp2);
    
    /* compute r1 */
    spinor_field_mulc_f(r1,chi[0],Mo);
    spinor_field_sub_f(r1,o,r1); 

    /* update delta and alpha[0] */
    oldalpha=alpha[0];
    _complex_mul(ctmp2,delta,chi[0]);
    delta = spinor_field_prod_f(o0,r1); /* in = omega al passo zero */
    _complex_minus(ctmp2,ctmp2);
    _complex_mul(ctmp3,beta[0],delta);
    _complex_div(alpha[0],ctmp3,ctmp2); /* alpha[0]=-beta[0]*delta/(deltam1*chi[0]) */

    /* compute new out[0] */
    _complex_minus(ctmp2,beta[0]);
    spinor_field_clc_add_assign_f(&out[0],ctmp2,&s[0],chi[0],o);

    /* compute new s[0] */
    _complex_mul(ctmp2,alpha[0],chi[0]);
    _complex_minus(ctmp2,ctmp2);
    spinor_field_clc_f(Mo,alpha[0],&s[0],ctmp2,Ms); /* use Mo as temporary storage */
    spinor_field_add_f(&s[0],r1,Mo);

    /* assign r<-r1 */
    sptmp=r;
    r=r1;
    r1=sptmp; /*exchange pointers */

    /* update phi */
    M(Ms,&s[0]);
    ctmp2 = spinor_field_prod_f(o0, Ms); /* in=o al passo 0 (see above) */
    _complex_div(phi,ctmp2,delta);

    _print_complex(delta);
    _print_complex(phi);
    _print_complex(alpha[0]);
    _print_complex(beta[0]);
    _print_complex(chi[0]);

    for (i=1; i<(par->n); ++i) { /* update shifted quantities i=1..n-1 */
      if(sflags[i]) {
	/* compute z3[i]/z2[i] */
	_compute_z(z3+i,z1+i,z2+i,&ctmp1,beta,&oldalpha,&(par->shift[i-1]));
	_print_complex(z3[i]);
	/*compute beta[i] */
	_complex_mul(beta[i],beta[0],z3[i]);
	_print_complex(beta[i]);
	/*compute chi[i] */
	_complex_mulr(ctmp2,-par->shift[i-1],chi[0]);
	_complex_add_1(ctmp2);
	_complex_mul(ctmp3,z3[i],z3[i]); /* needed for alpha[i] */
	_complex_div(chi[i],chi[0],ctmp2);
	_print_complex(chi[i]);
	/* compute alpha[i] */
	_complex_mul(alpha[i],alpha[0],ctmp3);
	_print_complex(alpha[i]);
	/* compute z3[i] */
	_complex_mul(ctmp2,z3[i],z2[i]);
	z3[i]=ctmp2;
	_print_complex(z3[i]);
	/* update solution */
	_complex_mul(ctmp2,z3[i],rho[i]);
	_complex_mul(ctmp3,ctmp2,chi[i]);
	_complex_minus(ctmp4,beta[i]);
	spinor_field_clc_add_assign_f(&out[i],ctmp4,&s[i],ctmp3,o);
	/*spinor_field_clc_add_assign_f(Mo,ctmp4,s[i],ctmp3,o);
	spinor_field_add_f(out[i],out[0],Mo);*/
	/* update s[i] */
	_complex_div(ctmp4,chi[i],beta[i]);
	_complex_mul(ctmp3,ctmp4,ctmp2);
	_complex_minus(ctmp3,ctmp3);
	_complex_mul(ctmp2,z2[i],rho[i]);
	_complex_mul(ctmp5,ctmp2,ctmp4);
	spinor_field_clc_add_assign_f(&s[i],ctmp3,o,ctmp5,r1); /* not done yet */
	
	ctmp3=rho[i];
	_complex_mulr(ctmp2,-par->shift[i-1],chi[0]);
	_complex_add_1(ctmp2);
	_complex_div(rho[i],ctmp3,ctmp2); /* update rho[i] */
	_print_complex(rho[i]);

	_complex_mul(ctmp2,z3[i],rho[i]);
	spinor_field_clc_f(Mo,ctmp2,r,alpha[i],&s[i]); /* use Mo as temporary storage */
	 spinor_field_copy_f(&s[i],Mo); 
	/* change pointers instead */
	/*
	  sptmp=s[i];
	  s[i]=Mo;
	  Mo=sptmp;
	*/

	 /* shift z2<-z3; z1<-z2 */
	 z2[i]=z3[i];
	 z1[i]=z2[i];

      }
    }

    rtmp1=spinor_field_sqnorm_f(r);

    if(rtmp1<par->err2)
      notconverged=0;
	
    printf("[ %d ] residuo=%e\n",cgiter,rtmp1);

    /* Uncomment this to print cg recursion parameters
       printf("[ %d ] alpha=%e\n",cgiter,alpha);
       printf("[ %d ] omega=%e\n",cgiter,omega);
       printf("[ %d ] still runnning=%d\n",cgiter,notconverged);
       for (i=0;i<par->n;++i) printf("z3[%d]=%e; ",i,z3[i]);
       printf("\n[ %d ] gamma=%e\n",cgiter,gamma);
       printf("[ %d ] delta=%e\n",cgiter,delta);
    */
      
  } while ((par->max_iter==0 || cgiter<par->max_iter) && notconverged);

  /* test results */
#ifndef NDEBUG
  for(i=0;i<par->n;++i){
    double norm;
    M(Ms,out[i]);
    if(i!=0) {
      spinor_field_mul_add_assign_f(Ms,-par->shift[i-1],out[i]);
    }
    spinor_field_mul_add_assign_f(Ms,-1.0,in);
    norm=spinor_field_sqnorm_f(Ms);
    if (fabs(norm)>5.*par->err2)
      printf("BiCGstab Failed: err2[%d] = %e\n",i,norm);
  }
#endif
   
  /* free memory */
  free_spinor_field_f(s);
  free(z1);
  free(sflags);

  /* return number of cg iter */
  return cgiter;
}
