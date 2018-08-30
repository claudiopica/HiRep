/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *   
* All rights reserved.                                                      * 
\***************************************************************************/

/*******************************************************************************
*
* File hirep_complex.h
*
* Type definitions and macros for complex numbers  
*
*******************************************************************************/

#ifndef HR_COMPLEX_H
#define HR_COMPLEX_H

#include <tgmath.h>
// tgmath includes math.h and complex.h
// and defines type-generic macros for math functions
// e.g: float complex fc; creal(fc) invokes crealf(fc)

/*******************************************************************************
*
* Definitions of type complex
*
*******************************************************************************/

typedef double complex hr_complex;
typedef float complex hr_complex_flt;

// USE STD C99 COMPLEX INSTEAD!!!

/*******************************************************************************
*
* Macros for complex numbers
*
* Arguments are variables of type complex
*
*******************************************************************************/

/*
* Re(a) (a complex)
*/
#define _complex_re(a) \
   creal(a)

/*
* Im(a) (a complex)
*/
#define _complex_im(a) \
   cimag(a)

/*
* a=0 (a complex)
*/
#define _complex_0(a) \
   (a)=0
/*
* a=1. (a complex)
*/
#define _complex_1(a) \
   (a)=1

/*
* a+=1. (a complex)
*/
#define _complex_add_1(a) \
   (a)+=1

/*
* a=i (a complex)
*/
#define _complex_i(a) \
   (a)=I

/*
* a=b^+ (a,b complex)
*/
#define _complex_star(a,b) \
   (a)=conj(b)
/*
* a=-b^+ (a,b complex)
*/
#define _complex_star_minus(a,b) \
   (a)=-conj(b)
/*
* a=a^+ (a complex)
*/
#define _complex_star_assign(a) \
   (a)=conj(a)

/*
* a=b*c (a,b,c complex)
*/
#define _complex_mul(a,b,c) \
   (a)=(b)*(c)
/*
* a=r*b (a,b complex; r real)
*/
#define _complex_mulr(a,r,b) \
   (a)=(r)*(b)

/*
* a=b+c (a,b,c complex)
*/
#define _complex_add(a,b,c) \
   (a)=(b)+(c)

/*
* a=b-c (a,b,c complex)
*/
#define _complex_sub(a,b,c) \
   (a)=(b)-(c)

/*
* a=b+c^(+) (a,b,c complex)
*/
#define _complex_add_star(a,b,c) \
   (a)=(b)+conj(c)
/*
* a=b-c^(+) (a,b,c complex)
*/
#define _complex_sub_star(a,b,c) \
   (a)=(b)-conj(c)

/*
* a=b/c (a,b,c complex)
*/
#define _complex_div(a,b,c) \
  (a)=(b)/(c)

/*
* a=1/b (a,b complex)
*/
#define _complex_inv(a,b) \
  (a)=1/(b)
  
/*
* a^*b (a,b complex)
*/
#define _complex_prod(a,b) \
   conj(a)*b

/*
* Re(a^*b) (a,b complex)
*/
#define _complex_prod_re(a,b) \
   creal(conj(a)*b)

/*
* Re((1-a)^*(1-b)) (a,b complex)
*/
#define _complex_prod_m1_re(a,b) \
   creal(conj(1-a)*b)

/*
* Im(a^*b) (a,b complex)
*/
#define _complex_prod_im(a,b) \
   cimag(conj(a)*b)

/*
* c+=(a^*b) (a,b,c complex)
*/
#define _complex_prod_assign(c,a,b) \
   (c)+=conj(a)*(b)

/*
* c+=(a^*b^) (a,b,c complex)
*/
#define _complex_mul_star_star_assign(c,a,b) \
   (c)+=conj((a)*(b))

/*
* a=-b (a complex)
*/
#define _complex_minus(a,b) \
   (a)=-(b)

/*
* a=-i*b (a,b complex)
*/
#define _complex_i_minus(a,b) \
   (a)=-I*(b)

/*
* a=i*b (a,b complex)
*/
#define _complex_i_plus(a,b) \
   (a)=I*(b)

/*
* a=b+i*c (a,b,c complex)
*/
#define _complex_i_add(a,b,c) \
   (a)=(b)+I*(c)

/*
* a=b-i*c (a,b,c complex)
*/
#define _complex_i_sub(a,b,c) \
   (a)=(b)-I*(c)

/*
* a+=b (a,b complex)
*/
#define _complex_add_assign(a,b) \
   (a)+=(b)

/*
* a-=b (a,b complex)
*/
#define _complex_sub_assign(a,b) \
   (a)-=(b)

/*
* a+=(b+c^*) (a,b complex)
*/
#define _complex_add_star_assign(a,b,c) \
   (a)+=(b)+conj(c)

/*
* a+=i*b (a,b complex)
*/
#define _complex_i_add_assign(a,b) \
   (a)+=I*(b)

/*
* a-=i*b (a,b complex)
*/
#define _complex_i_sub_assign(a,b) \
   (a)-=I*(b)

/*
* a+=b*c (a,b,c complex)
*/
#define _complex_mul_assign(a,b,c) \
   (a)+=(b)*(c)

/*
* a+=r*b*c (a,b,c complex, r real)
*/
#define _complex_mulcr_assign(a,r,b,c)	  \
  (a)+=(r)*(b)*(c)

/*
* a=b*(c^+) (a,b,c complex)
*/
#define _complex_mul_star(a,b,c) \
   (a)=(b)*conj(c)

/*
* a+=b*(c^+) (a,b,c complex)
*/
#define _complex_mul_star_assign(a,b,c) \
   (a)+=(b)*conj(c)

/*
* a+=Re[ b*(c^+) ] (a real ;b,c complex)
*/
#define _complex_mul_star_assign_re(a,b,c) \
   (a)+=creal((b)*conj(c))

/*
* a-=b*c (a,b,c complex)
*/
#define _complex_mul_sub_assign(a,b,c) \
   (a)-=(b)*(c)

/*
* a+=r*b (a,b complex; r real)
*/
#define _complex_mulr_assign(a,r,b) \
   (a)+=(r)*(b)


/*
* a=r1*c1+r2*c2 (a,c1,c2 complex; r1,r2 real)
*/
#define _complex_rlc(a,r1,c1,r2,c2) \
    (a)=(r1)*(c1)+(r2)*(c2)

/*
* a+=r1*c1+r2*c2 (a,c1,c2 complex; r1,r2 real)
*/
#define _complex_rlc_assign(a,r1,c1,r2,c2) \
    (a)+=(r1)*(c1)+(r2)*(c2)
/*
* a=z1*c1+z2*c2 (a,z1,c1,z2,c2 complex)
*/
#define _complex_clc(a,z1,c1,z2,c2) \
    (a)=(z1)*(c1)+(z2)*(c2)
/*
* a+=z1*c1+z2*c2 (a,z1,c1,z2,c2 complex)
*/
#define _complex_clc_assign(a,z1,c1,z2,c2) \
    (a)+=(z1)*(c1)+(z2)*(c2)


#endif
