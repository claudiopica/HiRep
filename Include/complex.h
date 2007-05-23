/*******************************************************************************
*
* File complex.h
*
* Type definitions and macros for complex numbers  
*
*******************************************************************************/

#ifndef COMPLEX_H
#define COMPLEX_H

/*******************************************************************************
*
* Definitions of type complex
*
*******************************************************************************/

typedef struct
{
   float re,im;
} complex;

typedef struct
{
   double re,im;
} complex_dble;

/*******************************************************************************
*
* Macros for complex numbers
*
* Arguments are variables of type complex
*
*******************************************************************************/


/*
* a=0 (a complex)
*/
#define _complex_0(a) \
   (a).re=0.;	\
   (a).im=0.

/*
* a=1. (a complex)
*/
#define _complex_1(a) \
   (a).re=1.;	\
   (a).im=0.

/*
* a+=1. (a complex)
*/
#define _complex_add_1(a) \
   (a).re+=1.

/*
* a=i (a complex)
*/
#define _complex_i(a) \
   (a).re=0.;	\
   (a).im=1.

/*
* a=b*c (a,b,c complex)
*/
#define _complex_mul(a,b,c) \
   (a).re=((b).re*(c).re-(b).im*(c).im); \
   (a).im=((b).re*(c).im+(b).im*(c).re)

/*
* a=r*b (a,b complex; r real)
*/
#define _complex_mulr(a,r,b) \
   (a).re=(r)*(b).re; \
   (a).im=(r)*(b).im

/*
* a=b+c (a,b,c complex)
*/
#define _complex_add(a,b,c) \
   (a).re=((b).re+(c).re); \
   (a).im=((b).im+(c).im)

/*
* a=b-c (a,b,c complex)
*/
#define _complex_sub(a,b,c) \
   (a).re=((b).re-(c).re); \
   (a).im=((b).im-(c).im)

/*
* a=b/c (a,b,c complex)
*/
#define _complex_div(a,b,c) \
  (a).im=((c).re*(c).re+(c).im*(c).im);		   \
  (a).re=(((b).re*(c).re+(b).im*(c).im)/(a).im);   \
  (a).im=(((b).im*(c).re-(b).re*(c).im)/(a).im)

/*
* a=1/b (a,b complex)
*/
#define _complex_inv(a,b) \
   (a).im=((b).re*(b).re+(b).im*(b).im); \
   (a).re=((b).re/(a).im); \
   (a).im=(-(b).im/(a).im)

/*
* Re(a^*b) (a,b,c complex)
*/
#define _complex_prod_re(a,b) \
   ((a).re*(b).re+(a).im*(b).im)

/*
* Im(a^*b) (a,b,c complex)
*/
#define _complex_prod_im(a,b) \
   ((a).re*(b).im-(a).im*(b).re)

/*
* a=-b (a complex)
*/
#define _complex_minus(a,b) \
   (a).re=-(b).re; \
   (a).im=-(b).im

/*
* a=b+i*c (a,b,c complex)
*/
#define _complex_i_add(a,b,c) \
   (a).re=((b).re-(c).im); \
   (a).im=((b).im+(c).re)

/*
* a=b-i*c (a,b,c complex)
*/
#define _complex_i_sub(a,b,c) \
   (a).re=((b).re+(c).im); \
   (a).im=((b).im-(c).re)

/*
* a+=b (a,b complex)
*/
#define _complex_add_assign(a,b) \
   (a).re+=(b).re; \
   (a).im+=(b).im

/*
* a-=b (a,b complex)
*/
#define _complex_sub_assign(a,b) \
   (a).re-=(b).re; \
   (a).im-=(b).im

/*
* a+=i*b (a,b complex)
*/
#define _complex_i_add_assign(a,b) \
   (a).re-=(b).im; \
   (a).im+=(b).re

/*
* a-=i*b (a,b complex)
*/
#define _complex_i_sub_assign(a,b) \
   (a).re+=(b).im; \
   (a).im-=(b).re

/*
* a+=b*c (a,b,c complex)
*/
#define _complex_mul_assign(a,b,c) \
   (a).re+=((b).re*(c).re-(b).im*(c).im); \
   (a).im+=((b).re*(c).im+(b).im*(c).re)

/*
* a+=b*(c^+) (a,b,c complex)
*/
#define _complex_mul_star_assign(a,b,c) \
   (a).re+=((b).re*(c).re+(b).im*(c).im); \
   (a).im+=((b).im*(c).re-(b).re*(c).im)

/*
* a+=Re[ b*(c^+) ] (a real ;b,c complex)
*/
#define _complex_mul_star_assign_re(a,b,c) \
   (a)+=((b).re*(c).re+(b).im*(c).im)

/*
* a-=b*c (a,b,c complex)
*/
#define _complex_mul_sub_assign(a,b,c) \
   (a).re-=((b).re*(c).re-(b).im*(c).im); \
   (a).im-=((b).re*(c).im+(b).im*(c).re)

/*
* a+=r*b (a,b complex; r real)
*/
#define _complex_mulr_assign(a,r,b) \
   (a).re+=(r)*(b).re; \
   (a).im+=(r)*(b).im


/*
* a=r1*c1+r2*c2 (a,c1,c2 complex; r1,r2 real)
*/
#define _complex_rlc(a,r1,c1,r2,c2) \
    _complex_mulr(a,r1,c1); \
    _complex_mulr_assign(a,r2,c2);

/*
* a+=r1*c1+r2*c2 (a,c1,c2 complex; r1,r2 real)
*/
#define _complex_rlc_assign(a,r1,c1,r2,c2) \
    _complex_mulr_assign(a,r1,c1); \
    _complex_mulr_assign(a,r2,c2);

/*
* a=z1*c1+z2*c2 (a,z1,c1,z2,c2 complex)
*/
#define _complex_clc(a,z1,c1,z2,c2) \
    _complex_mul(a,z1,c1); \
    _complex_mul_assign(a,z2,c2);

/*
* a+=z1*c1+z2*c2 (a,z1,c1,z2,c2 complex)
*/
#define _complex_clc_assign(a,z1,c1,z2,c2) \
    _complex_mul_assign(a,z1,c1); \
    _complex_mul_assign(a,z2,c2);


#endif
