/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *
* All rights reserved.                                                      *
\***************************************************************************/

/*******************************************************************************
*
* File gpu_complex.h
*
* Type definitions and macros for complex numbers used in C++ and CUDA
*
*******************************************************************************/

#ifndef COMPLEX_H
#define COMPLEX_H

/*******************************************************************************
*
* Definitions of type complex
*
*******************************************************************************/
struct hr_complex;
struct hr_complex_flt
{
  float re, im;

// Equal
  __host__ __device__ hr_complex_flt operator=(const int x)
  {
    re = (float)x;
    im = 0.0;
    return *this;
  }
  __host__ __device__ hr_complex_flt operator=(const hr_complex x);

// Negate
  __host__ __device__ hr_complex_flt operator-(void)
  {
    re = -re;
    im = -im;
    return *this;
  }

// Addition
  __host__ __device__ hr_complex_flt operator+(const float x)
  {
    re += x;
    return *this;
  }
  __host__ __device__ hr_complex_flt operator+(const hr_complex_flt x)
  {
    re += x.re;
    im += x.im;
    return *this;
  }
  __host__ __device__ hr_complex_flt operator+=(const float x)
  {
    re += x;
    return *this;
  }
  __host__ __device__ hr_complex_flt operator+=(const hr_complex_flt x)
  {
    re += x.re;
    im += x.im;
    return *this;
  }
  __host__ __device__ hr_complex_flt operator+=(const hr_complex x);

// Subtraction
  __host__ __device__ hr_complex_flt operator-(const float x)
  {
    re -= x;
    return *this;
  }
  __host__ __device__ hr_complex_flt operator-(const hr_complex_flt x)
  {
    re -= x.re;
    im -= x.im;
    return *this;
  }
  __host__ __device__ hr_complex_flt operator-=(const float x)
  {
    re -= x;
    return *this;
  }
  __host__ __device__ hr_complex_flt operator-=(const hr_complex_flt x)
  {
    re -= x.re;
    im -= x.im;
    return *this;
  }

// Multiplication
  __host__ __device__ hr_complex_flt operator*(const float x)
  {
    re *= x;
    return *this;
  }
  __host__ __device__ hr_complex_flt operator*(const hr_complex_flt x)
  {
    re = re*x.re - im*x.im;
    im = re*x.im + im*x.re;
    return *this;
  }
};
// Needed to make the operators work for both float O hr_complex_flt and hr_complex_flt O float
__host__ __device__ hr_complex_flt operator*(const float x, hr_complex_flt c){return c*x;}
__host__ __device__ hr_complex_flt operator+(const float x, hr_complex_flt c){return c+x;}
__host__ __device__ hr_complex_flt operator-(const float x, hr_complex_flt c){return c-x;}

struct hr_complex
{
  double re, im;
// Constructors
  __host__ __device__ hr_complex(void){}
  __host__ __device__ hr_complex(const double a, const double b)
  {
    re = a;
    im = b;
  }
  __host__ __device__ hr_complex(const int x)
  {
    re = (double)x;
    im = 0.0;
  }

// Equal
  __host__ __device__ hr_complex operator=(const int x)
  {
    re = (double)x;
    im = 0.0;
    return *this;
  }
  __host__ __device__ hr_complex operator=(const hr_complex_flt x)
  {
    re = (double)x.re;
    im = (double)x.im;
    return *this;
  }

// Negate
  __host__ __device__ hr_complex operator-(void)
  {
    re = -re;
    im = -im;
    return *this;
  }

// Addition
  __host__ __device__ hr_complex operator+(const double x)
  {
    re += x;
    return *this;
  }
  __host__ __device__ hr_complex operator+(const hr_complex x)
  {
    re += x.re;
    im += x.im;
    return *this;
  }
  __host__ __device__ hr_complex operator+=(const double x)
  {
    re += x;
    return *this;
  }
  __host__ __device__ hr_complex operator+=(const hr_complex x)
  {
    re += x.re;
    im += x.im;
    return *this;
  }
  __host__ __device__ hr_complex operator+=(const hr_complex_flt x)
  {
    re += (double)x.re;
    im += (double)x.im;
    return *this;
  }

// Subtraction
  __host__ __device__ hr_complex operator-(const double x)
  {
    re -= x;
    return *this;
  }
  __host__ __device__ hr_complex operator-(const hr_complex x)
  {
    re -= x.re;
    im -= x.im;
    return *this;
  }
  __host__ __device__ hr_complex operator-=(const double x)
  {
    re -= x;
    return *this;
  }
  __host__ __device__ hr_complex operator-=(const hr_complex x)
  {
    re -= x.re;
    im -= x.im;
    return *this;
  }

// Multiplication
  __host__ __device__ hr_complex operator*(const double x)
  {
    re *= x;
    return *this;
  }
  __host__ __device__ hr_complex operator*(const hr_complex x)
  {
    re = re*x.re - im*x.im;
    im = re*x.im + im*x.re;
    return *this;
  }
  __host__ __device__ hr_complex operator*(const hr_complex_flt x)
  {
    re = re*(double)x.re - im*(double)x.im;
    im = re*(double)x.im + im*(double)x.re;
    return *this;
  }
};
// Needed to make the operators work for both double O hr_complex and hr_complex O double
__host__ __device__ hr_complex operator*(const double x, hr_complex c){return c*x;}
__host__ __device__ hr_complex operator*(const hr_complex_flt x, hr_complex c){return c*x;}
__host__ __device__ hr_complex operator+(const double x, hr_complex c){return c+x;}
__host__ __device__ hr_complex operator-(const double x, hr_complex c){return c-x;}

// Operators for hr_complex_flt that references hr_complex
__host__ __device__ hr_complex_flt hr_complex_flt::operator=(const hr_complex x)
{
  re = (float)x.re;
  im = (float)x.im;
  return *this;
}
__host__ __device__ hr_complex_flt hr_complex_flt::operator+=(const hr_complex x)
{
  re += (float)x.re;
  im += (float)x.im;
  return *this;
}

#define I (hr_complex(0.0, 1.0))
#define creal(a) (a.re)
#define cimag(a) (a.im)
#define conj(a) ((a.im) = -(a.im))

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
   a.re

/*
* Im(a) (a complex)
*/
#define _complex_im(a) \
   a.im

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
