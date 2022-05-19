/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *
* All rights reserved.                                                      *
\***************************************************************************/

/*******************************************************************************
*
* File gpu_complex.c
*
* Type definitions and macros for complex numbers used in C++ and CUDA
*
*******************************************************************************/

/*******************************************************************************
*
* Definitions of type complex
* The following structs are generated using:
* https://github.com/erikkjellgren/SimplyComplex
* Do NOT change them manually
*
*******************************************************************************/
#include "hr_complex.h"
__host__ __device__ hr_complex_int::hr_complex_int(void){}
__host__ __device__ hr_complex_int::hr_complex_int(const int x)
{
  re = x;
  im = 0;
}
__host__ __device__ hr_complex_int::hr_complex_int(const float x)
{
  re = x;
  im = 0;
}
__host__ __device__ hr_complex_int::hr_complex_int(const double x)
{
  re = x;
  im = 0;
}
__host__ __device__ hr_complex_int::hr_complex_int(const int a, const int b)
{
  re = a;
  im = b;
}
__host__ __device__ hr_complex_int::hr_complex_int(const int a, const float b)
{
  re = a;
  im = b;
}
__host__ __device__ hr_complex_int::hr_complex_int(const int a, const double b)
{
  re = a;
  im = b;
}
__host__ __device__ hr_complex_int::hr_complex_int(const float a, const int b)
{
  re = a;
  im = b;
}
__host__ __device__ hr_complex_int::hr_complex_int(const float a, const float b)
{
  re = a;
  im = b;
}
__host__ __device__ hr_complex_int::hr_complex_int(const float a, const double b)
{
  re = a;
  im = b;
}
__host__ __device__ hr_complex_int::hr_complex_int(const double a, const int b)
{
  re = a;
  im = b;
}
__host__ __device__ hr_complex_int::hr_complex_int(const double a, const float b)
{
  re = a;
  im = b;
}
__host__ __device__ hr_complex_int::hr_complex_int(const double a, const double b)
{
  re = a;
  im = b;
}
__host__ __device__ hr_complex_int& hr_complex_int::operator=(const int x)
{
  re = x;
  im = 0;
  return *this;
}
__host__ __device__ hr_complex_int& hr_complex_int::operator=(const float x)
{
  re = x;
  im = 0;
  return *this;
}
__host__ __device__ hr_complex_int& hr_complex_int::operator=(const double x)
{
  re = x;
  im = 0;
  return *this;
}
__host__ __device__ hr_complex_int& hr_complex_int::operator=(const hr_complex_int &x)
{
  re = x.re;
  im = x.im;
  return *this;
}
__host__ __device__ hr_complex_int hr_complex_int::operator-(void) const
{
  return hr_complex_int(-re, -im);
}
__host__ __device__ hr_complex_int& hr_complex_int::operator+=(const int x)
{
  re += x;
  return *this;
}
__host__ __device__ hr_complex_int& hr_complex_int::operator+=(const float x)
{
  re += x;
  return *this;
}
__host__ __device__ hr_complex_int& hr_complex_int::operator+=(const double x)
{
  re += x;
  return *this;
}
__host__ __device__ hr_complex_int& hr_complex_int::operator+=(const hr_complex_int &x)
{
  re += x.re;
  im += x.im;
  return *this;
}
__host__ __device__ hr_complex_int& hr_complex_int::operator-=(const int x)
{
  re -= x;
  return *this;
}
__host__ __device__ hr_complex_int& hr_complex_int::operator-=(const float x)
{
  re -= x;
  return *this;
}
__host__ __device__ hr_complex_int& hr_complex_int::operator-=(const double x)
{
  re -= x;
  return *this;
}
__host__ __device__ hr_complex_int& hr_complex_int::operator-=(const hr_complex_int &x)
{
  re -= x.re;
  im -= x.im;
  return *this;
}
__host__ __device__ hr_complex_int& hr_complex_int::operator*=(const int x)
{
  re *= x;
  im *= x;
  return *this;
}
__host__ __device__ hr_complex_int& hr_complex_int::operator*=(const float x)
{
  re *= x;
  im *= x;
  return *this;
}
__host__ __device__ hr_complex_int& hr_complex_int::operator*=(const double x)
{
  re *= x;
  im *= x;
  return *this;
}
__host__ __device__ hr_complex_int& hr_complex_int::operator*=(const hr_complex_int &x)
{
  int re_tmp, im_tmp;
  re_tmp = re*(x.re) - im*(x.im);
  im_tmp = re*(x.im) + im*(x.re);
  re = re_tmp;
  im = im_tmp;
  return *this;
}
__host__ __device__ hr_complex_int& hr_complex_int::operator/=(const int x)
{
  re /= x;
  im /= x;
  return *this;
}
__host__ __device__ hr_complex_int& hr_complex_int::operator/=(const float x)
{
  re /= x;
  im /= x;
  return *this;
}
__host__ __device__ hr_complex_int& hr_complex_int::operator/=(const double x)
{
  re /= x;
  im /= x;
  return *this;
}
__host__ __device__ hr_complex_int& hr_complex_int::operator/=(const hr_complex_int &x)
{
  int re_tmp, im_tmp;
  re_tmp = (re*(x.re) + im*(x.im)) / ((x.re)*(x.re) + (x.im)*(x.im));
  im_tmp = (im*(x.re) - re*(x.im)) / ((x.re)*(x.re) + (x.im)*(x.im));
  re = re_tmp;
  im = im_tmp;
  return *this;
}
__host__ __device__ hr_complex_int hr_complex_int::operator+(const int x) const
{
  return hr_complex_int(re+x, im);
}
__host__ __device__ hr_complex_int hr_complex_int::operator+(const float x) const
{
  return hr_complex_int(re+x, im);
}
__host__ __device__ hr_complex_int hr_complex_int::operator+(const double x) const
{
  return hr_complex_int(re+x, im);
}
__host__ __device__ hr_complex_int hr_complex_int::operator+(const hr_complex_int &x) const
{
  return hr_complex_int(re+x.re, im+x.im);
}
__host__ __device__ hr_complex_int hr_complex_int::operator-(const int x) const
{
  return hr_complex_int(re-x, im);
}
__host__ __device__ hr_complex_int hr_complex_int::operator-(const float x) const
{
  return hr_complex_int(re-x, im);
}
__host__ __device__ hr_complex_int hr_complex_int::operator-(const double x) const
{
  return hr_complex_int(re-x, im);
}
__host__ __device__ hr_complex_int hr_complex_int::operator-(const hr_complex_int &x) const
{
  return hr_complex_int(re-x.re, im-x.im);
}
__host__ __device__ hr_complex_int hr_complex_int::operator*(const int x) const
{
  return hr_complex_int(re*x, im*x);
}
__host__ __device__ hr_complex_int hr_complex_int::operator*(const float x) const
{
  return hr_complex_int(re*x, im*x);
}
__host__ __device__ hr_complex_int hr_complex_int::operator*(const double x) const
{
  return hr_complex_int(re*x, im*x);
}
__host__ __device__ hr_complex_int hr_complex_int::operator*(const hr_complex_int &x) const
{
  return hr_complex_int(re*(x.re) - im*(x.im), re*(x.im) + im*(x.re));
}
__host__ __device__ hr_complex_int hr_complex_int::operator/(const int x) const
{
  return hr_complex_int(re/x, im/x);
}
__host__ __device__ hr_complex_int hr_complex_int::operator/(const float x) const
{
  return hr_complex_int(re/x, im/x);
}
__host__ __device__ hr_complex_int hr_complex_int::operator/(const double x) const
{
  return hr_complex_int(re/x, im/x);
}
__host__ __device__ hr_complex_int hr_complex_int::operator/(const hr_complex_int &x) const
{
  return hr_complex_int((re*(x.re) + im*(x.im)) / ((x.re)*(x.re) + (x.im)*(x.im)), (im*(x.re) - re*(x.im)) / ((x.re)*(x.re) + (x.im)*(x.im)));
}
__host__ __device__ hr_complex_flt::hr_complex_flt(void){}
__host__ __device__ hr_complex_flt::hr_complex_flt(const int x)
{
  re = x;
  im = 0;
}
__host__ __device__ hr_complex_flt::hr_complex_flt(const float x)
{
  re = x;
  im = 0;
}
__host__ __device__ hr_complex_flt::hr_complex_flt(const double x)
{
  re = x;
  im = 0;
}
__host__ __device__ hr_complex_flt::hr_complex_flt(const int a, const int b)
{
  re = a;
  im = b;
}
__host__ __device__ hr_complex_flt::hr_complex_flt(const int a, const float b)
{
  re = a;
  im = b;
}
__host__ __device__ hr_complex_flt::hr_complex_flt(const int a, const double b)
{
  re = a;
  im = b;
}
__host__ __device__ hr_complex_flt::hr_complex_flt(const float a, const int b)
{
  re = a;
  im = b;
}
__host__ __device__ hr_complex_flt::hr_complex_flt(const float a, const float b)
{
  re = a;
  im = b;
}
__host__ __device__ hr_complex_flt::hr_complex_flt(const float a, const double b)
{
  re = a;
  im = b;
}
__host__ __device__ hr_complex_flt::hr_complex_flt(const double a, const int b)
{
  re = a;
  im = b;
}
__host__ __device__ hr_complex_flt::hr_complex_flt(const double a, const float b)
{
  re = a;
  im = b;
}
__host__ __device__ hr_complex_flt::hr_complex_flt(const double a, const double b)
{
  re = a;
  im = b;
}
__host__ __device__ hr_complex_flt::hr_complex_flt(const hr_complex_int &x)
{
  re = x.re;
  im = x.im;
}
__host__ __device__ hr_complex_flt& hr_complex_flt::operator=(const int x)
{
  re = x;
  im = 0;
  return *this;
}
__host__ __device__ hr_complex_flt& hr_complex_flt::operator=(const float x)
{
  re = x;
  im = 0;
  return *this;
}
__host__ __device__ hr_complex_flt& hr_complex_flt::operator=(const double x)
{
  re = x;
  im = 0;
  return *this;
}
__host__ __device__ hr_complex_flt& hr_complex_flt::operator=(const hr_complex_int &x)
{
  re = x.re;
  im = x.im;
  return *this;
}
__host__ __device__ hr_complex_flt& hr_complex_flt::operator=(const hr_complex_flt &x)
{
  re = x.re;
  im = x.im;
  return *this;
}
__host__ __device__ hr_complex_flt hr_complex_flt::operator-(void) const
{
  return hr_complex_flt(-re, -im);
}
__host__ __device__ hr_complex_flt& hr_complex_flt::operator+=(const int x)
{
  re += x;
  return *this;
}
__host__ __device__ hr_complex_flt& hr_complex_flt::operator+=(const float x)
{
  re += x;
  return *this;
}
__host__ __device__ hr_complex_flt& hr_complex_flt::operator+=(const double x)
{
  re += x;
  return *this;
}
__host__ __device__ hr_complex_flt& hr_complex_flt::operator+=(const hr_complex_int &x)
{
  re += x.re;
  im += x.im;
  return *this;
}
__host__ __device__ hr_complex_flt& hr_complex_flt::operator+=(const hr_complex_flt &x)
{
  re += x.re;
  im += x.im;
  return *this;
}
__host__ __device__ hr_complex_flt& hr_complex_flt::operator-=(const int x)
{
  re -= x;
  return *this;
}
__host__ __device__ hr_complex_flt& hr_complex_flt::operator-=(const float x)
{
  re -= x;
  return *this;
}
__host__ __device__ hr_complex_flt& hr_complex_flt::operator-=(const double x)
{
  re -= x;
  return *this;
}
__host__ __device__ hr_complex_flt& hr_complex_flt::operator-=(const hr_complex_int &x)
{
  re -= x.re;
  im -= x.im;
  return *this;
}
__host__ __device__ hr_complex_flt& hr_complex_flt::operator-=(const hr_complex_flt &x)
{
  re -= x.re;
  im -= x.im;
  return *this;
}
__host__ __device__ hr_complex_flt& hr_complex_flt::operator*=(const int x)
{
  re *= x;
  im *= x;
  return *this;
}
__host__ __device__ hr_complex_flt& hr_complex_flt::operator*=(const float x)
{
  re *= x;
  im *= x;
  return *this;
}
__host__ __device__ hr_complex_flt& hr_complex_flt::operator*=(const double x)
{
  re *= x;
  im *= x;
  return *this;
}
__host__ __device__ hr_complex_flt& hr_complex_flt::operator*=(const hr_complex_int &x)
{
  float re_tmp, im_tmp;
  re_tmp = re*(x.re) - im*(x.im);
  im_tmp = re*(x.im) + im*(x.re);
  re = re_tmp;
  im = im_tmp;
  return *this;
}
__host__ __device__ hr_complex_flt& hr_complex_flt::operator*=(const hr_complex_flt &x)
{
  float re_tmp, im_tmp;
  re_tmp = re*(x.re) - im*(x.im);
  im_tmp = re*(x.im) + im*(x.re);
  re = re_tmp;
  im = im_tmp;
  return *this;
}
__host__ __device__ hr_complex_flt& hr_complex_flt::operator/=(const int x)
{
  re /= x;
  im /= x;
  return *this;
}
__host__ __device__ hr_complex_flt& hr_complex_flt::operator/=(const float x)
{
  re /= x;
  im /= x;
  return *this;
}
__host__ __device__ hr_complex_flt& hr_complex_flt::operator/=(const double x)
{
  re /= x;
  im /= x;
  return *this;
}
__host__ __device__ hr_complex_flt& hr_complex_flt::operator/=(const hr_complex_int &x)
{
  float re_tmp, im_tmp;
  re_tmp = (re*(x.re) + im*(x.im)) / ((x.re)*(x.re) + (x.im)*(x.im));
  im_tmp = (im*(x.re) - re*(x.im)) / ((x.re)*(x.re) + (x.im)*(x.im));
  re = re_tmp;
  im = im_tmp;
  return *this;
}
__host__ __device__ hr_complex_flt& hr_complex_flt::operator/=(const hr_complex_flt &x)
{
  float re_tmp, im_tmp;
  re_tmp = (re*(x.re) + im*(x.im)) / ((x.re)*(x.re) + (x.im)*(x.im));
  im_tmp = (im*(x.re) - re*(x.im)) / ((x.re)*(x.re) + (x.im)*(x.im));
  re = re_tmp;
  im = im_tmp;
  return *this;
}
__host__ __device__ hr_complex_flt hr_complex_flt::operator+(const int x) const
{
  return hr_complex_flt(re+x, im);
}
__host__ __device__ hr_complex_flt hr_complex_flt::operator+(const float x) const
{
  return hr_complex_flt(re+x, im);
}
__host__ __device__ hr_complex_flt hr_complex_flt::operator+(const double x) const
{
  return hr_complex_flt(re+x, im);
}
__host__ __device__ hr_complex_flt hr_complex_flt::operator+(const hr_complex_int &x) const
{
  return hr_complex_flt(re+x.re, im+x.im);
}
__host__ __device__ hr_complex_flt hr_complex_flt::operator+(const hr_complex_flt &x) const
{
  return hr_complex_flt(re+x.re, im+x.im);
}
__host__ __device__ hr_complex_flt hr_complex_flt::operator-(const int x) const
{
  return hr_complex_flt(re-x, im);
}
__host__ __device__ hr_complex_flt hr_complex_flt::operator-(const float x) const
{
  return hr_complex_flt(re-x, im);
}
__host__ __device__ hr_complex_flt hr_complex_flt::operator-(const double x) const
{
  return hr_complex_flt(re-x, im);
}
__host__ __device__ hr_complex_flt hr_complex_flt::operator-(const hr_complex_int &x) const
{
  return hr_complex_flt(re-x.re, im-x.im);
}
__host__ __device__ hr_complex_flt hr_complex_flt::operator-(const hr_complex_flt &x) const
{
  return hr_complex_flt(re-x.re, im-x.im);
}
__host__ __device__ hr_complex_flt hr_complex_flt::operator*(const int x) const
{
  return hr_complex_flt(re*x, im*x);
}
__host__ __device__ hr_complex_flt hr_complex_flt::operator*(const float x) const
{
  return hr_complex_flt(re*x, im*x);
}
__host__ __device__ hr_complex_flt hr_complex_flt::operator*(const double x) const
{
  return hr_complex_flt(re*x, im*x);
}
__host__ __device__ hr_complex_flt hr_complex_flt::operator*(const hr_complex_int &x) const
{
  return hr_complex_flt(re*(x.re) - im*(x.im), re*(x.im) + im*(x.re));
}
__host__ __device__ hr_complex_flt hr_complex_flt::operator*(const hr_complex_flt &x) const
{
  return hr_complex_flt(re*(x.re) - im*(x.im), re*(x.im) + im*(x.re));
}
__host__ __device__ hr_complex_flt hr_complex_flt::operator/(const int x) const
{
  return hr_complex_flt(re/x, im/x);
}
__host__ __device__ hr_complex_flt hr_complex_flt::operator/(const float x) const
{
  return hr_complex_flt(re/x, im/x);
}
__host__ __device__ hr_complex_flt hr_complex_flt::operator/(const double x) const
{
  return hr_complex_flt(re/x, im/x);
}
__host__ __device__ hr_complex_flt hr_complex_flt::operator/(const hr_complex_int &x) const
{
  return hr_complex_flt((re*(x.re) + im*(x.im)) / ((x.re)*(x.re) + (x.im)*(x.im)), (im*(x.re) - re*(x.im)) / ((x.re)*(x.re) + (x.im)*(x.im)));
}
__host__ __device__ hr_complex_flt hr_complex_flt::operator/(const hr_complex_flt &x) const
{
  return hr_complex_flt((re*(x.re) + im*(x.im)) / ((x.re)*(x.re) + (x.im)*(x.im)), (im*(x.re) - re*(x.im)) / ((x.re)*(x.re) + (x.im)*(x.im)));
}
__host__ __device__ hr_complex::hr_complex(void){}
__host__ __device__ hr_complex::hr_complex(const int x)
{
  re = x;
  im = 0;
}
__host__ __device__ hr_complex::hr_complex(const float x)
{
  re = x;
  im = 0;
}
__host__ __device__ hr_complex::hr_complex(const double x)
{
  re = x;
  im = 0;
}
__host__ __device__ hr_complex::hr_complex(const int a, const int b)
{
  re = a;
  im = b;
}
__host__ __device__ hr_complex::hr_complex(const int a, const float b)
{
  re = a;
  im = b;
}
__host__ __device__ hr_complex::hr_complex(const int a, const double b)
{
  re = a;
  im = b;
}
__host__ __device__ hr_complex::hr_complex(const float a, const int b)
{
  re = a;
  im = b;
}
__host__ __device__ hr_complex::hr_complex(const float a, const float b)
{
  re = a;
  im = b;
}
__host__ __device__ hr_complex::hr_complex(const float a, const double b)
{
  re = a;
  im = b;
}
__host__ __device__ hr_complex::hr_complex(const double a, const int b)
{
  re = a;
  im = b;
}
__host__ __device__ hr_complex::hr_complex(const double a, const float b)
{
  re = a;
  im = b;
}
__host__ __device__ hr_complex::hr_complex(const double a, const double b)
{
  re = a;
  im = b;
}
__host__ __device__ hr_complex::hr_complex(const hr_complex_int &x)
{
  re = x.re;
  im = x.im;
}
__host__ __device__ hr_complex::hr_complex(const hr_complex_flt &x)
{
  re = x.re;
  im = x.im;
}
__host__ __device__ hr_complex& hr_complex::operator=(const int x)
{
  re = x;
  im = 0;
  return *this;
}
__host__ __device__ hr_complex& hr_complex::operator=(const float x)
{
  re = x;
  im = 0;
  return *this;
}
__host__ __device__ hr_complex& hr_complex::operator=(const double x)
{
  re = x;
  im = 0;
  return *this;
}
__host__ __device__ hr_complex& hr_complex::operator=(const hr_complex_int &x)
{
  re = x.re;
  im = x.im;
  return *this;
}
__host__ __device__ hr_complex& hr_complex::operator=(const hr_complex_flt &x)
{
  re = x.re;
  im = x.im;
  return *this;
}
__host__ __device__ hr_complex& hr_complex::operator=(const hr_complex &x)
{
  re = x.re;
  im = x.im;
  return *this;
}
__host__ __device__ hr_complex hr_complex::operator-(void) const
{
  return hr_complex(-re, -im);
}
__host__ __device__ hr_complex& hr_complex::operator+=(const int x)
{
  re += x;
  return *this;
}
__host__ __device__ hr_complex& hr_complex::operator+=(const float x)
{
  re += x;
  return *this;
}
__host__ __device__ hr_complex& hr_complex::operator+=(const double x)
{
  re += x;
  return *this;
}
__host__ __device__ hr_complex& hr_complex::operator+=(const hr_complex_int &x)
{
  re += x.re;
  im += x.im;
  return *this;
}
__host__ __device__ hr_complex& hr_complex::operator+=(const hr_complex_flt &x)
{
  re += x.re;
  im += x.im;
  return *this;
}
__host__ __device__ hr_complex& hr_complex::operator+=(const hr_complex &x)
{
  re += x.re;
  im += x.im;
  return *this;
}
__host__ __device__ hr_complex& hr_complex::operator-=(const int x)
{
  re -= x;
  return *this;
}
__host__ __device__ hr_complex& hr_complex::operator-=(const float x)
{
  re -= x;
  return *this;
}
__host__ __device__ hr_complex& hr_complex::operator-=(const double x)
{
  re -= x;
  return *this;
}
__host__ __device__ hr_complex& hr_complex::operator-=(const hr_complex_int &x)
{
  re -= x.re;
  im -= x.im;
  return *this;
}
__host__ __device__ hr_complex& hr_complex::operator-=(const hr_complex_flt &x)
{
  re -= x.re;
  im -= x.im;
  return *this;
}
__host__ __device__ hr_complex& hr_complex::operator-=(const hr_complex &x)
{
  re -= x.re;
  im -= x.im;
  return *this;
}
__host__ __device__ hr_complex& hr_complex::operator*=(const int x)
{
  re *= x;
  im *= x;
  return *this;
}
__host__ __device__ hr_complex& hr_complex::operator*=(const float x)
{
  re *= x;
  im *= x;
  return *this;
}
__host__ __device__ hr_complex& hr_complex::operator*=(const double x)
{
  re *= x;
  im *= x;
  return *this;
}
__host__ __device__ hr_complex& hr_complex::operator*=(const hr_complex_int &x)
{
  double re_tmp, im_tmp;
  re_tmp = re*(x.re) - im*(x.im);
  im_tmp = re*(x.im) + im*(x.re);
  re = re_tmp;
  im = im_tmp;
  return *this;
}
__host__ __device__ hr_complex& hr_complex::operator*=(const hr_complex_flt &x)
{
  double re_tmp, im_tmp;
  re_tmp = re*(x.re) - im*(x.im);
  im_tmp = re*(x.im) + im*(x.re);
  re = re_tmp;
  im = im_tmp;
  return *this;
}
__host__ __device__ hr_complex& hr_complex::operator*=(const hr_complex &x)
{
  double re_tmp, im_tmp;
  re_tmp = re*(x.re) - im*(x.im);
  im_tmp = re*(x.im) + im*(x.re);
  re = re_tmp;
  im = im_tmp;
  return *this;
}
__host__ __device__ hr_complex& hr_complex::operator/=(const int x)
{
  re /= x;
  im /= x;
  return *this;
}
__host__ __device__ hr_complex& hr_complex::operator/=(const float x)
{
  re /= x;
  im /= x;
  return *this;
}
__host__ __device__ hr_complex& hr_complex::operator/=(const double x)
{
  re /= x;
  im /= x;
  return *this;
}
__host__ __device__ hr_complex& hr_complex::operator/=(const hr_complex_int &x)
{
  double re_tmp, im_tmp;
  re_tmp = (re*(x.re) + im*(x.im)) / ((x.re)*(x.re) + (x.im)*(x.im));
  im_tmp = (im*(x.re) - re*(x.im)) / ((x.re)*(x.re) + (x.im)*(x.im));
  re = re_tmp;
  im = im_tmp;
  return *this;
}
__host__ __device__ hr_complex& hr_complex::operator/=(const hr_complex_flt &x)
{
  double re_tmp, im_tmp;
  re_tmp = (re*(x.re) + im*(x.im)) / ((x.re)*(x.re) + (x.im)*(x.im));
  im_tmp = (im*(x.re) - re*(x.im)) / ((x.re)*(x.re) + (x.im)*(x.im));
  re = re_tmp;
  im = im_tmp;
  return *this;
}
__host__ __device__ hr_complex& hr_complex::operator/=(const hr_complex &x)
{
  double re_tmp, im_tmp;
  re_tmp = (re*(x.re) + im*(x.im)) / ((x.re)*(x.re) + (x.im)*(x.im));
  im_tmp = (im*(x.re) - re*(x.im)) / ((x.re)*(x.re) + (x.im)*(x.im));
  re = re_tmp;
  im = im_tmp;
  return *this;
}
__host__ __device__ hr_complex hr_complex::operator+(const int x) const
{
  return hr_complex(re+x, im);
}
__host__ __device__ hr_complex hr_complex::operator+(const float x) const
{
  return hr_complex(re+x, im);
}
__host__ __device__ hr_complex hr_complex::operator+(const double x) const
{
  return hr_complex(re+x, im);
}
__host__ __device__ hr_complex hr_complex::operator+(const hr_complex_int &x) const
{
  return hr_complex(re+x.re, im+x.im);
}
__host__ __device__ hr_complex hr_complex::operator+(const hr_complex_flt &x) const
{
  return hr_complex(re+x.re, im+x.im);
}
__host__ __device__ hr_complex hr_complex::operator+(const hr_complex &x) const
{
  return hr_complex(re+x.re, im+x.im);
}
__host__ __device__ hr_complex hr_complex::operator-(const int x) const
{
  return hr_complex(re-x, im);
}
__host__ __device__ hr_complex hr_complex::operator-(const float x) const
{
  return hr_complex(re-x, im);
}
__host__ __device__ hr_complex hr_complex::operator-(const double x) const
{
  return hr_complex(re-x, im);
}
__host__ __device__ hr_complex hr_complex::operator-(const hr_complex_int &x) const
{
  return hr_complex(re-x.re, im-x.im);
}
__host__ __device__ hr_complex hr_complex::operator-(const hr_complex_flt &x) const
{
  return hr_complex(re-x.re, im-x.im);
}
__host__ __device__ hr_complex hr_complex::operator-(const hr_complex &x) const
{
  return hr_complex(re-x.re, im-x.im);
}
__host__ __device__ hr_complex hr_complex::operator*(const int x) const
{
  return hr_complex(re*x, im*x);
}
__host__ __device__ hr_complex hr_complex::operator*(const float x) const
{
  return hr_complex(re*x, im*x);
}
__host__ __device__ hr_complex hr_complex::operator*(const double x) const
{
  return hr_complex(re*x, im*x);
}
__host__ __device__ hr_complex hr_complex::operator*(const hr_complex_int &x) const
{
  return hr_complex(re*(x.re) - im*(x.im), re*(x.im) + im*(x.re));
}
__host__ __device__ hr_complex hr_complex::operator*(const hr_complex_flt &x) const
{
  return hr_complex(re*(x.re) - im*(x.im), re*(x.im) + im*(x.re));
}
__host__ __device__ hr_complex hr_complex::operator*(const hr_complex &x) const
{
  return hr_complex(re*(x.re) - im*(x.im), re*(x.im) + im*(x.re));
}
__host__ __device__ hr_complex hr_complex::operator/(const int x) const
{
  return hr_complex(re/x, im/x);
}
__host__ __device__ hr_complex hr_complex::operator/(const float x) const
{
  return hr_complex(re/x, im/x);
}
__host__ __device__ hr_complex hr_complex::operator/(const double x) const
{
  return hr_complex(re/x, im/x);
}
__host__ __device__ hr_complex hr_complex::operator/(const hr_complex_int &x) const
{
  return hr_complex((re*(x.re) + im*(x.im)) / ((x.re)*(x.re) + (x.im)*(x.im)), (im*(x.re) - re*(x.im)) / ((x.re)*(x.re) + (x.im)*(x.im)));
}
__host__ __device__ hr_complex hr_complex::operator/(const hr_complex_flt &x) const
{
  return hr_complex((re*(x.re) + im*(x.im)) / ((x.re)*(x.re) + (x.im)*(x.im)), (im*(x.re) - re*(x.im)) / ((x.re)*(x.re) + (x.im)*(x.im)));
}
__host__ __device__ hr_complex hr_complex::operator/(const hr_complex &x) const
{
  return hr_complex((re*(x.re) + im*(x.im)) / ((x.re)*(x.re) + (x.im)*(x.im)), (im*(x.re) - re*(x.im)) / ((x.re)*(x.re) + (x.im)*(x.im)));
}
__host__ __device__ hr_complex_int operator+(const int x, const hr_complex_int &c){return c+x;}
__host__ __device__ hr_complex_int operator+(const float x, const hr_complex_int &c){return c+x;}
__host__ __device__ hr_complex_int operator+(const double x, const hr_complex_int &c){return c+x;}
__host__ __device__ hr_complex_int operator-(const int x, const hr_complex_int &c){return -c+x;}
__host__ __device__ hr_complex_int operator-(const float x, const hr_complex_int &c){return -c+x;}
__host__ __device__ hr_complex_int operator-(const double x, const hr_complex_int &c){return -c+x;}
__host__ __device__ hr_complex_int operator*(const int x, const hr_complex_int &c){return c*x;}
__host__ __device__ hr_complex_int operator*(const float x, const hr_complex_int &c){return c*x;}
__host__ __device__ hr_complex_int operator*(const double x, const hr_complex_int &c){return c*x;}
__host__ __device__ hr_complex_int operator/(const int x, const hr_complex_int &c){return ((hr_complex_int)x)/c;}
__host__ __device__ hr_complex_int operator/(const float x, const hr_complex_int &c){return ((hr_complex_int)x)/c;}
__host__ __device__ hr_complex_int operator/(const double x, const hr_complex_int &c){return ((hr_complex_int)x)/c;}
__host__ __device__ hr_complex_flt operator+(const int x, const hr_complex_flt &c){return c+x;}
__host__ __device__ hr_complex_flt operator+(const float x, const hr_complex_flt &c){return c+x;}
__host__ __device__ hr_complex_flt operator+(const double x, const hr_complex_flt &c){return c+x;}
__host__ __device__ hr_complex_flt operator+(const hr_complex_int &x, const hr_complex_flt &c){return c+x;}
__host__ __device__ hr_complex_flt operator-(const int x, const hr_complex_flt &c){return -c+x;}
__host__ __device__ hr_complex_flt operator-(const float x, const hr_complex_flt &c){return -c+x;}
__host__ __device__ hr_complex_flt operator-(const double x, const hr_complex_flt &c){return -c+x;}
__host__ __device__ hr_complex_flt operator-(const hr_complex_int &x, const hr_complex_flt &c){return -c+x;}
__host__ __device__ hr_complex_flt operator*(const int x, const hr_complex_flt &c){return c*x;}
__host__ __device__ hr_complex_flt operator*(const float x, const hr_complex_flt &c){return c*x;}
__host__ __device__ hr_complex_flt operator*(const double x, const hr_complex_flt &c){return c*x;}
__host__ __device__ hr_complex_flt operator*(const hr_complex_int &x, const hr_complex_flt &c){return c*x;}
__host__ __device__ hr_complex_flt operator/(const int x, const hr_complex_flt &c){return ((hr_complex_flt)x)/c;}
__host__ __device__ hr_complex_flt operator/(const float x, const hr_complex_flt &c){return ((hr_complex_flt)x)/c;}
__host__ __device__ hr_complex_flt operator/(const double x, const hr_complex_flt &c){return ((hr_complex_flt)x)/c;}
__host__ __device__ hr_complex_flt operator/(const hr_complex_int &x, const hr_complex_flt &c){return ((hr_complex_flt)x)/c;}
__host__ __device__ hr_complex operator+(const int x, const hr_complex &c){return c+x;}
__host__ __device__ hr_complex operator+(const float x, const hr_complex &c){return c+x;}
__host__ __device__ hr_complex operator+(const double x, const hr_complex &c){return c+x;}
__host__ __device__ hr_complex operator+(const hr_complex_int &x, const hr_complex &c){return c+x;}
__host__ __device__ hr_complex operator+(const hr_complex_flt &x, const hr_complex &c){return c+x;}
__host__ __device__ hr_complex operator-(const int x, const hr_complex &c){return -c+x;}
__host__ __device__ hr_complex operator-(const float x, const hr_complex &c){return -c+x;}
__host__ __device__ hr_complex operator-(const double x, const hr_complex &c){return -c+x;}
__host__ __device__ hr_complex operator-(const hr_complex_int &x, const hr_complex &c){return -c+x;}
__host__ __device__ hr_complex operator-(const hr_complex_flt &x, const hr_complex &c){return -c+x;}
__host__ __device__ hr_complex operator*(const int x, const hr_complex &c){return c*x;}
__host__ __device__ hr_complex operator*(const float x, const hr_complex &c){return c*x;}
__host__ __device__ hr_complex operator*(const double x, const hr_complex &c){return c*x;}
__host__ __device__ hr_complex operator*(const hr_complex_int &x, const hr_complex &c){return c*x;}
__host__ __device__ hr_complex operator*(const hr_complex_flt &x, const hr_complex &c){return c*x;}
__host__ __device__ hr_complex operator/(const int x, const hr_complex &c){return ((hr_complex)x)/c;}
__host__ __device__ hr_complex operator/(const float x, const hr_complex &c){return ((hr_complex)x)/c;}
__host__ __device__ hr_complex operator/(const double x, const hr_complex &c){return ((hr_complex)x)/c;}
__host__ __device__ hr_complex operator/(const hr_complex_int &x, const hr_complex &c){return ((hr_complex)x)/c;}
__host__ __device__ hr_complex operator/(const hr_complex_flt &x, const hr_complex &c){return ((hr_complex)x)/c;}
__host__ __device__  hr_complex_int::hr_complex_int(const hr_complex_flt &x)
{
  re = x.re;
  im = x.im;
}
__host__ __device__  hr_complex_int::hr_complex_int(const hr_complex &x)
{
  re = x.re;
  im = x.im;
}
__host__ __device__ hr_complex_int& hr_complex_int::operator=(const hr_complex_flt &x)
{
  re = x.re;
  im = x.im;
  return *this;
}
__host__ __device__ hr_complex_int& hr_complex_int::operator=(const hr_complex &x)
{
  re = x.re;
  im = x.im;
  return *this;
}
__host__ __device__ hr_complex_int& hr_complex_int::operator+=(const hr_complex_flt &x)
{
  re += x.re;
  im += x.im;
  return *this;
}
__host__ __device__ hr_complex_int& hr_complex_int::operator+=(const hr_complex &x)
{
  re += x.re;
  im += x.im;
  return *this;
}
__host__ __device__ hr_complex_int& hr_complex_int::operator-=(const hr_complex_flt &x)
{
  re -= x.re;
  im -= x.im;
  return *this;
}
__host__ __device__ hr_complex_int& hr_complex_int::operator-=(const hr_complex &x)
{
  re -= x.re;
  im -= x.im;
  return *this;
}
__host__ __device__ hr_complex_int& hr_complex_int::operator*=(const hr_complex_flt &x)
{
  int re_tmp, im_tmp;
  re_tmp = re*(x.re) - im*(x.im);
  im_tmp = re*(x.im) + im*(x.re);
  re = re_tmp;
  im = im_tmp;
  return *this;
}
__host__ __device__ hr_complex_int& hr_complex_int::operator*=(const hr_complex &x)
{
  int re_tmp, im_tmp;
  re_tmp = re*(x.re) - im*(x.im);
  im_tmp = re*(x.im) + im*(x.re);
  re = re_tmp;
  im = im_tmp;
  return *this;
}
__host__ __device__ hr_complex_int& hr_complex_int::operator/=(const hr_complex_flt &x)
{
  int re_tmp, im_tmp;
  re_tmp = (re*(x.re) + im*(x.im)) / ((x.re)*(x.re) + (x.im)*(x.im));
  im_tmp = (im*(x.re) - re*(x.im)) / ((x.re)*(x.re) + (x.im)*(x.im));
  re = re_tmp;
  im = im_tmp;
  return *this;
}
__host__ __device__ hr_complex_int& hr_complex_int::operator/=(const hr_complex &x)
{
  int re_tmp, im_tmp;
  re_tmp = (re*(x.re) + im*(x.im)) / ((x.re)*(x.re) + (x.im)*(x.im));
  im_tmp = (im*(x.re) - re*(x.im)) / ((x.re)*(x.re) + (x.im)*(x.im));
  re = re_tmp;
  im = im_tmp;
  return *this;
}
__host__ __device__  hr_complex_flt::hr_complex_flt(const hr_complex &x)
{
  re = x.re;
  im = x.im;
}
__host__ __device__ hr_complex_flt& hr_complex_flt::operator=(const hr_complex &x)
{
  re = x.re;
  im = x.im;
  return *this;
}
__host__ __device__ hr_complex_flt& hr_complex_flt::operator+=(const hr_complex &x)
{
  re += x.re;
  im += x.im;
  return *this;
}
__host__ __device__ hr_complex_flt& hr_complex_flt::operator-=(const hr_complex &x)
{
  re -= x.re;
  im -= x.im;
  return *this;
}
__host__ __device__ hr_complex_flt& hr_complex_flt::operator*=(const hr_complex &x)
{
  float re_tmp, im_tmp;
  re_tmp = re*(x.re) - im*(x.im);
  im_tmp = re*(x.im) + im*(x.re);
  re = re_tmp;
  im = im_tmp;
  return *this;
}
__host__ __device__ hr_complex_flt& hr_complex_flt::operator/=(const hr_complex &x)
{
  float re_tmp, im_tmp;
  re_tmp = (re*(x.re) + im*(x.im)) / ((x.re)*(x.re) + (x.im)*(x.im));
  im_tmp = (im*(x.re) - re*(x.im)) / ((x.re)*(x.re) + (x.im)*(x.im));
  re = re_tmp;
  im = im_tmp;
  return *this;
}
