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
* The following structs are generated using:
* https://github.com/erikkjellgren/SimplyComplex
* Do NOT change them manually
*
*******************************************************************************/
struct hr_complex_int;
struct hr_complex_flt;
struct hr_complex;
struct hr_complex_int{
  int re, im;
  __host__ __device__ hr_complex_int(void){}
  __host__ __device__ hr_complex_int(const int x)
  {
    re = (int)x;
    im = (int)0;
  }
  __host__ __device__ hr_complex_int(const float x)
  {
    re = (int)x;
    im = (int)0;
  }
  __host__ __device__ hr_complex_int(const double x)
  {
    re = (int)x;
    im = (int)0;
  }
  __host__ __device__ hr_complex_int(const int a, const int b)
  {
    re = (int)a;
    im = (int)b;
  }
  __host__ __device__ hr_complex_int(const int a, const float b)
  {
    re = (int)a;
    im = (int)b;
  }
  __host__ __device__ hr_complex_int(const int a, const double b)
  {
    re = (int)a;
    im = (int)b;
  }
  __host__ __device__ hr_complex_int(const float a, const int b)
  {
    re = (int)a;
    im = (int)b;
  }
  __host__ __device__ hr_complex_int(const float a, const float b)
  {
    re = (int)a;
    im = (int)b;
  }
  __host__ __device__ hr_complex_int(const float a, const double b)
  {
    re = (int)a;
    im = (int)b;
  }
  __host__ __device__ hr_complex_int(const double a, const int b)
  {
    re = (int)a;
    im = (int)b;
  }
  __host__ __device__ hr_complex_int(const double a, const float b)
  {
    re = (int)a;
    im = (int)b;
  }
  __host__ __device__ hr_complex_int(const double a, const double b)
  {
    re = (int)a;
    im = (int)b;
  }
  __host__ __device__ hr_complex_int(const hr_complex_flt x);
  __host__ __device__ hr_complex_int(const hr_complex x);
  __host__ __device__ hr_complex_int operator=(const int x)
  {
    re = (int)x;
    im = (int)0;
    return *this;
  }
  __host__ __device__ hr_complex_int operator=(const float x)
  {
    re = (int)x;
    im = (int)0;
    return *this;
  }
  __host__ __device__ hr_complex_int operator=(const double x)
  {
    re = (int)x;
    im = (int)0;
    return *this;
  }
  __host__ __device__ hr_complex_int operator=(const hr_complex_int x)
  {
    re = (int)x.re;
    im = (int)x.im;
    return *this;
  }
  __host__ __device__ hr_complex_int operator=(const hr_complex_flt x);
  __host__ __device__ hr_complex_int operator=(const hr_complex x);
  __host__ __device__ hr_complex_int operator-(void)
  {
    return hr_complex_int(-re, -im);
  }
  __host__ __device__ hr_complex_int operator+=(const int x)
  {
    re += (int)x;
    return *this;
  }
  __host__ __device__ hr_complex_int operator+=(const float x)
  {
    re += (int)x;
    return *this;
  }
  __host__ __device__ hr_complex_int operator+=(const double x)
  {
    re += (int)x;
    return *this;
  }
  __host__ __device__ hr_complex_int operator+=(const hr_complex_int x)
  {
    re += (int)x.re;
    im += (int)x.im;
    return *this;
  }
  __host__ __device__ hr_complex_int operator+=(const hr_complex_flt x);
  __host__ __device__ hr_complex_int operator+=(const hr_complex x);
  __host__ __device__ hr_complex_int operator-=(const int x)
  {
    re -= (int)x;
    return *this;
  }
  __host__ __device__ hr_complex_int operator-=(const float x)
  {
    re -= (int)x;
    return *this;
  }
  __host__ __device__ hr_complex_int operator-=(const double x)
  {
    re -= (int)x;
    return *this;
  }
  __host__ __device__ hr_complex_int operator-=(const hr_complex_int x)
  {
    re -= (int)x.re;
    im -= (int)x.im;
    return *this;
  }
  __host__ __device__ hr_complex_int operator-=(const hr_complex_flt x);
  __host__ __device__ hr_complex_int operator-=(const hr_complex x);
  __host__ __device__ hr_complex_int operator*=(const int x)
  {
    re *= (int)x;
    im *= (int)x;
    return *this;
  }
  __host__ __device__ hr_complex_int operator*=(const float x)
  {
    re *= (int)x;
    im *= (int)x;
    return *this;
  }
  __host__ __device__ hr_complex_int operator*=(const double x)
  {
    re *= (int)x;
    im *= (int)x;
    return *this;
  }
  __host__ __device__ hr_complex_int operator*=(const hr_complex_int x)
  {
    int re_tmp, im_tmp;
    re_tmp = re*((int)x.re) - im*((int)x.im);
    im_tmp = re*((int)x.im) + im*((int)x.re);
    re = re_tmp;
    im = im_tmp;
    return *this;
  }
  __host__ __device__ hr_complex_int operator*=(const hr_complex_flt x);
  __host__ __device__ hr_complex_int operator*=(const hr_complex x);
  __host__ __device__ hr_complex_int operator/=(const int x)
  {
    re /= (int)x;
    im /= (int)x;
    return *this;
  }
  __host__ __device__ hr_complex_int operator/=(const float x)
  {
    re /= (int)x;
    im /= (int)x;
    return *this;
  }
  __host__ __device__ hr_complex_int operator/=(const double x)
  {
    re /= (int)x;
    im /= (int)x;
    return *this;
  }
  __host__ __device__ hr_complex_int operator/=(const hr_complex_int x)
  {
    int re_tmp, im_tmp;
    re_tmp = (re*((int)x.re) + im*((int)x.im)) / (((int)x.re)*((int)x.re) + ((int)x.im)*((int)x.im));
    im_tmp = (im*((int)x.re) - re*((int)x.im)) / (((int)x.re)*((int)x.re) + ((int)x.im)*((int)x.im));
    re = re_tmp;
    im = im_tmp;
    return *this;
  }
  __host__ __device__ hr_complex_int operator/=(const hr_complex_flt x);
  __host__ __device__ hr_complex_int operator/=(const hr_complex x);
  __host__ __device__ hr_complex_int operator+(const int x)
  {
    return hr_complex_int(re+(int)x, im);
  }
  __host__ __device__ hr_complex_int operator+(const float x)
  {
    return hr_complex_int(re+(int)x, im);
  }
  __host__ __device__ hr_complex_int operator+(const double x)
  {
    return hr_complex_int(re+(int)x, im);
  }
  __host__ __device__ hr_complex_int operator+(const hr_complex_int x)
  {
    return hr_complex_int(re+(int)x.re, im+(int)x.im);
  }
  __host__ __device__ hr_complex_int operator-(const int x)
  {
    return hr_complex_int(re-(int)x, im);
  }
  __host__ __device__ hr_complex_int operator-(const float x)
  {
    return hr_complex_int(re-(int)x, im);
  }
  __host__ __device__ hr_complex_int operator-(const double x)
  {
    return hr_complex_int(re-(int)x, im);
  }
  __host__ __device__ hr_complex_int operator-(const hr_complex_int x)
  {
    return hr_complex_int(re-(int)x.re, im-(int)x.im);
  }
  __host__ __device__ hr_complex_int operator*(const int x)
  {
    return hr_complex_int(re*(int)x, im*(int)x);
  }
  __host__ __device__ hr_complex_int operator*(const float x)
  {
    return hr_complex_int(re*(int)x, im*(int)x);
  }
  __host__ __device__ hr_complex_int operator*(const double x)
  {
    return hr_complex_int(re*(int)x, im*(int)x);
  }
  __host__ __device__ hr_complex_int operator*(const hr_complex_int x)
  {
    return hr_complex_int(re*((int)x.re) - im*((int)x.im), re*((int)x.im) + im*((int)x.re));
  }
  __host__ __device__ hr_complex_int operator/(const int x)
  {
    return hr_complex_int(re/(int)x, im/(int)x);
  }
  __host__ __device__ hr_complex_int operator/(const float x)
  {
    return hr_complex_int(re/(int)x, im/(int)x);
  }
  __host__ __device__ hr_complex_int operator/(const double x)
  {
    return hr_complex_int(re/(int)x, im/(int)x);
  }
  __host__ __device__ hr_complex_int operator/(const hr_complex_int x)
  {
    return hr_complex_int((re*((int)x.re) + im*((int)x.im)) / (((int)x.re)*((int)x.re) + ((int)x.im)*((int)x.im)), (im*((int)x.re) - re*((int)x.im)) / (((int)x.re)*((int)x.re) + ((int)x.im)*((int)x.im)));
  }
};
struct hr_complex_flt{
  float re, im;
  __host__ __device__ hr_complex_flt(void){}
  __host__ __device__ hr_complex_flt(const int x)
  {
    re = (float)x;
    im = (float)0;
  }
  __host__ __device__ hr_complex_flt(const float x)
  {
    re = (float)x;
    im = (float)0;
  }
  __host__ __device__ hr_complex_flt(const double x)
  {
    re = (float)x;
    im = (float)0;
  }
  __host__ __device__ hr_complex_flt(const int a, const int b)
  {
    re = (float)a;
    im = (float)b;
  }
  __host__ __device__ hr_complex_flt(const int a, const float b)
  {
    re = (float)a;
    im = (float)b;
  }
  __host__ __device__ hr_complex_flt(const int a, const double b)
  {
    re = (float)a;
    im = (float)b;
  }
  __host__ __device__ hr_complex_flt(const float a, const int b)
  {
    re = (float)a;
    im = (float)b;
  }
  __host__ __device__ hr_complex_flt(const float a, const float b)
  {
    re = (float)a;
    im = (float)b;
  }
  __host__ __device__ hr_complex_flt(const float a, const double b)
  {
    re = (float)a;
    im = (float)b;
  }
  __host__ __device__ hr_complex_flt(const double a, const int b)
  {
    re = (float)a;
    im = (float)b;
  }
  __host__ __device__ hr_complex_flt(const double a, const float b)
  {
    re = (float)a;
    im = (float)b;
  }
  __host__ __device__ hr_complex_flt(const double a, const double b)
  {
    re = (float)a;
    im = (float)b;
  }
  __host__ __device__ hr_complex_flt(const hr_complex_int x)
  {
    re = (float)x.re;
    im = (float)x.im;
  }
  __host__ __device__ hr_complex_flt(const hr_complex x);
  __host__ __device__ hr_complex_flt operator=(const int x)
  {
    re = (float)x;
    im = (float)0;
    return *this;
  }
  __host__ __device__ hr_complex_flt operator=(const float x)
  {
    re = (float)x;
    im = (float)0;
    return *this;
  }
  __host__ __device__ hr_complex_flt operator=(const double x)
  {
    re = (float)x;
    im = (float)0;
    return *this;
  }
  __host__ __device__ hr_complex_flt operator=(const hr_complex_int x)
  {
    re = (float)x.re;
    im = (float)x.im;
    return *this;
  }
  __host__ __device__ hr_complex_flt operator=(const hr_complex_flt x)
  {
    re = (float)x.re;
    im = (float)x.im;
    return *this;
  }
  __host__ __device__ hr_complex_flt operator=(const hr_complex x);
  __host__ __device__ hr_complex_flt operator-(void)
  {
    return hr_complex_flt(-re, -im);
  }
  __host__ __device__ hr_complex_flt operator+=(const int x)
  {
    re += (float)x;
    return *this;
  }
  __host__ __device__ hr_complex_flt operator+=(const float x)
  {
    re += (float)x;
    return *this;
  }
  __host__ __device__ hr_complex_flt operator+=(const double x)
  {
    re += (float)x;
    return *this;
  }
  __host__ __device__ hr_complex_flt operator+=(const hr_complex_int x)
  {
    re += (float)x.re;
    im += (float)x.im;
    return *this;
  }
  __host__ __device__ hr_complex_flt operator+=(const hr_complex_flt x)
  {
    re += (float)x.re;
    im += (float)x.im;
    return *this;
  }
  __host__ __device__ hr_complex_flt operator+=(const hr_complex x);
  __host__ __device__ hr_complex_flt operator-=(const int x)
  {
    re -= (float)x;
    return *this;
  }
  __host__ __device__ hr_complex_flt operator-=(const float x)
  {
    re -= (float)x;
    return *this;
  }
  __host__ __device__ hr_complex_flt operator-=(const double x)
  {
    re -= (float)x;
    return *this;
  }
  __host__ __device__ hr_complex_flt operator-=(const hr_complex_int x)
  {
    re -= (float)x.re;
    im -= (float)x.im;
    return *this;
  }
  __host__ __device__ hr_complex_flt operator-=(const hr_complex_flt x)
  {
    re -= (float)x.re;
    im -= (float)x.im;
    return *this;
  }
  __host__ __device__ hr_complex_flt operator-=(const hr_complex x);
  __host__ __device__ hr_complex_flt operator*=(const int x)
  {
    re *= (float)x;
    im *= (float)x;
    return *this;
  }
  __host__ __device__ hr_complex_flt operator*=(const float x)
  {
    re *= (float)x;
    im *= (float)x;
    return *this;
  }
  __host__ __device__ hr_complex_flt operator*=(const double x)
  {
    re *= (float)x;
    im *= (float)x;
    return *this;
  }
  __host__ __device__ hr_complex_flt operator*=(const hr_complex_int x)
  {
    float re_tmp, im_tmp;
    re_tmp = re*((float)x.re) - im*((float)x.im);
    im_tmp = re*((float)x.im) + im*((float)x.re);
    re = re_tmp;
    im = im_tmp;
    return *this;
  }
  __host__ __device__ hr_complex_flt operator*=(const hr_complex_flt x)
  {
    float re_tmp, im_tmp;
    re_tmp = re*((float)x.re) - im*((float)x.im);
    im_tmp = re*((float)x.im) + im*((float)x.re);
    re = re_tmp;
    im = im_tmp;
    return *this;
  }
  __host__ __device__ hr_complex_flt operator*=(const hr_complex x);
  __host__ __device__ hr_complex_flt operator/=(const int x)
  {
    re /= (float)x;
    im /= (float)x;
    return *this;
  }
  __host__ __device__ hr_complex_flt operator/=(const float x)
  {
    re /= (float)x;
    im /= (float)x;
    return *this;
  }
  __host__ __device__ hr_complex_flt operator/=(const double x)
  {
    re /= (float)x;
    im /= (float)x;
    return *this;
  }
  __host__ __device__ hr_complex_flt operator/=(const hr_complex_int x)
  {
    float re_tmp, im_tmp;
    re_tmp = (re*((float)x.re) + im*((float)x.im)) / (((float)x.re)*((float)x.re) + ((float)x.im)*((float)x.im));
    im_tmp = (im*((float)x.re) - re*((float)x.im)) / (((float)x.re)*((float)x.re) + ((float)x.im)*((float)x.im));
    re = re_tmp;
    im = im_tmp;
    return *this;
  }
  __host__ __device__ hr_complex_flt operator/=(const hr_complex_flt x)
  {
    float re_tmp, im_tmp;
    re_tmp = (re*((float)x.re) + im*((float)x.im)) / (((float)x.re)*((float)x.re) + ((float)x.im)*((float)x.im));
    im_tmp = (im*((float)x.re) - re*((float)x.im)) / (((float)x.re)*((float)x.re) + ((float)x.im)*((float)x.im));
    re = re_tmp;
    im = im_tmp;
    return *this;
  }
  __host__ __device__ hr_complex_flt operator/=(const hr_complex x);
  __host__ __device__ hr_complex_flt operator+(const int x)
  {
    return hr_complex_flt(re+(float)x, im);
  }
  __host__ __device__ hr_complex_flt operator+(const float x)
  {
    return hr_complex_flt(re+(float)x, im);
  }
  __host__ __device__ hr_complex_flt operator+(const double x)
  {
    return hr_complex_flt(re+(float)x, im);
  }
  __host__ __device__ hr_complex_flt operator+(const hr_complex_int x)
  {
    return hr_complex_flt(re+(float)x.re, im+(float)x.im);
  }
  __host__ __device__ hr_complex_flt operator+(const hr_complex_flt x)
  {
    return hr_complex_flt(re+(float)x.re, im+(float)x.im);
  }
  __host__ __device__ hr_complex_flt operator-(const int x)
  {
    return hr_complex_flt(re-(float)x, im);
  }
  __host__ __device__ hr_complex_flt operator-(const float x)
  {
    return hr_complex_flt(re-(float)x, im);
  }
  __host__ __device__ hr_complex_flt operator-(const double x)
  {
    return hr_complex_flt(re-(float)x, im);
  }
  __host__ __device__ hr_complex_flt operator-(const hr_complex_int x)
  {
    return hr_complex_flt(re-(float)x.re, im-(float)x.im);
  }
  __host__ __device__ hr_complex_flt operator-(const hr_complex_flt x)
  {
    return hr_complex_flt(re-(float)x.re, im-(float)x.im);
  }
  __host__ __device__ hr_complex_flt operator*(const int x)
  {
    return hr_complex_flt(re*(float)x, im*(float)x);
  }
  __host__ __device__ hr_complex_flt operator*(const float x)
  {
    return hr_complex_flt(re*(float)x, im*(float)x);
  }
  __host__ __device__ hr_complex_flt operator*(const double x)
  {
    return hr_complex_flt(re*(float)x, im*(float)x);
  }
  __host__ __device__ hr_complex_flt operator*(const hr_complex_int x)
  {
    return hr_complex_flt(re*((float)x.re) - im*((float)x.im), re*((float)x.im) + im*((float)x.re));
  }
  __host__ __device__ hr_complex_flt operator*(const hr_complex_flt x)
  {
    return hr_complex_flt(re*((float)x.re) - im*((float)x.im), re*((float)x.im) + im*((float)x.re));
  }
  __host__ __device__ hr_complex_flt operator/(const int x)
  {
    return hr_complex_flt(re/(float)x, im/(float)x);
  }
  __host__ __device__ hr_complex_flt operator/(const float x)
  {
    return hr_complex_flt(re/(float)x, im/(float)x);
  }
  __host__ __device__ hr_complex_flt operator/(const double x)
  {
    return hr_complex_flt(re/(float)x, im/(float)x);
  }
  __host__ __device__ hr_complex_flt operator/(const hr_complex_int x)
  {
    return hr_complex_flt((re*((float)x.re) + im*((float)x.im)) / (((float)x.re)*((float)x.re) + ((float)x.im)*((float)x.im)), (im*((float)x.re) - re*((float)x.im)) / (((float)x.re)*((float)x.re) + ((float)x.im)*((float)x.im)));
  }
  __host__ __device__ hr_complex_flt operator/(const hr_complex_flt x)
  {
    return hr_complex_flt((re*((float)x.re) + im*((float)x.im)) / (((float)x.re)*((float)x.re) + ((float)x.im)*((float)x.im)), (im*((float)x.re) - re*((float)x.im)) / (((float)x.re)*((float)x.re) + ((float)x.im)*((float)x.im)));
  }
};
struct hr_complex{
  double re, im;
  __host__ __device__ hr_complex(void){}
  __host__ __device__ hr_complex(const int x)
  {
    re = (double)x;
    im = (double)0;
  }
  __host__ __device__ hr_complex(const float x)
  {
    re = (double)x;
    im = (double)0;
  }
  __host__ __device__ hr_complex(const double x)
  {
    re = (double)x;
    im = (double)0;
  }
  __host__ __device__ hr_complex(const int a, const int b)
  {
    re = (double)a;
    im = (double)b;
  }
  __host__ __device__ hr_complex(const int a, const float b)
  {
    re = (double)a;
    im = (double)b;
  }
  __host__ __device__ hr_complex(const int a, const double b)
  {
    re = (double)a;
    im = (double)b;
  }
  __host__ __device__ hr_complex(const float a, const int b)
  {
    re = (double)a;
    im = (double)b;
  }
  __host__ __device__ hr_complex(const float a, const float b)
  {
    re = (double)a;
    im = (double)b;
  }
  __host__ __device__ hr_complex(const float a, const double b)
  {
    re = (double)a;
    im = (double)b;
  }
  __host__ __device__ hr_complex(const double a, const int b)
  {
    re = (double)a;
    im = (double)b;
  }
  __host__ __device__ hr_complex(const double a, const float b)
  {
    re = (double)a;
    im = (double)b;
  }
  __host__ __device__ hr_complex(const double a, const double b)
  {
    re = (double)a;
    im = (double)b;
  }
  __host__ __device__ hr_complex(const hr_complex_int x)
  {
    re = (double)x.re;
    im = (double)x.im;
  }
  __host__ __device__ hr_complex(const hr_complex_flt x)
  {
    re = (double)x.re;
    im = (double)x.im;
  }
  __host__ __device__ hr_complex operator=(const int x)
  {
    re = (double)x;
    im = (double)0;
    return *this;
  }
  __host__ __device__ hr_complex operator=(const float x)
  {
    re = (double)x;
    im = (double)0;
    return *this;
  }
  __host__ __device__ hr_complex operator=(const double x)
  {
    re = (double)x;
    im = (double)0;
    return *this;
  }
  __host__ __device__ hr_complex operator=(const hr_complex_int x)
  {
    re = (double)x.re;
    im = (double)x.im;
    return *this;
  }
  __host__ __device__ hr_complex operator=(const hr_complex_flt x)
  {
    re = (double)x.re;
    im = (double)x.im;
    return *this;
  }
  __host__ __device__ hr_complex operator=(const hr_complex x)
  {
    re = (double)x.re;
    im = (double)x.im;
    return *this;
  }
  __host__ __device__ hr_complex operator-(void)
  {
    return hr_complex(-re, -im);
  }
  __host__ __device__ hr_complex operator+=(const int x)
  {
    re += (double)x;
    return *this;
  }
  __host__ __device__ hr_complex operator+=(const float x)
  {
    re += (double)x;
    return *this;
  }
  __host__ __device__ hr_complex operator+=(const double x)
  {
    re += (double)x;
    return *this;
  }
  __host__ __device__ hr_complex operator+=(const hr_complex_int x)
  {
    re += (double)x.re;
    im += (double)x.im;
    return *this;
  }
  __host__ __device__ hr_complex operator+=(const hr_complex_flt x)
  {
    re += (double)x.re;
    im += (double)x.im;
    return *this;
  }
  __host__ __device__ hr_complex operator+=(const hr_complex x)
  {
    re += (double)x.re;
    im += (double)x.im;
    return *this;
  }
  __host__ __device__ hr_complex operator-=(const int x)
  {
    re -= (double)x;
    return *this;
  }
  __host__ __device__ hr_complex operator-=(const float x)
  {
    re -= (double)x;
    return *this;
  }
  __host__ __device__ hr_complex operator-=(const double x)
  {
    re -= (double)x;
    return *this;
  }
  __host__ __device__ hr_complex operator-=(const hr_complex_int x)
  {
    re -= (double)x.re;
    im -= (double)x.im;
    return *this;
  }
  __host__ __device__ hr_complex operator-=(const hr_complex_flt x)
  {
    re -= (double)x.re;
    im -= (double)x.im;
    return *this;
  }
  __host__ __device__ hr_complex operator-=(const hr_complex x)
  {
    re -= (double)x.re;
    im -= (double)x.im;
    return *this;
  }
  __host__ __device__ hr_complex operator*=(const int x)
  {
    re *= (double)x;
    im *= (double)x;
    return *this;
  }
  __host__ __device__ hr_complex operator*=(const float x)
  {
    re *= (double)x;
    im *= (double)x;
    return *this;
  }
  __host__ __device__ hr_complex operator*=(const double x)
  {
    re *= (double)x;
    im *= (double)x;
    return *this;
  }
  __host__ __device__ hr_complex operator*=(const hr_complex_int x)
  {
    double re_tmp, im_tmp;
    re_tmp = re*((double)x.re) - im*((double)x.im);
    im_tmp = re*((double)x.im) + im*((double)x.re);
    re = re_tmp;
    im = im_tmp;
    return *this;
  }
  __host__ __device__ hr_complex operator*=(const hr_complex_flt x)
  {
    double re_tmp, im_tmp;
    re_tmp = re*((double)x.re) - im*((double)x.im);
    im_tmp = re*((double)x.im) + im*((double)x.re);
    re = re_tmp;
    im = im_tmp;
    return *this;
  }
  __host__ __device__ hr_complex operator*=(const hr_complex x)
  {
    double re_tmp, im_tmp;
    re_tmp = re*((double)x.re) - im*((double)x.im);
    im_tmp = re*((double)x.im) + im*((double)x.re);
    re = re_tmp;
    im = im_tmp;
    return *this;
  }
  __host__ __device__ hr_complex operator/=(const int x)
  {
    re /= (double)x;
    im /= (double)x;
    return *this;
  }
  __host__ __device__ hr_complex operator/=(const float x)
  {
    re /= (double)x;
    im /= (double)x;
    return *this;
  }
  __host__ __device__ hr_complex operator/=(const double x)
  {
    re /= (double)x;
    im /= (double)x;
    return *this;
  }
  __host__ __device__ hr_complex operator/=(const hr_complex_int x)
  {
    double re_tmp, im_tmp;
    re_tmp = (re*((double)x.re) + im*((double)x.im)) / (((double)x.re)*((double)x.re) + ((double)x.im)*((double)x.im));
    im_tmp = (im*((double)x.re) - re*((double)x.im)) / (((double)x.re)*((double)x.re) + ((double)x.im)*((double)x.im));
    re = re_tmp;
    im = im_tmp;
    return *this;
  }
  __host__ __device__ hr_complex operator/=(const hr_complex_flt x)
  {
    double re_tmp, im_tmp;
    re_tmp = (re*((double)x.re) + im*((double)x.im)) / (((double)x.re)*((double)x.re) + ((double)x.im)*((double)x.im));
    im_tmp = (im*((double)x.re) - re*((double)x.im)) / (((double)x.re)*((double)x.re) + ((double)x.im)*((double)x.im));
    re = re_tmp;
    im = im_tmp;
    return *this;
  }
  __host__ __device__ hr_complex operator/=(const hr_complex x)
  {
    double re_tmp, im_tmp;
    re_tmp = (re*((double)x.re) + im*((double)x.im)) / (((double)x.re)*((double)x.re) + ((double)x.im)*((double)x.im));
    im_tmp = (im*((double)x.re) - re*((double)x.im)) / (((double)x.re)*((double)x.re) + ((double)x.im)*((double)x.im));
    re = re_tmp;
    im = im_tmp;
    return *this;
  }
  __host__ __device__ hr_complex operator+(const int x)
  {
    return hr_complex(re+(double)x, im);
  }
  __host__ __device__ hr_complex operator+(const float x)
  {
    return hr_complex(re+(double)x, im);
  }
  __host__ __device__ hr_complex operator+(const double x)
  {
    return hr_complex(re+(double)x, im);
  }
  __host__ __device__ hr_complex operator+(const hr_complex_int x)
  {
    return hr_complex(re+(double)x.re, im+(double)x.im);
  }
  __host__ __device__ hr_complex operator+(const hr_complex_flt x)
  {
    return hr_complex(re+(double)x.re, im+(double)x.im);
  }
  __host__ __device__ hr_complex operator+(const hr_complex x)
  {
    return hr_complex(re+(double)x.re, im+(double)x.im);
  }
  __host__ __device__ hr_complex operator-(const int x)
  {
    return hr_complex(re-(double)x, im);
  }
  __host__ __device__ hr_complex operator-(const float x)
  {
    return hr_complex(re-(double)x, im);
  }
  __host__ __device__ hr_complex operator-(const double x)
  {
    return hr_complex(re-(double)x, im);
  }
  __host__ __device__ hr_complex operator-(const hr_complex_int x)
  {
    return hr_complex(re-(double)x.re, im-(double)x.im);
  }
  __host__ __device__ hr_complex operator-(const hr_complex_flt x)
  {
    return hr_complex(re-(double)x.re, im-(double)x.im);
  }
  __host__ __device__ hr_complex operator-(const hr_complex x)
  {
    return hr_complex(re-(double)x.re, im-(double)x.im);
  }
  __host__ __device__ hr_complex operator*(const int x)
  {
    return hr_complex(re*(double)x, im*(double)x);
  }
  __host__ __device__ hr_complex operator*(const float x)
  {
    return hr_complex(re*(double)x, im*(double)x);
  }
  __host__ __device__ hr_complex operator*(const double x)
  {
    return hr_complex(re*(double)x, im*(double)x);
  }
  __host__ __device__ hr_complex operator*(const hr_complex_int x)
  {
    return hr_complex(re*((double)x.re) - im*((double)x.im), re*((double)x.im) + im*((double)x.re));
  }
  __host__ __device__ hr_complex operator*(const hr_complex_flt x)
  {
    return hr_complex(re*((double)x.re) - im*((double)x.im), re*((double)x.im) + im*((double)x.re));
  }
  __host__ __device__ hr_complex operator*(const hr_complex x)
  {
    return hr_complex(re*((double)x.re) - im*((double)x.im), re*((double)x.im) + im*((double)x.re));
  }
  __host__ __device__ hr_complex operator/(const int x)
  {
    return hr_complex(re/(double)x, im/(double)x);
  }
  __host__ __device__ hr_complex operator/(const float x)
  {
    return hr_complex(re/(double)x, im/(double)x);
  }
  __host__ __device__ hr_complex operator/(const double x)
  {
    return hr_complex(re/(double)x, im/(double)x);
  }
  __host__ __device__ hr_complex operator/(const hr_complex_int x)
  {
    return hr_complex((re*((double)x.re) + im*((double)x.im)) / (((double)x.re)*((double)x.re) + ((double)x.im)*((double)x.im)), (im*((double)x.re) - re*((double)x.im)) / (((double)x.re)*((double)x.re) + ((double)x.im)*((double)x.im)));
  }
  __host__ __device__ hr_complex operator/(const hr_complex_flt x)
  {
    return hr_complex((re*((double)x.re) + im*((double)x.im)) / (((double)x.re)*((double)x.re) + ((double)x.im)*((double)x.im)), (im*((double)x.re) - re*((double)x.im)) / (((double)x.re)*((double)x.re) + ((double)x.im)*((double)x.im)));
  }
  __host__ __device__ hr_complex operator/(const hr_complex x)
  {
    return hr_complex((re*((double)x.re) + im*((double)x.im)) / (((double)x.re)*((double)x.re) + ((double)x.im)*((double)x.im)), (im*((double)x.re) - re*((double)x.im)) / (((double)x.re)*((double)x.re) + ((double)x.im)*((double)x.im)));
  }
};
__host__ __device__ hr_complex_int operator+(const int x, hr_complex_int c){return c+x;}
__host__ __device__ hr_complex_int operator+(const float x, hr_complex_int c){return c+x;}
__host__ __device__ hr_complex_int operator+(const double x, hr_complex_int c){return c+x;}
__host__ __device__ hr_complex_int operator-(const int x, hr_complex_int c){return -c+x;}
__host__ __device__ hr_complex_int operator-(const float x, hr_complex_int c){return -c+x;}
__host__ __device__ hr_complex_int operator-(const double x, hr_complex_int c){return -c+x;}
__host__ __device__ hr_complex_int operator*(const int x, hr_complex_int c){return c*x;}
__host__ __device__ hr_complex_int operator*(const float x, hr_complex_int c){return c*x;}
__host__ __device__ hr_complex_int operator*(const double x, hr_complex_int c){return c*x;}
__host__ __device__ hr_complex_int operator/(const int x, hr_complex_int c){return ((hr_complex_int)x)/c;}
__host__ __device__ hr_complex_int operator/(const float x, hr_complex_int c){return ((hr_complex_int)x)/c;}
__host__ __device__ hr_complex_int operator/(const double x, hr_complex_int c){return ((hr_complex_int)x)/c;}
__host__ __device__ hr_complex_flt operator+(const int x, hr_complex_flt c){return c+x;}
__host__ __device__ hr_complex_flt operator+(const float x, hr_complex_flt c){return c+x;}
__host__ __device__ hr_complex_flt operator+(const double x, hr_complex_flt c){return c+x;}
__host__ __device__ hr_complex_flt operator+(const hr_complex_int x, hr_complex_flt c){return c+x;}
__host__ __device__ hr_complex_flt operator-(const int x, hr_complex_flt c){return -c+x;}
__host__ __device__ hr_complex_flt operator-(const float x, hr_complex_flt c){return -c+x;}
__host__ __device__ hr_complex_flt operator-(const double x, hr_complex_flt c){return -c+x;}
__host__ __device__ hr_complex_flt operator-(const hr_complex_int x, hr_complex_flt c){return -c+x;}
__host__ __device__ hr_complex_flt operator*(const int x, hr_complex_flt c){return c*x;}
__host__ __device__ hr_complex_flt operator*(const float x, hr_complex_flt c){return c*x;}
__host__ __device__ hr_complex_flt operator*(const double x, hr_complex_flt c){return c*x;}
__host__ __device__ hr_complex_flt operator*(const hr_complex_int x, hr_complex_flt c){return c*x;}
__host__ __device__ hr_complex_flt operator/(const int x, hr_complex_flt c){return ((hr_complex_flt)x)/c;}
__host__ __device__ hr_complex_flt operator/(const float x, hr_complex_flt c){return ((hr_complex_flt)x)/c;}
__host__ __device__ hr_complex_flt operator/(const double x, hr_complex_flt c){return ((hr_complex_flt)x)/c;}
__host__ __device__ hr_complex_flt operator/(const hr_complex_int x, hr_complex_flt c){return ((hr_complex_flt)x)/c;}
__host__ __device__ hr_complex operator+(const int x, hr_complex c){return c+x;}
__host__ __device__ hr_complex operator+(const float x, hr_complex c){return c+x;}
__host__ __device__ hr_complex operator+(const double x, hr_complex c){return c+x;}
__host__ __device__ hr_complex operator+(const hr_complex_int x, hr_complex c){return c+x;}
__host__ __device__ hr_complex operator+(const hr_complex_flt x, hr_complex c){return c+x;}
__host__ __device__ hr_complex operator-(const int x, hr_complex c){return -c+x;}
__host__ __device__ hr_complex operator-(const float x, hr_complex c){return -c+x;}
__host__ __device__ hr_complex operator-(const double x, hr_complex c){return -c+x;}
__host__ __device__ hr_complex operator-(const hr_complex_int x, hr_complex c){return -c+x;}
__host__ __device__ hr_complex operator-(const hr_complex_flt x, hr_complex c){return -c+x;}
__host__ __device__ hr_complex operator*(const int x, hr_complex c){return c*x;}
__host__ __device__ hr_complex operator*(const float x, hr_complex c){return c*x;}
__host__ __device__ hr_complex operator*(const double x, hr_complex c){return c*x;}
__host__ __device__ hr_complex operator*(const hr_complex_int x, hr_complex c){return c*x;}
__host__ __device__ hr_complex operator*(const hr_complex_flt x, hr_complex c){return c*x;}
__host__ __device__ hr_complex operator/(const int x, hr_complex c){return ((hr_complex)x)/c;}
__host__ __device__ hr_complex operator/(const float x, hr_complex c){return ((hr_complex)x)/c;}
__host__ __device__ hr_complex operator/(const double x, hr_complex c){return ((hr_complex)x)/c;}
__host__ __device__ hr_complex operator/(const hr_complex_int x, hr_complex c){return ((hr_complex)x)/c;}
__host__ __device__ hr_complex operator/(const hr_complex_flt x, hr_complex c){return ((hr_complex)x)/c;}
__host__ __device__  hr_complex_int::hr_complex_int(const hr_complex_flt x)
{
  re = (int)x.re;
  im = (int)x.im;
}
__host__ __device__  hr_complex_int::hr_complex_int(const hr_complex x)
{
  re = (int)x.re;
  im = (int)x.im;
}
__host__ __device__ hr_complex_int hr_complex_int::operator=(const hr_complex_flt x)
{
  re = (int)x.re;
  im = (int)x.im;
  return *this;
}
__host__ __device__ hr_complex_int hr_complex_int::operator=(const hr_complex x)
{
  re = (int)x.re;
  im = (int)x.im;
  return *this;
}
__host__ __device__ hr_complex_int hr_complex_int::operator+=(const hr_complex_flt x)
{
  re += (int)x.re;
  im += (int)x.im;
  return *this;
}
__host__ __device__ hr_complex_int hr_complex_int::operator+=(const hr_complex x)
{
  re += (int)x.re;
  im += (int)x.im;
  return *this;
}
__host__ __device__ hr_complex_int hr_complex_int::operator-=(const hr_complex_flt x)
{
  re -= (int)x.re;
  im -= (int)x.im;
  return *this;
}
__host__ __device__ hr_complex_int hr_complex_int::operator-=(const hr_complex x)
{
  re -= (int)x.re;
  im -= (int)x.im;
  return *this;
}
__host__ __device__ hr_complex_int hr_complex_int::operator*=(const hr_complex_flt x)
{
  int re_tmp, im_tmp;
  re_tmp = re*((int)x.re) - im*((int)x.im);
  im_tmp = re*((int)x.im) + im*((int)x.re);
  re = re_tmp;
  im = im_tmp;
  return *this;
}
__host__ __device__ hr_complex_int hr_complex_int::operator*=(const hr_complex x)
{
  int re_tmp, im_tmp;
  re_tmp = re*((int)x.re) - im*((int)x.im);
  im_tmp = re*((int)x.im) + im*((int)x.re);
  re = re_tmp;
  im = im_tmp;
  return *this;
}
__host__ __device__ hr_complex_int hr_complex_int::operator/=(const hr_complex_flt x)
{
  int re_tmp, im_tmp;
  re_tmp = (re*((int)x.re) + im*((int)x.im)) / (((int)x.re)*((int)x.re) + ((int)x.im)*((int)x.im));
  im_tmp = (im*((int)x.re) - re*((int)x.im)) / (((int)x.re)*((int)x.re) + ((int)x.im)*((int)x.im));
  re = re_tmp;
  im = im_tmp;
  return *this;
}
__host__ __device__ hr_complex_int hr_complex_int::operator/=(const hr_complex x)
{
  int re_tmp, im_tmp;
  re_tmp = (re*((int)x.re) + im*((int)x.im)) / (((int)x.re)*((int)x.re) + ((int)x.im)*((int)x.im));
  im_tmp = (im*((int)x.re) - re*((int)x.im)) / (((int)x.re)*((int)x.re) + ((int)x.im)*((int)x.im));
  re = re_tmp;
  im = im_tmp;
  return *this;
}
__host__ __device__  hr_complex_flt::hr_complex_flt(const hr_complex x)
{
  re = (float)x.re;
  im = (float)x.im;
}
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
__host__ __device__ hr_complex_flt hr_complex_flt::operator-=(const hr_complex x)
{
  re -= (float)x.re;
  im -= (float)x.im;
  return *this;
}
__host__ __device__ hr_complex_flt hr_complex_flt::operator*=(const hr_complex x)
{
  float re_tmp, im_tmp;
  re_tmp = re*((float)x.re) - im*((float)x.im);
  im_tmp = re*((float)x.im) + im*((float)x.re);
  re = re_tmp;
  im = im_tmp;
  return *this;
}
__host__ __device__ hr_complex_flt hr_complex_flt::operator/=(const hr_complex x)
{
  float re_tmp, im_tmp;
  re_tmp = (re*((float)x.re) + im*((float)x.im)) / (((float)x.re)*((float)x.re) + ((float)x.im)*((float)x.im));
  im_tmp = (im*((float)x.re) - re*((float)x.im)) / (((float)x.re)*((float)x.re) + ((float)x.im)*((float)x.im));
  re = re_tmp;
  im = im_tmp;
  return *this;
}
/*******************************************************************************
*
* End of generated code
*
*******************************************************************************/
#define I (hr_complex_int(0, 1))
#define creal(a) ((a).re)
#define cimag(a) ((a).im)
#define conj(a) (hr_complex((a).re, -(a).im))

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
   creal(conj(a)*(b))

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
