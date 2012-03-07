/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *   
* All rights reserved.                                                      * 
\***************************************************************************/

/*******************************************************************************
*
* File single_double_utils.c
*
* Functions for conversion from single to double precision and viceversa
*
*******************************************************************************/

#include <stdlib.h>
#include "utils.h"
#include "suN.h"
#include "error.h"
#include "global.h"
#include "spinor_field.h"

void assign_u2ud_cpu(void)
{
  _DECLARE_INT_ITERATOR(ix);
  int i,mu;
  complex *r;
  complex_flt *rf;

  _MASTER_FOR(&glattice,ix){
    for (mu=0;mu<4;mu++)
    {
      r=(complex*)(pu_gauge(ix,mu));
      rf=(complex_flt*)(pu_gauge_flt(ix,mu));
      
#ifdef WITH_QUATERNIONS
#if (NG != 2)
#error "WITH_QUATERNIONS not defined with NG!=2"
#else
      r[0].re=(double)(rf[0].re);
      r[0].im=(double)(rf[0].im);
      r[1].re=(double)(rf[1].re);
      r[1].im=(double)(rf[1].im);
#endif
#else
      for (i=0;i<(NG*NG);++i)
      {
        r[i].re=(double)(rf[i].re);
        r[i].im=(double)(rf[i].im);
      }
#endif
      
      project_to_suNg(pu_gauge(ix,mu));
    }
  }
}


void assign_ud2u_cpu(void)
{
  _DECLARE_INT_ITERATOR(ix);
  int i,mu;
  complex *r;
  complex_flt *rf;

  _MASTER_FOR(&glattice,ix){
    for (mu=0;mu<4;mu++)
    {
      r=(complex*)(pu_gauge(ix,mu));
      rf=(complex_flt*)(pu_gauge_flt(ix,mu));

#ifdef WITH_QUATERNIONS
#if (NG != 2)
#error "WITH_QUATERNIONS not defined with NG!=2"
#else
      rf[0].re=(float)(r[0].re);
      rf[0].im=(float)(r[0].im);
      rf[1].re=(float)(r[1].re);
      rf[1].im=(float)(r[1].im);
#endif
#else
      for (i=0;i<(NG*NG);++i)
      {
        rf[i].re=(float)(r[i].re);
        rf[i].im=(float)(r[i].im);
      }
#endif
      
    }
  }
}

void assign_ud2u_f_cpu(void)
{
#ifdef WITH_QUATERNIONS
  assign_ud2u_cpu();
#else
  _DECLARE_INT_ITERATOR(ix);
  int i,mu;

#ifdef REPR_ADJOINT
  double *r;
  float *rf;
#else
  complex *r;
  complex_flt *rf;
#endif

  _MASTER_FOR(&glattice,ix){
    for (mu=0;mu<4;mu++)
    {
#ifdef REPR_ADJOINT
      r=(double*)(pu_gauge_f(ix,mu));
      rf=(float*)(pu_gauge_f_flt(ix,mu));

      for (i=0;i<(NF*NF);++i)
        rf[i]=(float)(r[i]);
#else
      r=(complex*)(pu_gauge_f(ix,mu));
      rf=(complex_flt*)(pu_gauge_f_flt(ix,mu));

      for (i=0;i<(NF*NF);++i)
      {
        rf[i].re=(float)(r[i].re);
        rf[i].im=(float)(r[i].im);
      }
#endif
    }
  }
#endif //WITH_QUATERNIONS
}

void assign_u2ud_f_cpu(void)
{
#ifdef WITH_QUATERNIONS
  assign_u2ud_cpu();
#else
  _DECLARE_INT_ITERATOR(ix);
  int i,mu;
#ifdef REPR_ADJOINT
  double *r;
  float *rf;
#else
  complex *r;
  complex_flt *rf;
#endif
  
  _MASTER_FOR(&glattice,ix){
    for (mu=0;mu<4;mu++)
    {
#ifdef REPR_ADJOINT
      r=(double*)(pu_gauge_f(ix,mu));
      rf=(float*)(pu_gauge_f_flt(ix,mu));
      
      for (i=0;i<(NF*NF);++i)
        r[i]=(double)(rf[i]);
#else
      r=(complex*)(pu_gauge_f(ix,mu));
      rf=(complex_flt*)(pu_gauge_f_flt(ix,mu));
      
      for (i=0;i<(NF*NF);++i)
      {
        r[i].re=(double)(rf[i].re);
        r[i].im=(double)(rf[i].im);
      }
#endif
    }
  }
#endif //WITH_QUATERNIONS
}


void assign_s2sd_cpu(spinor_field *out, spinor_field_flt *in) {

  int n;
  float *i;
  double *o;

  _DECLARE_INT_ITERATOR(ix);
  _DECLARE_SPINOR_ITERATOR(out);
  _DECLARE_SPINOR_FLT_ITERATOR(in);

  _TWO_SPINORS_MATCHING(in,out);

  _TWO_SPINORS_FOR(out,in,ix) {
    o=(double*)_SPINOR_PTR(out);
    i=(float*)_SPINOR_PTR(in);
    for (n=0;n<(8*NF);n++)
    {
      *(o++) = (double) *(i++);
    }
  }

}

void assign_sd2s_cpu(spinor_field_flt *out, spinor_field *in) {

  int n;
  float *o;
  double *i;

  _DECLARE_INT_ITERATOR(ix); 
  _DECLARE_SPINOR_FLT_ITERATOR(out);
  _DECLARE_SPINOR_ITERATOR(in);

  _TWO_SPINORS_MATCHING(in,out);

  _TWO_SPINORS_FOR(out,in,ix) {
    o=(float*)_SPINOR_PTR(out);
    i=(double*)_SPINOR_PTR(in);
    for (n=0;n<(8*NF);n++)
    {
      *(o++) = (float) *(i++);
    }
  }
}

#ifndef WITH_GPU
void (*assign_s2sd)(spinor_field *out, spinor_field_flt *in)=assign_s2sd_cpu;
void (*assign_sd2s)(spinor_field_flt *out, spinor_field *in)=assign_sd2s_cpu;
#endif
/* Gauge fields not yet on GPU */
void (*assign_u2ud) (void)=assign_u2ud_cpu;
void (*assign_ud2u) (void)=assign_ud2u_cpu;
void (*assign_u2ud_f)(void)=assign_u2ud_f_cpu;
void (*assign_ud2u_f) (void)=assign_ud2u_f_cpu;
