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

/*
void assign_u2ud(void)
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

      for (i=0;i<(NG*NG);++i)
      {
        r[i].re=(double)(rf[i].re);
        r[i].im=(double)(rf[i].im);
      }

      project_to_suNg(pu_gauge(ix,mu));
    }
  }
}


void assign_ud2u(void)
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

      for (i=0;i<(NG*NG);++i)
      {
        rf[i].re=(float)(r[i].re);
        rf[i].im=(float)(r[i].im);
      }
    }
  }
}
*/

void assign_ud2u_f(void)
{
  if(u_gauge_f_flt!=NULL){
    double *d;
    float *f;
  
    d = (double*)(u_gauge_f->ptr);
    f = (float*)(u_gauge_f_flt->ptr);

_OMP_PRAGMA( _omp_parallel )
_OMP_PRAGMA( _omp_for )
    for(int i=0; i<4*glattice.gsize_gauge*(sizeof(suNf)/sizeof(double)); i++) {
      *(f+i) = (float)(*(d+i));
    }
  }
}

void assign_s2sd(spinor_field *out, spinor_field_flt *in) {

  _TWO_SPINORS_FOR(out,in) {
    double *o=(double*)_SPINOR_PTR(out);
    float *i=(float*)_SPINOR_PTR(in);
    for (int n=0;n<(8*NF);n++) { *(o++) = (double) *(i++); }
  }

}

void assign_sd2s(spinor_field_flt *out, spinor_field *in) {

  _TWO_SPINORS_FOR(out,in) {
    float *o=(float*)_SPINOR_PTR(out);
    double *i=(double*)_SPINOR_PTR(in);
    for (int n=0;n<(8*NF);n++) { *(o++) = (float) *(i++); }
  }

}
