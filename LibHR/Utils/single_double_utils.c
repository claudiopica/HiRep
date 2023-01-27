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

#include "utils.h"
#include "libhr_core.h"
#include "memory.h"

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
  #ifdef WITH_GPU
    copy_from_gpu_gfield_f(u_gauge_f);
  #endif

  if (u_gauge_f_flt != NULL)
  {
    double *d;
    float *f;

    d = (double *)(u_gauge_f->ptr);
    f = (float *)(u_gauge_f_flt->ptr);

    _OMP_PRAGMA(_omp_parallel)
    _OMP_PRAGMA(_omp_for)
    for (int i = 0; i < 4 * glattice.gsize_gauge * (sizeof(suNf) / sizeof(double)); i++)
    {
      *(f + i) = (float)(*(d + i));
    }
  }

  #ifdef WITH_GPU
    copy_to_gpu_gfield_f_flt(u_gauge_f_flt);
    start_sendrecv_gfield_f_flt(u_gauge_f_flt);
    complete_sendrecv_gfield_f_flt(u_gauge_f_flt);
  #endif
}

void assign_s2sd(spinor_field *out, spinor_field_flt *in)
{
  #ifdef WITH_GPU
    copy_from_gpu_spinor_field_f_flt(in);
  #endif

  _TWO_SPINORS_FOR(out, in)
  {
    double *o = (double *)_SPINOR_PTR(out);
    float *i = (float *)_SPINOR_PTR(in);
    for (int n = 0; n < (8 * NF); n++)
    {
      *(o++) = (double)*(i++);
    }
  }

  #ifdef WITH_GPU
    copy_to_gpu_spinor_field_f(out);
    start_sendrecv_spinor_field_f(out);
    complete_sendrecv_spinor_field_f(out);
  #endif
}

void assign_sd2s(spinor_field_flt *out, spinor_field *in)
{
  #ifdef WITH_GPU
    copy_from_gpu_spinor_field_f(in);
  #endif

  _TWO_SPINORS_FOR(out, in)
  {
    float *o = (float *)_SPINOR_PTR(out);
    double *i = (double *)_SPINOR_PTR(in);
    for (int n = 0; n < (8 * NF); n++)
    {
      *(o++) = (float)*(i++);
    }
  }

  #ifdef WITH_GPU
    copy_to_gpu_spinor_field_f_flt(out);
    start_sendrecv_spinor_field_f_flt(out);
    complete_sendrecv_spinor_field_f_flt(out);
  #endif
}
