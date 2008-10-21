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

void assign_u2ud(void)
{
   int i,ix,mu;
   complex *r;
   complex_flt *rf;

   for (ix=0;ix<VOLUME;ix++)
   {
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
   int i,ix,mu;
   complex *r;
   complex_flt *rf;

   for (ix=0;ix<VOLUME;ix++)
   {
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


void assign_s2sd(int len, suNf_spinor *out, suNf_spinor_flt *in) {
  int ix;
  float *i;
  double *o;
  
  ix=8*NF*len;
  o=(double*)(out)+(ix-1);
  i=(float*)(in)+(ix-1);

  for (;ix>0;--ix) {
    *o = (double) *i;
    --o;
    --i;
  }
  
}

void assign_sd2s(int len, suNf_spinor_flt *out, suNf_spinor *in) {
  int ix;
  float *o;
  double *i;
  
  ix=8*NF*len;
  o=(float*)(out);
  i=(double*)(in);

  for (;ix>0;--ix) {
    *o = (float) *i;
    ++o;
    ++i;
  }
  
}
