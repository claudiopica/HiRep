/*******************************************************************************
*
* Test of modules
*
*******************************************************************************/

#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "io.h"
#include "random.h"
#include "error.h"
#include "global.h"
#include "suN.h"
#include "suN_types.h"
#include "linear_algebra.h"
#include "inverters.h"
#include "representation.h"
#include "utils.h"
#include "update.h"

#if (NG!=3)
#error: check_algebra_2 only works for SU(3)
#endif

#ifndef REPR_ANTISYMMETRIC
#error: check_algebra_2 only works for 2AS representation
#endif


#ifdef WITH_MPI
#error: check_algebra_2 only works only on serial jobs
#endif

int main(int argc,char *argv[])
{
   suNg A;
   suNf a,s,tmp;
   int i,j;
   complex zp,zm;
   
   printf("Gauge group: SU(%d)\n",NG);
   printf("Fermion representation: dim = %d\n",NF);
   printf("Check 2AS = fund* for SU(3)\n");
   printf("\n");

   random_suNg(&A);
   _group_represent2(&a,&A);

   _complex_1(zp);
   _complex_mulr(zm,-1.0,zp);

   _suNf_zero(s);
   (s).c[2]=zp;
   (s).c[4]=zm;
   (s).c[6]=zp;
   _suNf_times_suNf(tmp,a,s);
   _suNf_times_suNf(a,s,tmp);

   printf("fundamental representation:\n");
   for (i=0;i<NF;i++)
   {
      for (j=0;j<NF;j++)
      {
         printf("%.4f + i %.4f   ",A.c[i*NF+j].re,A.c[i*NF+j].im);
      }
      printf("\n");
   }
   
   printf("2AS representation:\n");
   for (i=0;i<NF;i++)
   {
      for (j=0;j<NF;j++)
      {
         printf("%.4f + i %.4f   ",a.c[i*NF+j].re,a.c[i*NF+j].im);
      }
      printf("\n");
   }
   
   exit(0);
}
