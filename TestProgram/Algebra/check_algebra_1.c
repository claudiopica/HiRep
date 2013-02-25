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

#ifdef REPR_FUNDAMENTAL
static float C2=(float)(NG*NG-1)/(float)(2*NG);
static float Tr=0.5;
#endif


#ifdef REPR_ADJOINT
static float C2=(float)NG;
static float Tr=(float)NG;
#endif

#ifdef REPR_ANTISYMMETRIC
static float C2=(float)(NG-2)*(NG+1)/(float)NG;
static float Tr=(float)(NG-2)/2;
#endif

#ifdef REPR_SYMMETRIC
static float C2=(float)(NG+2)*(NG-1)/(float)NG;
static float Tr=(float)(NG+2)/2;
#endif

#ifdef WITH_MPI
#error: check_algebra_1 only works only on serial jobs
#endif


static int dAdj=NG*NG-1;
static float fund=(float)(NG*NG-1)/(2*(float)(NG));

int main(int argc,char *argv[])
{
   suNg_algebra_vector f[dAdj];
   suNg A,B,TMP,CAS;
   suNf a,b,tmp,cas;
   double tau,trace;
   int i,j;
   
   
   printf("Gauge group: SU(%d)\n",NG);
   printf("Fermion representation: dim = %d\n",NF);
   printf("\n");

   for (i=0;i<dAdj;i++)
   {
      _algebra_vector_zero_g(f[i]);
      f[i].c[i]=1.;
   }

   for (i=0;i<dAdj;i++)
   {
      for (j=0;j<dAdj;j++)
      {
         _algebra_represent(a,f[i]);
         _algebra_represent(b,f[j]);
         
         _suNf_times_suNf(tmp,a,b);
         _suNf_trace_re(trace,tmp);
         printf("tr_R (T[%d] T[%d]): %.4f ", i, j,trace);
         if (i==j)
            printf("  [should be: %.4f]\n",-Tr);
         else
            printf("  [should be: 0.00]\n");

         _fund_algebra_represent(A,f[i]);
         _fund_algebra_represent(B,f[j]);
         
         _suNg_times_suNg(TMP,A,B);
         _suNg_trace_re(trace,TMP);
         printf("tr_f (T[%d] T[%d]): %.4f ", i, j,trace);
         if (i==j)
            printf("  [should be: %.4f]\n",-0.5);
         else
            printf("  [should be: 0.00]\n");
      }
   }
   
	 _algebra_represent(a,f[0]);
	 _fund_algebra_represent(A,f[0]);
	 _suNf_times_suNf(cas,a,a);
	 _suNg_times_suNg(CAS,A,A);
         
   for (i=1;i<dAdj;i++)
   {
      _algebra_represent(a,f[i]);
      _fund_algebra_represent(A,f[i]);

      _suNf_times_suNf(tmp,a,a);
      _suNf_add_assign(cas,tmp);
      _suNg_times_suNg(TMP,A,A);
      _suNg_add_assign(CAS,TMP);
   }

   _suNf_unit(tmp);
   _suNf_mul(tmp,C2,tmp);
   _suNf_add_assign(cas,tmp);
   _suNf_sqnorm(tau,cas);
   printf("casimir check: %.3f\n",tau);
   printf("(should be 0.00)\n");

   _suNg_unit(TMP);
   _suNg_mul(TMP,fund,TMP);
   _suNg_add_assign(CAS,TMP);
   _suNg_sqnorm(tau,CAS);
   printf("casimir check: %.3f\n",tau);
   printf("(should be 0.00)\n");

   exit(0);
}
