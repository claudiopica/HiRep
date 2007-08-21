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

#ifndef REPR_ADJOINT
#error : check8 is written for REPR_ADJOINT
#endif

static int dAdj=NG*NG-1;

int main(int argc,char *argv[])
{
   suNg_algebra_vector f[dAdj];
   suNg A,B,TMP,CAS;
   suNf a,b,tmp,cas;
   double tau,fund,sym,trace;
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
         printf("tr_R (T[%d] T[%d]): %.2f ", i, j,trace);
         if (i==j)
            printf("  [should be: 2.50]\n");
         else
            printf("  [should be: 0.00]\n");

         _fund_algebra_represent(A,f[i]);
         _fund_algebra_represent(B,f[j]);
         
         _suNg_times_suNg(TMP,A,B);
				 _suNg_trace_re(trace,TMP);
         printf("tr_f (T[%d] T[%d]): %.2f ", i, j,trace);
         if (i==j)
            printf("  [should be: 0.50]\n");
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

   sym=-4.0*((double)NG-1.0)*((double)NG+2.0)/((double)NG);
   fund=-4.0*((double)NG*(double)NG-1.0)/(2.0*(double)NG);
      
   _suNf_unit(tmp);
   _suNf_mul(tmp,sym,tmp);
   _suNf_sub_assign(cas,tmp);
   _suNf_sqnorm(tau,cas);
   printf("casimir check: %.3f\n",tau);
   printf("(should be 0.00)\n");

   _suNg_unit(TMP);
   _suNg_mul(TMP,fund,TMP);
   _suNg_sub_assign(CAS,TMP);
   _suNg_sqnorm(tau,CAS);
   printf("casimir check: %.3f\n",tau);
   printf("(should be 0.00)\n");


   exit(0);
}
