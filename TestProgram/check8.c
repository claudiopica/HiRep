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

#if NG!=3
#error : check8 is written for NG=3
#endif
#if NF!=8
#error : check8 is written for REPR_ADJOINT
#endif

static int dAdj=NG*NG-1;

int main(int argc,char *argv[])
{
   suNg_algebra_vector f[dAdj];
   suNg A,B,TMP,CAS;
   suNf a,b,tmp,cas;
   float tau,fund,sym;
   int i,j;
   
   
   printf("Gauge group: SU(%d)\n",NG);
   printf("Fermion representation: dim = %d\n",NF);
   printf("\n");

   for (i=0;i<dAdj;i++)
   {
      _algebra_vector_zero_g(f[i]);
   }

   f[0].c1=1.0;
   f[1].c2=1.0;
   f[2].c3=1.0;
   f[3].c4=1.0;
   f[4].c5=1.0;
   f[5].c6=1.0;
   f[6].c7=1.0;
   f[7].c8=1.0;

   for (i=0;i<dAdj;i++)
   {
      for (j=0;j<dAdj;j++)
      {
         _algebra_represent(a,f[i]);
         _algebra_represent(b,f[j]);
         
         _suNf_times_suNf(tmp,a,b);
         printf("tr_R (T[%d] T[%d]): %.2f ", i, j,
                _suNf_trace_re(tmp));
         if (i==j)
            printf("  [should be: 2.50]\n");
         else
            printf("  [should be: 0.00]\n");

         _fund_algebra_represent(A,f[i]);
         _fund_algebra_represent(B,f[j]);
         
         _suNg_times_suNg(TMP,A,B);
         printf("tr_f (T[%d] T[%d]): %.2f ", i, j,
                _suNg_trace_re(TMP));
         if (i==j)
            printf("  [should be: 0.50]\n");
         else
            printf("  [should be: 0.00]\n");
      }
   }
   
   for (i=0;i<dAdj;i++)
   {
      _algebra_represent(a,f[i]);
      _fund_algebra_represent(A,f[i]);

      if (i==0)
      {
         _suNf_times_suNf(cas,a,a);
         _suNg_times_suNg(CAS,A,A);
      }
      else
      {
         _suNf_times_suNf(tmp,a,a);
         _suNf_add_assign(cas,tmp);
         _suNg_times_suNg(TMP,A,A);
         _suNg_add_assign(CAS,TMP);
      }
   }

   sym=-4.0*((float)NG-1.0)*((float)NG+2.0)/((float)NG);
   fund=-4.0*((float)NG*(float)NG-1.0)/(2.0*(float)NG);
      
   _suNf_unit(tmp);
   _suNf_mul(tmp,sym,tmp);
   _suNf_sub_assign(cas,tmp);
   tau=_suNf_sqnorm(cas);
   printf("casimir check: %.3f\n",tau);
   printf("(should be 0.00)\n");

   _suNg_unit(TMP);
   _suNg_mul(TMP,fund,TMP);
   _suNg_sub_assign(CAS,TMP);
   tau=_suNg_sqnorm(CAS);
   printf("casimir check: %.3f\n",tau);
   printf("(should be 0.00)\n");


   exit(0);
}
