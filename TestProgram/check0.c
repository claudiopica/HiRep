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


int main(int argc,char *argv[])
{
   suNg A,B,C,E;
   suNf a,b,c,e,tmp;
   int level,seed;
   double tau;
   
   printf("Gauge group: SU(%d)\n",NG);
   printf("Fermion representation: dim = %d\n",NF);
   printf("\n");
   
   level=0;
   seed=123;
   rlxs_init(level,seed);
   printf("ranlux: level = %d, seed = %d\n\n",level,seed); 
   fflush(stdout);
 
   _suNg_unit(E);
   _suNf_unit(e);
   
   _group_represent(tmp,E);
   _suNf_sub_assign(e,tmp);

   tau=_suNf_sqnorm(e);
   printf("checking that _group_represent works on E: %.3f\n",tau);
   printf("(should be 0.00)\n");
   

   printf("Generating random matrices A and B... ");
   fflush(stdout);
   random_suNg(&A);
	 random_suNg(&B);
   printf("done.\n");

   _suNg_times_suNg(C,A,B);

   _group_represent(a,A);
   _group_represent(b,B);

   _suNf_times_suNf(c,a,b);
   _group_represent(tmp,C);

   _suNf_sub_assign(c,tmp);

   tau=_suNf_sqnorm(c);
   printf("checking that _group_represent is a homo: %.3f\n",tau);
   printf("(should be 0.00)\n");
   
   exit(0);
}
