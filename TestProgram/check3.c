/******************************************************************************
*
* Test of modules
*
******************************************************************************/

#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "io.h"
#include "random.h"
#include "error.h"
#include "geometry.h"
#include "memory.h"
#include "statistics.h"
#include "update.h"
#include "global.h"
#include "observables.h"
#include "suN.h"
#include "suN_types.h"
#include "dirac.h"
#include "linear_algebra.h"
#include "inverters.h"
#include "representation.h"
#include "utils.h"

int nhb,nor,nit,nth,nms,level,seed;
double beta;

static double hmass=0.1;


int main(int argc,char *argv[])
{
   double tau;
   spinor_field *s1,*s2,*s3,*s4;

   printf("Gauge group: SU(%d)\n",NG);
   printf("Fermion representation: dim = %d\n",NF);
   geometry_eo_lexi();
   test_geometry();
   printf("The lattice size is %dx%dx%dx%d\n",T,X,Y,Z);
   printf("\n");
   
   level=1;
   seed=123;
   printf("ranlux: level = %d, seed = %d\n\n",level,seed); 
   fflush(stdout);
   
   rlxd_init(level,seed);

   u_gauge=alloc_gfield();
#ifndef REPR_FUNDAMENTAL
   u_gauge_f=alloc_gfield_f();
#endif
   represent_gauge_field();
	
   printf("Generating a random gauge field... ");fflush(stdout);
   random_u();
   printf("done.\n");
   represent_gauge_field();

   set_spinor_len(VOLUME);
   s1=alloc_spinor_field_f(4);
	 s2=s1+1;
	 s3=s2+1;
	 s4=s3+1;

   gaussian_spinor_field(&(s1[0]));
   gaussian_spinor_field(&(s2[0]));
   
   tau = 1./sqrt(spinor_field_sqnorm_f(s1));
   spinor_field_mul_f(s1,tau,s1);
   tau = 1./sqrt(spinor_field_sqnorm_f(s2));
   spinor_field_mul_f(s2,tau,s2);

   
   printf("Test new Dirac implementation: ");
   
   g5Dphi_old(hmass,_SPINOR_ADDR(s3),_SPINOR_ADDR(s1));
   g5Dphi(hmass,s4,s1);

   spinor_field_sub_assign_f(s3,s4);
   tau=sqrt(spinor_field_sqnorm_f(s3));
   if (fabs(tau)>1.e-15) 
     printf("FAILED ");
   else 
     printf("OK ");
   printf("[norm = %e]\n",tau);

   free_field(u_gauge);
#ifndef REPR_FUNDAMENTAL
   free_field(u_gauge_f);
#endif
	free_spinor_field(s1);

   exit(0);
}
