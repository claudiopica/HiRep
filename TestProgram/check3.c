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
float beta;

static float hmass=0.1;


int main(int argc,char *argv[])
{
   double tau;
   suNf_spinor s1[VOLUME],s2[VOLUME],s3[VOLUME],s4[VOLUME];

   printf("Gauge group: SU(%d)\n",NG);
   printf("Fermion representation: dim = %d\n",NF);
   printf("The lattice size is %dx%d^3\n",T,L);
   printf("\n");
   
   level=0;
   seed=123;
   printf("ranlux: level = %d, seed = %d\n\n",level,seed); 
   fflush(stdout);
   
   rlxs_init(level,seed);

   geometry_eo_lexi();
   test_geometry();
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
   
   gaussian_spinor_field(&(s1[0]));
   gaussian_spinor_field(&(s2[0]));
   
   tau = 1./sqrt(spinor_field_sqnorm_f(s1));
   spinor_field_mul_f(s1,tau,s1);
   tau = 1./sqrt(spinor_field_sqnorm_f(s2));
   spinor_field_mul_f(s2,tau,s2);

   
   printf("Test new Dirac implementation: ");
   
   g5Dphi_old(hmass,s3,s1);
   g5Dphi(hmass,s4,s1);

   spinor_field_mul_add_assign_f(s3,-1.0,s4);
   tau=spinor_field_sqnorm_f(s3);
   if (fabs(tau)>1.e-7) 
     printf("FAILED ");
   else 
     printf("OK ");
   printf("[norm = %e]\n",tau);

   exit(0);
}
