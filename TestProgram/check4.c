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


void D(suNf_spinor *out, suNf_spinor *in){
   Dphi(hmass,out,in);
}

void H(suNf_spinor *out, suNf_spinor *in){
   g5Dphi(-hmass,out,in);
}

void M(suNf_spinor *out, suNf_spinor *in){
   suNf_spinor *tmp=alloc_spinor_field_f(1);
   g5Dphi(-hmass,tmp,in); 
   g5Dphi(-hmass,out,tmp);
	 free_field(tmp);
}

void test_herm(spinor_operator S, char *name){
   suNf_spinor *s1, *s2, *s3, *s4;
   double tau;

	 s1=alloc_spinor_field_f(4);
	 s2=s1+VOLUME;
	 s3=s2+VOLUME;
	 s4=s3+VOLUME;

   printf("Test if %s is hermitean: ",name);

   gaussian_spinor_field(s1);
   gaussian_spinor_field(s2);
   S(s3,s1);
   S(s4,s2);

   tau=spinor_field_prod_re_f(s2,s3);
   tau-=spinor_field_prod_re_f(s4,s1);
   tau+=spinor_field_prod_im_f(s2,s3);
   tau-=spinor_field_prod_im_f(s4,s1);
   tau/=sqrt(spinor_field_sqnorm_f(s1));
   tau/=sqrt(spinor_field_sqnorm_f(s2));
   if (fabs(tau)>1.e-7) 
     printf("FAILED ");
   else 
     printf("OK ");
   printf("[norm = %e]\n",tau);

	 free_field(s1);

}


int main(int argc,char *argv[])
{
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
   
	 /*
   gaussian_spinor_field(&(s1[0]));
   gaussian_spinor_field(&(s2[0]));
   
   tau = 1./sqrt(spinor_field_sqnorm_f(s1));
   spinor_field_mul_f(s1,tau,s1);
   tau = 1./sqrt(spinor_field_sqnorm_f(s2));
   spinor_field_mul_f(s2,tau,s2);
*/

   printf("Test hermiticity of the Dirac operator\n");
   printf("--------------------------------------\n");
   
   test_herm(&M,"M");
   test_herm(&H,"H");

   exit(0);
}
