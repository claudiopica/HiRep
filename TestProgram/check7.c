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
#include "logger.h"

int nhb,nor,nit,nth,nms,level,seed;
double beta;

static double hmass=0.1;


void D(spinor_field *out, spinor_field *in){
   Dphi(hmass,out,in);
}

void H(spinor_field *out, spinor_field *in){
   g5Dphi(hmass,out,in);
}

static spinor_field *tmp;
void M(spinor_field *out, spinor_field *in){
   g5Dphi(-hmass,tmp,in); 
   g5Dphi(-hmass,out,tmp);
}


int main(int argc,char *argv[])
{
   int i;
   double tau;
   spinor_field *s1, *s2;
   spinor_field *res;

   mshift_par par;

   int cgiters;

   printf("Gauge group: SU(%d)\n",NG);
   printf("Fermion representation: dim = %d\n",NF);
   geometry_eo_lexi();
   test_geometry();
   printf("The lattice size is %dx%dx%dx%d\n",T,X,Y,Z);
   printf("\n");
   
   level=0;
   seed=123;
   printf("ranlux: level = %d, seed = %d\n\n",level,seed); 
   fflush(stdout);
   
   rlxs_init(level,seed);

	 /* logger_setlevel(0,10015); */

   u_gauge=alloc_gfield();
#ifndef REPR_FUNDAMENTAL
   u_gauge_f=alloc_gfield_f();
#endif
   represent_gauge_field();

   printf("Generating a random gauge field... ");
   fflush(stdout);
   random_u();
   printf("done.\n");
   represent_gauge_field();

   set_spinor_len(VOLUME);
   
   par.n = 6;
   par.shift=(double*)malloc(sizeof(double)*(par.n));
   par.err2=1.e-28;
	 par.max_iter=0;
   res=alloc_spinor_field_f(par.n+3);
	s1=res+par.n;
	s2=s1+1;
	tmp=s2+1;


   par.shift[0]=+0.1;
   par.shift[1]=-0.21;
   par.shift[2]=+0.05;
   par.shift[3]=-0.01;
   par.shift[4]=-0.15;
   par.shift[5]=-0.05;
   

	 gaussian_spinor_field(s1); 
	 tau=spinor_field_sqnorm_f(s1);
	 printf("Norma iniziale: %e\n",tau);

   /* TEST g5QMR_M */

   printf("\n");
   printf("Testing g5QMR multishift\n");
   printf("------------------------\n");

   cgiters=g5QMR_mshift(&par, &D, s1, res);
   printf("Converged in %d iterations\n",cgiters);

   for(i=0;i<par.n;++i){
      D(s2,&res[i]);
			spinor_field_mul_add_assign_f(s2,-par.shift[i],&res[i]);
      spinor_field_sub_assign_f(s2,s1);
      tau=spinor_field_sqnorm_f(s2)/spinor_field_sqnorm_f(s1);
      printf("test g5QMR[%d] = %e (req. %e)\n",i,tau,par.err2);
   }
	 
	free_spinor_field(res);
	 free(par.shift);

   exit(0);
}
