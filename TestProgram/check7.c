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


void D(suNf_spinor *out, suNf_spinor *in){
   Dphi(hmass,out,in);
}

void H(suNf_spinor *out, suNf_spinor *in){
   g5Dphi(hmass,out,in);
}

void M(suNf_spinor *out, suNf_spinor *in){
   static suNf_spinor tmp[VOLUME];
   g5Dphi(-hmass,tmp,in); 
   g5Dphi(-hmass,out,tmp);
}


int main(int argc,char *argv[])
{
   int i;
   double tau;
   suNf_spinor s1[VOLUME],s2[VOLUME];
   suNf_spinor **res;
   suNf_spinor_dble s3[VOLUME],s4[VOLUME];
   suNf_spinor_dble **resd;

   mshift_par par;
   mshift_par QMRpar;

   int cgiters;

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
   
   printf("Generating a random gauge field... ");
   fflush(stdout);
   random_u();
   printf("done.\n");
   represent_gauge_field();

   set_spinor_len(VOLUME);
   
   par.n = 6;
   par.shift=(double*)malloc(sizeof(double)*(par.n));
   par.err2=1e-8;
   par.max_iter=0;
   res=(suNf_spinor**)malloc(sizeof(suNf_spinor*)*(par.n));
   res[0]=(suNf_spinor*)malloc(sizeof(suNf_spinor)*par.n*VOLUME);
   for(i=1;i<par.n;++i)
     res[i]=res[i-1]+VOLUME;

   resd=(suNf_spinor_dble**)malloc(sizeof(suNf_spinor_dble*)*(par.n));
   resd[0]=(suNf_spinor_dble*)malloc(sizeof(suNf_spinor_dble)*par.n*VOLUME);
   for(i=1;i<par.n;++i)
     resd[i]=resd[i-1]+VOLUME;


   par.shift[0]=+0.1;
   par.shift[1]=-0.21;
   par.shift[2]=+0.05;
   par.shift[3]=-0.01;
   par.shift[4]=-0.15;
   par.shift[5]=-0.05;

   gaussian_spinor_field(&(s1[0]));

   /* TEST g5QMR_M */

   printf("\n");
   printf("Testing g5QMR multishift\n");
   printf("------------------------\n");

   QMRpar.n = 6;
   QMRpar.shift=par.shift;
   QMRpar.err2=1.e-7;
   QMRpar.max_iter=0;
   
   cgiters=g5QMR_mshift(&QMRpar, &D, s1, resd);
   printf("Converged in %d iterations\n",cgiters);

   for(i=0;i<QMRpar.n;++i){
     assign_sd2s(VOLUME,res[i],resd[i]);
      D(s2,res[i]);
      assign_s2sd(VOLUME,s3,s2);
      if(i!=0)
	spinor_field_mul_add_assign_dble_f(s3,-par.shift[i-1],resd[i]);
      assign_s2sd(VOLUME,s4,s1);
      spinor_field_mul_add_assign_dble_f(s3,-1.0,s4);
      tau=spinor_field_sqnorm_dble_f(s3)/(double)spinor_field_sqnorm_f(s1);
      printf("test g5QMR[%d] = %e\n",i,tau);
   }

   exit(0);
}
