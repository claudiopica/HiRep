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

void read_cmdline(int argc, char*argv[])
{
  int i, ai=0, ao=0;
  FILE *in=NULL, *out=NULL;
  for (i=1;i<argc;i++){
      if (strcmp(argv[i],"-i")==0)
         ai=i+1;
      else if (strcmp(argv[i],"-o")==0)
         ao=i+1;
   }

   if (ao!=0)
      out=freopen(argv[ao],"w",stdout);
   else
      out=stdout;
   
   error(ai==0,1,"suN.c",
         "Syntax: suN -i <input file> [-o <output file>]");

   in=freopen(argv[ai],"r",stdin);
   error(in==NULL,1,"run1.c","Cannot open input file");
   
   scanf("beta %f nhb %d nor %d nit %d nth %d nms %d level %d seed %d",
         &beta,&nhb,&nor,&nit,&nth,&nms,&level,&seed);
   fclose(in);
 
}

void D(suNf_spinor *out, suNf_spinor *in){
   Dphi(0.1,out,in);
}

void H(suNf_spinor *out, suNf_spinor *in){
   g5Dphi(hmass,out,in);
}

void M(suNf_spinor *out, suNf_spinor *in){
   static suNf_spinor tmp[VOLUME];
   g5Dphi(-0.1,tmp,in); 
   g5Dphi(-0.1,out,tmp);
}

void test_herm(spinor_operator S, char *name){
   static suNf_spinor s1[VOLUME],s2[VOLUME],s3[VOLUME],s4[VOLUME];
   double tau;
   printf("Test if %s is hermitean: ",name);

   gaussian_spinor_field(&(s1[0]));
   gaussian_spinor_field(&(s2[0]));
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


}

int main(int argc,char *argv[])
{
   int i;
   double tau;
   suNf_spinor s1[VOLUME],s2[VOLUME],s3[VOLUME],s4[VOLUME];
   suNf_spinor s8[VOLUME];
   suNf_spinor_dble s9[VOLUME],s10[VOLUME];

   mshift_par par;
   MINRES_par MINRESpar2;
   suNf_spinor **res;
   suNf_spinor_dble **resd;
   int cgiters;

   spinor_field_zero_f(s2);
   
   /* read_cmdline(argc, argv); */

   printf("Gauge group: SU(%d)\n",NG);
   printf("Fermion representation: dim = %d\n",NF);
   printf("The lattice size is %dx%d^3\n",T,L);
   printf("\n");
   /*
   printf("beta = %2.4f\n",beta);
   printf("nth  = %d\tNumber of thermalization cycles\n",nth);
   printf("nms  = %d\tNumber of measure cycles\n",nms);
   printf("nit  = %d\tNumber of hb-or iterations per cycle\n",nit);
   printf("nhb  = %d\tNumber of heatbaths per iteration\n",nhb);   
   printf("nor  = %d\tNumber or overrelaxations per iteration\n",nor);
   */
   
   level=0;
   seed=123;
   printf("ranlux: level = %d, seed = %d\n\n",level,seed); 
   fflush(stdout);
   
   rlxs_init(level,seed);
   /*   rlxd_init(level+1,seed); */

   geometry_eo_lexi();
   /*geometry_blocked();*/
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

/* print a matrix */
/*
   av.c1=0.34;
   av.c2=0.5;
   av.c3=0.34;
   av.c4=0.4;
   av.c5=0.;
   av.c6=0.;
   av.c7=0.;
   av.c8=0.;
   _suNg_unit(mg3);
   ExpX(1.,&av,&mg3);
   _suNg_copy(mg,mg3);
   _suNg_dagger(mg2,mg);
   _suNg_times_suNg(mg3,mg,mg2);

  
      printf("(%e,%e)(%e,%e)(%e,%e)\n",mg3.c1_1.re,mg3.c1_1.im,mg3.c1_2.re,mg3.c1_2.im,mg3.c1_3.re,mg3.c1_3.im);
      printf("(%e,%e)(%e,%e)(%e,%e)\n",mg3.c2_1.re,mg3.c2_1.im,mg3.c2_2.re,mg3.c2_2.im,mg3.c2_3.re,mg3.c2_3.im);
      printf("(%e,%e)(%e,%e)(%e,%e)\n\n",mg3.c3_1.re,mg3.c3_1.im,mg3.c3_2.re,mg3.c3_2.im,mg3.c3_3.re,mg3.c3_3.im);
      tau=_suNg_trace_im(mg);
      printf("traccia %e\n",tau);
*/
   
   set_spinor_len(VOLUME);
   
   gaussian_spinor_field(&(s1[0]));
   gaussian_spinor_field(&(s2[0]));
   
   tau = 1./sqrt(spinor_field_sqnorm_f(s1));
   spinor_field_mul_f(s1,tau,s1);
   tau = 1./sqrt(spinor_field_sqnorm_f(s2));
   spinor_field_mul_f(s2,tau,s2);

   
   printf("Test new Dirac implementation: ");
   
   g5Dphi_old(0.1,s3,s1);
   g5Dphi(0.1,s4,s1);

   spinor_field_mul_add_assign_f(s3,-1.0,s4);
   tau=spinor_field_sqnorm_f(s3);
   if (fabs(tau)>1.e-7) 
     printf("FAILED ");
   else 
     printf("OK ");
   printf("[norm = %e]\n",tau);

   test_herm(&M,"M");
   test_herm(&H,"H");

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
   

   /* TEST CG_M */

   printf("Testing CG multishift\n");
   
   cgiters = cg_mshift(&par, &M, s1, res);
   printf("Converged in %d iterations\n",cgiters);

   for(i=0;i<par.n;++i){
      M(s8,res[i]);
      spinor_field_mul_add_assign_f(s8,-par.shift[i],res[i]);
      spinor_field_mul_add_assign_f(s8,-1.0,s1);
      tau=spinor_field_sqnorm_f(s8)/spinor_field_sqnorm_f(s1);
      printf("test cg[%d] = %e\n",i,tau);
   }


   /* TEST MINRES_M */

   hmass=0.1;

   printf("Testing MINRES multishift\n");

	 /*
   MINRESpar.spinorlen=VOLUME;
   MINRESpar.n = 6;
   MINRESpar.shift=par.shift;
   MINRESpar.err2=1.e-8;
   MINRESpar.max_iter=0;
   */

   cgiters=MINRES_mshift(&par, &H, s1, res);
   printf("Converged in %d iterations\n",cgiters);

   for(i=0;i<par.n;++i){
      H(s8,res[i]);
      if(i!=0)
	spinor_field_mul_add_assign_f(s8,-par.shift[i-1],res[i]);
      spinor_field_mul_add_assign_f(s8,-1.0,s1);
      tau=spinor_field_sqnorm_f(s8)/spinor_field_sqnorm_f(s1);
      printf("test MINRES[%d] = %e\n",i,tau);
   }
   
   /* TEST MINRES_M */

   printf("Testing MINRES \n");

   MINRESpar2.err2=1.e-8;
   MINRESpar2.max_iter=0;
   
   cgiters=MINRES(&MINRESpar2, &H, s1, res[0],0);
   for(i=1;i<par.n;++i){
     hmass=0.1-par.shift[i-1];
     cgiters+=MINRES(&MINRESpar2, &H, s1, res[i],res[i-1]);
   }
   printf("Converged in %d iterations\n",cgiters);

   hmass=0.1;
   for(i=0;i<par.n;++i){
     if(i!=0)
       hmass=0.1-par.shift[i-1];
     H(s8,res[i]);
      spinor_field_mul_add_assign_f(s8,-1.0,s1);
      tau=spinor_field_sqnorm_f(s8)/spinor_field_sqnorm_f(s1);
      printf("test MINRES[%d] = %e\n",i,tau);
   }

   /* TEST g5QMR_M */

   printf("Testing g5QMR multishift\n");

	 /*
   QMRpar.spinorlen=VOLUME;
   QMRpar.n = 6;
   QMRpar.shift=par.shift;
   QMRpar.err2=1.e-7;
   QMRpar.max_iter=0;
   */

   cgiters=g5QMR_mshift(&par, &D, s1, resd);
   printf("Converged in %d iterations\n",cgiters);

   for(i=0;i<par.n;++i){
     assign_sd2s(VOLUME,res[i],resd[i]);
      D(s8,res[i]);
      assign_s2sd(VOLUME,s9,s8);
      if(i!=0)
	spinor_field_mul_add_assign_dble_f(s9,-par.shift[i-1],resd[i]);
      assign_s2sd(VOLUME,s10,s1);
      spinor_field_mul_add_assign_dble_f(s9,-1.0,s10);
      tau=spinor_field_sqnorm_dble_f(s9)/(double)spinor_field_sqnorm_f(s1);
      printf("test g5QMR[%d] = %e\n",i,tau);
   }

   exit(0); /* biCGStab not working... yet */

   /* TEST BICGSTAB_M */

   printf("Testing BiCGstab multishift\n");

	 /*
   BiCGpar.n = 3;
   BiCGpar.shift=par.shift;
   BiCGpar.err2=1.e-8;
   BiCGpar.max_iter=0;
   */
   cgiters=HBiCGstab_mshift(&par, &M,s1, res);
   printf("Converged in %d iterations\n",cgiters);

   for(i=0;i<par.n;++i){
      M(s8,res[i]);
      if(i!=0)
	spinor_field_mul_add_assign_f(s8,-par.shift[i-1],res[i]);
      spinor_field_mul_add_assign_f(s8,-1.0,s1);
      tau=spinor_field_sqnorm_f(s8)/spinor_field_sqnorm_f(s1);
      printf("test BiCGStab[%d] = %e\n",i,tau);
   }
   


   free_field(u_gauge);
#ifndef REPR_FUNDAMENTAL
   free_field(u_gauge_f);
#endif

   exit(0);
}
