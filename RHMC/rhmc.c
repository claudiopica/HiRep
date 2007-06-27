/*******************************************************************************
*
* Main RHMC program
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

int main(int argc,char *argv[])
{
   int i,acc;
   rhmc_par rpar;
	 int_par t_par;
   
   /* read_cmdline(argc, argv); */

   printf("Gauge group: SU(%d)\n",NG);
   printf("Fermion representation: " REPR_NAME " [dim=%d]\n",NF);
   printf("The lattice size is %dx%d^3\n",T,L);
   /*
   printf("beta = %2.4f\n",beta);
   printf("nth  = %d\tNumber of thermalization cycles\n",nth);
   printf("nms  = %d\tNumber of measure cycles\n",nms);
   printf("nit  = %d\tNumber of hb-or iterations per cycle\n",nit);
   printf("nhb  = %d\tNumber of heatbaths per iteration\n",nhb);   
   printf("nor  = %d\tNumber or overrelaxations per iteration\n",nor);
   */
   
   level=0;
   seed=666;
   printf("ranlux: level = %d, seed = %d\n\n",level,seed); 
   fflush(stdout);
   
   rlxs_init(level,seed);
   /*   rlxd_init(level+1,seed); */

   /* test gauss */
   /* {
     int l;
     double *r;
     r=malloc(sizeof(double)*100000);
     gauss_dble(r,100000);
     for(l=0;l<100000;++l) r[l]/=sqrt(2.);
     printf("Avr: %e", average(100000,r));
     printf(" s0: %e\n",sqrt(100000.)*sigma0(100000,r));
     free(r);
   } */

   geometry_eo_lexi();
   /* geometry_blocked(); */
   test_geometry();
   u_gauge=alloc_gfield();
#ifndef REPR_FUNDAMENTAL
   u_gauge_f=alloc_gfield_f();
#endif

   /* 
   printf("Generating a random gauge field... ");fflush(stdout);
   random_u();
   printf("done.\n");
   */

   project_gauge_field();
   represent_gauge_field();

   rpar.beta = 5.6;
   rpar.mass = -1.4;
   rpar.nf = 2;
	 rpar.MT_prec = 1.e-10;
	 rpar.MD_prec = 1.e-6;
	 rpar.HB_prec = 1.e-10;
	 rpar.force_prec = 1.e-10;
	 rpar.n_pf = 2;
	 rpar.integrator=&O2MN_multistep;
	 rpar.MD_par=&t_par;
	 rpar.mshift_solver=&cg_mshift;

	 t_par.tlen = 5.;
	 t_par.nsteps = 30;
	 t_par.gsteps = 3;


	 /*
   printf("Thermalizing with CM.\n");
   for(i=0;i<50;i++) {
     update(rpar.beta,1,10);
     if(i%10) {
       printf(".");
     } else {
       printf("%d",i);
     }
     fflush(stdout);
   }
   represent_gauge_field();
   printf("Thermalization done.\n");
   */
	 /*
   read_gauge_field_single("therm_conf");
	 represent_gauge_field();
	 */
   test_staples();

	 printf("Initializing RHMC...\n");
   init_rhmc(&rpar);
	 printf("Initializing RHMC... done.\n");

   acc=0;
	 for(i=0;i<200;i++) {
		 float perc=(acc==0)?0.:(float)(100*acc)/(float)i;
		 printf("[Trajectory #%d: %d/%d (%3.4f%%)]\n",i,acc,i,perc);
     printf("[Plaq: %1.8e]\n",avr_plaquette());fflush(stdout);
     acc+=update_rhmc();
/*
     if(i%10==0)
       write_gauge_field_single("therm_conf_5.6"); 
*/  
   }
   write_gauge_field_single("therm_conf"); 



   free_rhmc();
	 printf("RHMC ended. \n");
   free_field(u_gauge);
#ifndef REPR_FUNDAMENTAL
   free_field(u_gauge_f);
#endif

   exit(0);
}
