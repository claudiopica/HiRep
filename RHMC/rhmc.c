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
#include "logger.h"

int nhb,nor,nit,nth,nms,level,seed;
double beta;

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
   
   scanf("beta %lf nhb %d nor %d nit %d nth %d nms %d level %d seed %d",
         &beta,&nhb,&nor,&nit,&nth,&nms,&level,&seed);
   fclose(in);
 
}

int main(int argc,char *argv[])
{
   int i,acc;
   rhmc_par rpar;
	 int_par t_par;

	 logger_setlevel(0,10000);

	 /* the following are the logger ID used in the code */
	 /*
			logger_map("ERROR",""); 
			logger_map("TESTING",""); 
			logger_map("MAIN",""); 
			logger_map("SPECLIMITS","");
			logger_map("MaxH2","");
			logger_map("MD_INT","");
			logger_map("INVERTER","");
			logger_map("EVA","");
			logger_map("RAPPROX","");
			logger_map("RHMC","");
			logger_map("FORCE0","");
			logger_map("FORCE_RHMC","");
			logger_stdout("");
		*/

   /* read_cmdline(argc, argv); */

   lprintf("MAIN",0,"Gauge group: SU(%d)\n",NG);
   lprintf("MAIN",0,"Fermion representation: " REPR_NAME " [dim=%d]\n",NF);
   lprintf("MAIN",0,"The lattice size is %dx%d^3\n",T,L);
   
   level=0;
   seed=666;
   rlxs_init(level,seed);
   /*   rlxd_init(level+1,seed); */
   lprintf("MAIN",0,"ranlux: level = %d, seed = %d\n",level,seed); 

   geometry_eo_lexi();
   /* geometry_blocked(); */
   test_geometry();
   u_gauge=alloc_gfield();
#ifndef REPR_FUNDAMENTAL
   u_gauge_f=alloc_gfield_f();
#endif

   rpar.beta = 4.9;
   rpar.mass = -1.3045822102425876010; /* k=0.1855 */
   rpar.nf = 2;
	 rpar.MT_prec = 1.e-15;
	 rpar.MD_prec = 1.e-15;
	 rpar.HB_prec = 1.e-15;
	 rpar.force_prec = 1.e-20;
	 rpar.n_pf = 2;
	 rpar.integrator=&O2MN_multistep;
	 /*rpar.integrator=&leapfrog;*/
	 rpar.MD_par=&t_par;
	 rpar.mshift_solver=&cg_mshift; /* this is not used in the code now */

	 t_par.tlen = 1.;
	 t_par.nsteps = 7;
	 t_par.gsteps = 5;

   printf("Generating a random gauge field... ");fflush(stdout);
   random_u();
   printf("done.\n");
   
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
   printf("Thermalization done.\n");
	 */
/*	 
   read_gauge_field("therm_conf");
*/	 
   
   project_gauge_field();
   represent_gauge_field();

   test_staples();

   init_rhmc(&rpar);

	 lprintf("MAIN",0,"MVM during RHMC initialzation: %ld\n",getMVM());
	 lprintf("MAIN",0,"Initial plaquette: %1.8e\n",avr_plaquette());
   acc=0;
	 for(i=1;i<200;i++) {
		 int rr;
		 float perc;
		 lprintf("MAIN",0,"Trajectory #%d...\n",i);
     rr=update_rhmc();
		 if(rr<0) {
			 lprintf("MAIN",0,"Error in updating the gauge field!!\n");
			 return 1;
		 } else {
			 acc+=rr;
		 }
		 perc=(acc==0)?0.:(float)(100*acc)/(float)i;

		 lprintf("MAIN",0,"Trajectory #%d: %d/%d (%3.4f%%) MVM = %ld\n",i,acc,i,perc,getMVM());
     lprintf("MAIN",0,"Plaquette: %1.8e\n",avr_plaquette());

     if(i%10==0)
       write_gauge_field("therm_conf"); 
  
   }
   write_gauge_field("therm_conf"); 

   free_rhmc();

   free_field(u_gauge);
#ifndef REPR_FUNDAMENTAL
   free_field(u_gauge_f);
#endif

   exit(0);
}
