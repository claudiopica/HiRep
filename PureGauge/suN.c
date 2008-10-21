/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *   
* All rights reserved.                                                      * 
\***************************************************************************/

/*******************************************************************************
*
* Computation of the average plaquette
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
#include "logger.h"
#include "observables.h"


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
   int i,n,flag;
   double *p,pbar,sig,tau;

   read_cmdline(argc, argv);

	 logger_setlevel(0,10000);

   printf("Gauge group: SU(%d)\n",NG);

	 /*
   geometry_blocked();
	 */
	 read_input("test_input");
	 setup_process();
	 define_geometry();
   /*test_geometry();*/

   printf("The lattice size is %dx%dx%dx%d\n",T,X,Y,Z);
   printf("\n");
   printf("beta = %2.4f\n",beta);
   printf("nth  = %d\tNumber of thermalization cycles\n",nth);
   printf("nms  = %d\tNumber of measure cycles\n",nms);
   printf("nit  = %d\tNumber of hb-or iterations per cycle\n",nit);
   printf("nhb  = %d\tNumber of heatbaths per iteration\n",nhb);   
   printf("nor  = %d\tNumber or overrelaxations per iteration\n",nor);
   printf("ranlux: level = %d, seed = %d\n\n",level,seed); 
   fflush(stdout);
   
   rlxd_init(level,seed);
   /*   rlxd_init(level+1,seed); */


   u_gauge=alloc_gfield();

   /*read_gauge_field_single_biagio("Config_biagio"); */

   /* random_u(); */

   p=malloc(nms*sizeof(double)); /* array per le misure */

   /* Termalizzazione */
   for (i=0;i<nth;++i) {
      update(beta,nhb,nor);
      if ((i%10)==0)
              printf("%d",i);
      else
              printf(".");
      fflush(stdout);
   }
   if(i) printf("%d\nThemalization done.\n",i);

   
   /* Misure */
   for (i=0;i<nms;i++){ /* nms misure */
     p[i]=avr_plaquette();
     printf("[%d] <p> = %1.6f\n",i,p[i]);
     fflush(stdout);

     for (n=0;n<nit;n++) /* nit updates */
       update(beta,nhb,nor);
     
     /* Plaquette media */
   }

   printf("[%d] <p> = %1.6f\n",i,avr_plaquette());
   fflush(stdout);

   pbar=average(nms,p);
   sig=sigma(nms,p,&tau,&flag);

   if (flag!=0)
   {
      printf("Warning: unreliable error estimation ");
      printf("(data series is too short)\n\n");
   }
   
   printf("<p>   = %1.6f\n",pbar);
   printf("sigma = %1.6f\n",sig);
   printf("tau   = %1.2f [iterations]\n\n",tau);


   write_gauge_field("Config");

   free_field(u_gauge);


   return 0;
}
