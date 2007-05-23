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

float mass;

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


void zero_picorr(double *c){
  int i;
  for(i=0;i<T;++i)
   c[i]=0.;
}

void add_picorr(double *c1, float *c2){
  int i;
  for(i=0;i<T;++i)
   c1[i]+=c2[i];
}

void H(suNf_spinor *out, suNf_spinor *in){
   g5Dphi(mass,out,in);
}

int main(int argc,char *argv[])
{
  int i,k,n;

   suNf_spinor **quark_prop;
   float pi_corr[T];
   char propname[256];
   FILE *propfile;

   float m[]={-1.0,-1.6,-1.9};
   int nm=sizeof(m)/sizeof(float);

   read_cmdline(argc, argv); 

   printf("Gauge group: SU(%d)\n",NG);
   printf("Fermion representation: dim = %d\n",NF);
   printf("The lattice size is %dx%d^3\n",T,L);
   printf("beta = %2.4f\n",beta);
   printf("nth  = %d\tNumber of thermalization cycles\n",nth);
   printf("nms  = %d\tNumber of measure cycles\n",nms);
   printf("nit  = %d\tNumber of hb-or iterations per cycle\n",nit);
   printf("nhb  = %d\tNumber of heatbaths per iteration\n",nhb);   
   printf("nor  = %d\tNumber or overrelaxations per iteration\n",nor);
   printf("ranlux: level = %d, seed = %d\n\n",level,seed); 
   printf("Computing quark prop for %d masses: ",nm);
   for(i=0;i<nm;++i)
    printf("%e ",m[i]);
   printf("\n\n");

   fflush(stdout);

   rlxs_init(level,seed);

   geometry_eo_lexi();
   /*geometry_blocked();*/
   test_geometry();

   u_gauge=alloc_gfield();
#ifndef REPR_FUNDAMENTAL
   u_gauge_f=alloc_gfield_f();
#endif

   represent_gauge_field();
   
   /* setup for quark propagators measures */
   quark_prop=(suNf_spinor**)malloc(sizeof(suNf_spinor*)*(nm));
   
   quark_prop[0] = (suNf_spinor*)malloc(sizeof(suNf_spinor)*(nm*VOLUME));
   for (i=1;i<nm;++i) {
     quark_prop[i]=quark_prop[i-1]+VOLUME;
   }
   /* oppure...
   for (i=0;i<4*NF;++i) {
     quark_prop[i]=(suNf_spinor*)malloc(sizeof(suNf_spinor)*(VOLUME));
   }
   */
   sprintf(propname,"quark_prop_%3.5f_%d_%d_%d_%d",beta,T,L,L,L);
   error((propfile = fopen(propname, "wb"))==NULL,1,"Main",
        "Failed to open propagator file for writing");
   fwrite(&nm,(size_t) sizeof(int),1,propfile);
   for(i=0;i<nm;++i)
     fwrite(m+i,(size_t) sizeof(float),1,propfile);

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
   represent_gauge_field();
   
   /* Misure */
   for (i=0;i<nms;++i){ /* nms misure */
     double picorr[nm*T];
     printf("[%d] <p> = %1.6f\n",i,avr_plaquette());
     fflush(stdout);
     
     for (n=0;n<4*NF;++n){
       quark_propagator(n,nm,m,quark_prop);
       for (k=0;k<nm;++k){
         if(n==0) zero_picorr(picorr+k*T);
         pi_correlator(pi_corr, quark_prop[k]);
         add_picorr(picorr+k*T,pi_corr);
       }
       
       /* write propagator on file */
       for (k=0;k<nm;++k) {
	 error(fwrite(quark_prop[k],(size_t) sizeof(suNf_spinor),(size_t)(VOLUME),propfile)!=(VOLUME),1,"Main",
	       "Failed to write quark propagator to file");
       }   
     }
     for(n=0;n<nm;++n){
       printf("[%d] mass=%2.4f pi_corr= ",i,m[n]);
       for(k=0;k<T;++k) {
	 printf("%e ",picorr[k+n*T]);
       }
       printf("\n");
       fflush(stdout);
     }

     for (n=0;n<nit;n++) /* nit updates */
       update(beta,nhb,nor);
     represent_gauge_field();
     
   }


   fclose(propfile);


   free(quark_prop[0]);
   free(quark_prop);

   free_field(u_gauge);

#ifndef REPR_FUNDAMENTAL
   free_field(u_gauge_f);
#endif

   return 0;
}
