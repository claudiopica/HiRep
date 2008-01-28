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

#define NAME_SIZE 128

int ai=0,ao=0,ic=0,ia=0;
int nf,n_pf,nsteps,gsteps,level,seed,ntraj,wrt_cnfg;
double beta,mass,MT_prec,MD_prec,HB_prec,force_prec,tlen;
char log_dir[NAME_SIZE],dat_dir[NAME_SIZE],cnfg_dir[NAME_SIZE];
char log_file[NAME_SIZE],dat_file[NAME_SIZE],log_save[NAME_SIZE],dat_save[NAME_SIZE];
char cnfg_file[NAME_SIZE],end_file[NAME_SIZE];
char base[NAME_SIZE],cnfg[NAME_SIZE];


int copy_file(char *in,char *out)
{
   int c;
   FILE *fin,*fout;

   fin=fopen(in,"rb");

   if (fin==NULL)
   {
      error(1,1,"copy_file [rhmc.c]","Unable to open input file");
      return 1;
   }

   fout=fopen(out,"wb");

   if (fout==NULL)
   {
      error(1,1,"copy_file [rhmc.c]","Unable to open output file");
      fclose(fin);
      return 1;
   }

   c=getc(fin);

   while (feof(fin)==0)
   {
      putc(c,fout);
      c=getc(fin);
   }

   if ((ferror(fin)==0)&&(ferror(fout)==0))
   {
      fclose(fin);
      fclose(fout);
      return 0;
   }
   else
   {
      error(1,1,"copy_file [rhmc.c]","Read or write error");
      fclose(fin);
      fclose(fout);
      return 1;
   }
}


static int find_opt(int argc,char *argv[],char *opt)
{
   int k;

   for (k=1;k<argc;k++)
      if (strcmp(argv[k],opt)==0)
         return (k+1);

   return 0;
}

void read_cmdline(int argc, char*argv[])
{
  FILE *in=NULL,*out=NULL;

  ai=find_opt(argc,argv,"-i");
  ao=find_opt(argc,argv,"-o");
  ic=find_opt(argc,argv,"-c");
  ia=find_opt(argc,argv,"-a");

  if (ao!=0)
    out=freopen(argv[ao],"w",stdout);
  else
    out=stdout;
  
  error(ai==0,1,"rhmc.c",
	"Syntax: rhmc -i <input file> [-o <output file>] [-c <configuration file> [-a]]");
  error(((ic>0)&&(ic==(argc-1)))||((ia>0)&&(ic==0)),1,"main [run3.c]",
	"Syntax: rhmc -i <input file> [-o <output file>] [-c <configuration file> [-a]]");

  in=freopen(argv[ai],"r",stdin);
  error(in==NULL,1,"rhmc.c","Cannot open input file");
  
  scanf("beta %lf mass %lf nf %d MT_prec %lf MD_prec %lf HB_prec %lf\n",
	&beta,&mass,&nf,&MT_prec,&MD_prec,&HB_prec);
  scanf("force_prec %lf n_pf %d tlen %lf  nsteps %d gsteps %d level %d seed %d ntraj %d\n",
	&force_prec,&n_pf,&tlen,&nsteps,&gsteps,&level,&seed,&ntraj);
  scanf("wrt_cnfg %d\n",
	&wrt_cnfg);
  scanf("log_dir %s dat_dir %s cnfg_dir %s\n",
	log_dir,dat_dir,cnfg_dir);
  
  fclose(in);
  
  if (ic>0)
    strcpy(cnfg,argv[ic]);
}

int main(int argc,char *argv[])
{
   int i,acc,icnfg,n,iend;
   rhmc_par rpar;
   int_par t_par;
   float kappa;
   FILE *end,*dat,*log;


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

   read_cmdline(argc, argv);

   kappa=1.0/(2*mass+8.0);
   sprintf(base,"%dx%dx%dx%dNc%dNf%db%1.2fk%.5f",
           T,L,L,L,NG,nf,beta,kappa);

   sprintf(log_file,"%s/%s.log",log_dir,base);
   sprintf(dat_file,"%s/%s.dat",dat_dir,base);
   sprintf(end_file,"%s/%s.end",log_dir,base);
   sprintf(log_save,"%s~",log_file);
   sprintf(dat_save,"%s~",dat_file);

   geometry_eo_lexi();
   test_geometry();
   u_gauge=alloc_gfield();
#ifndef REPR_FUNDAMENTAL
   u_gauge_f=alloc_gfield_f();
#endif
   
   if (ia==0)
     {
       end=fopen(log_file,"r");
       dat=fopen(dat_file,"rb");
       
       error((end!=NULL)||(dat!=NULL),1,"main [rhmc.c]",
	     "Attempt to overwrite old .log or .dat file");
       icnfg=1;
     }
   else
     {
       error(strstr(cnfg,base)==NULL,1,"main [rhmc.c]",
	     "New and old run parameters do not match");
       n=strlen(base);
       sscanf(cnfg+n,"n%d",&icnfg);
       icnfg+=1;
     }
   
   if (ia==0)
     log=freopen(log_file,"w",stdout);
   else
     log=freopen(log_file,"a",stdout);
   
   error(log==NULL,1,"main [rhmc.c]","Unable to open output file");

   if (ic>0)
     {
       sprintf(cnfg_file,"%s/%s",cnfg_dir,cnfg);
       read_gauge_field(cnfg_file);
       
       if (ia>0)
	 {
	   sprintf(cnfg_file,"%s/%sn%d",cnfg_dir,base,icnfg);
	   dat=fopen(cnfg_file,"rb");
	   error(dat!=NULL,1,"main [rhmc.c]",
		 "Initial configuration is not the last from the previous run");
	 }
     }
   else
     {
       printf("Generating a random gauge field... ");fflush(stdout);
       random_u();
       printf("done.\n");
     }
   
   if (ia==0)
     {
       printf("\n");
       printf("Generation of gauge field configurations\n");
       printf("----------------------------------------\n\n");
       
       if (ic>0)
	 printf("New run, start from configuration %s\n\n",cnfg);
       else
	 printf("New run, start from random configuration\n\n");
     }
   else
     printf("Continuation run, start from configuration %s\n\n",cnfg);

   lprintf("MAIN",0,"Gauge group: SU(%d)\n",NG);
   lprintf("MAIN",0,"Fermion representation: " REPR_NAME " [dim=%d]\n",NF);
   lprintf("MAIN",0,"The lattice size is %dx%d^3\n",T,L);
   
   rlxd_init(level,seed);
   lprintf("MAIN",0,"ranlux: level = %d, seed = %d\n",level,seed); 

   rpar.beta = beta;
   rpar.mass = mass;
   rpar.nf = nf;
   rpar.MT_prec = MT_prec;
   rpar.MD_prec = MD_prec;
   rpar.HB_prec = HB_prec;
   rpar.force_prec = force_prec;
   rpar.n_pf = n_pf;
   rpar.integrator=&O2MN_multistep;
   rpar.MD_par=&t_par;
   rpar.mshift_solver=&cg_mshift; /* this is not used in the code now */
   
   t_par.tlen = tlen;
   t_par.nsteps = nsteps;
   t_par.gsteps = gsteps;

   project_gauge_field();
   represent_gauge_field();

   test_staples();

   init_rhmc(&rpar);
   
   lprintf("MAIN",0,"MVM during RHMC initialzation: %ld\n",getMVM());
   lprintf("MAIN",0,"Initial plaquette: %1.8e\n",avr_plaquette());
   acc=0;
   for(i=icnfg;i<ntraj;i++) {
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
     perc=(acc==0)?0.:(float)(100*acc)/(float)(i-icnfg+1);
     
     lprintf("MAIN",0,"Trajectory #%d: %d/%d (%3.4f%%) MVM = %ld\n",i,acc,i,perc,getMVM());
     lprintf("MAIN",0,"Plaquette: %1.8e\n",avr_plaquette());
     
		 if((i%wrt_cnfg)==0)
		 {
			 sprintf(cnfg_file,"%s/%sn%d",cnfg_dir,base,i);
			 write_gauge_field(cnfg_file); 

			 printf("Configuration no %d saved and exported\n\n",i);
			 end=fopen(end_file,"r");

			 if (end!=NULL)
			 {
				 fclose(end);
				 remove(end_file);
				 iend=1;
				 printf("End flag set, run stopped\n\n");
			 }

			 fflush(log);
			 copy_file(log_file,log_save);
			 /* copy_file(dat_file,dat_save); */
		 }
	 }

	 free_rhmc();
	 free_field(u_gauge);
#ifndef REPR_FUNDAMENTAL
	 free_field(u_gauge_f);
#endif

	 exit(0);
}
