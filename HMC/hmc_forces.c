/****************************************************************************
* Copyright (c) 2014, Ari Hietanen
* All rights reserved.                                                      * 
\***************************************************************************/

/*******************************************************************************
*
* Main HMC program
*
*******************************************************************************/

#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "io.h"
#include "ranlux.h"
#include "geometry.h"
#include "update.h"
#include "global.h"
#include "observables.h"
#include "dirac.h"
#include "logger.h"
#include "hmc_utils.h"
#include "memory.h"
#include "communications.h"
#include "observables.h"
#include "utils.h"
#include "spectrum.h"
#include "cinfo.c"
#include "representation.h"
#include "linear_algebra.h"




/* flow control variable */
hmc_flow flow=init_hmc_flow(flow);


char input_filename[256] = "input_file";
char output_filename[256] = "out_0";
char error_filename[256] = "err_0";
char cnfg_filename[256]="";
char list_filename[256]="";

static void read_cmdline(int argc, char* argv[]) {
  int i, ai=0, ao=0, am=0, requested=1,ac=0,al=0;
  FILE *list=NULL;

  for (i=1;i<argc;i++) {
    if (strcmp(argv[i],"-i")==0) {ai=i+1;requested+=2;}
    else if (strcmp(argv[i],"-o")==0) {ao=i+1;requested+=2;}
    else if (strcmp(argv[i],"-m")==0) {am=i;requested+=1;}
    else if (strcmp(argv[i],"-c")==0) {ac=i+1;requested+=2;}
    else if (strcmp(argv[i],"-l")==0) {al=i+1;;requested+=2;}
  }

  if (am != 0) {
    print_compiling_info();
    exit(0);
  }

  error(argc!=requested,1,"read_cmdline [hmc.c]",
      "Arguments: -l <list file> | -c <cnfg file> [-i <input file>] [-o <output file>]  [-m]");

  if (ao!=0) strcpy(output_filename,argv[ao]);
  if (ai!=0) strcpy(input_filename,argv[ai]);
  if(ac != 0) {
    strcpy(cnfg_filename,argv[ac]);
    strcpy(list_filename,"");
  } else if(al != 0) {
    strcpy(list_filename,argv[al]);
    error((list=fopen(list_filename,"r"))==NULL,1,"parse_cmdline [mk_mesons.c]" ,
          "Failed to open list file\n");
    error(fscanf(list,"%s",cnfg_filename)==0,1,"parse_cmdline [mk_mesons.c]" ,
          "Empty list file\n");
    fclose(list);
  }
}


enum { UNKNOWN_CNFG, DYNAMICAL_CNFG, QUENCHED_CNFG };

typedef struct {
  char string[256];
  int t, x, y, z;
  int nc, nf;
  double b, m;
  int n;
  int type;
} filename_t;


int parse_cnfg_filename(char* filename, filename_t* fn) {
  int hm;
  char *tmp = NULL;
  char *basename;

  basename = filename;
  while ((tmp = strchr(basename, '/')) != NULL) {
    basename = tmp+1;
  }            

#ifdef REPR_FUNDAMENTAL
#define repr_name "FUN"
#elif defined REPR_SYMMETRIC
#define repr_name "SYM"
#elif defined REPR_ANTISYMMETRIC
#define repr_name "ASY"
#elif defined REPR_ADJOINT
#define repr_name "ADJ"
#endif
  hm=sscanf(basename,"%*[^_]_%dx%dx%dx%d%*[Nn]c%dr" repr_name "%*[Nn]f%db%lfm%lfn%d",
      &(fn->t),&(fn->x),&(fn->y),&(fn->z),&(fn->nc),&(fn->nf),&(fn->b),&(fn->m),&(fn->n));
  if(hm==9) {
    fn->m=-fn->m; /* invert sign of mass */
    fn->type=DYNAMICAL_CNFG;
    return DYNAMICAL_CNFG;
  }
#undef repr_name

  double kappa;
  hm=sscanf(basename,"%dx%dx%dx%d%*[Nn]c%d%*[Nn]f%db%lfk%lfn%d",
      &(fn->t),&(fn->x),&(fn->y),&(fn->z),&(fn->nc),&(fn->nf),&(fn->b),&kappa,&(fn->n));
  if(hm==9) {
    fn->m = .5/kappa-4.;
    fn->type=DYNAMICAL_CNFG;
    return DYNAMICAL_CNFG;
  }

  hm=sscanf(basename,"%dx%dx%dx%d%*[Nn]c%db%lfn%d",
      &(fn->t),&(fn->x),&(fn->y),&(fn->z),&(fn->nc),&(fn->b),&(fn->n));
  if(hm==7) {
    fn->type=QUENCHED_CNFG;
    return QUENCHED_CNFG;
  }

  hm=sscanf(basename,"%*[^_]_%dx%dx%dx%d%*[Nn]c%db%lfn%d",
      &(fn->t),&(fn->x),&(fn->y),&(fn->z),&(fn->nc),&(fn->b),&(fn->n));
  if(hm==7) {
    fn->type=QUENCHED_CNFG;
    return QUENCHED_CNFG;
  }

  fn->type=UNKNOWN_CNFG;
  return UNKNOWN_CNFG;
}


int main(int argc,char *argv[]) {
  char sbuf[128];
  
#ifndef MEASURE_FORCE
  read_cmdline(argc,argv);
  setup_process(&argc,&argv);
  read_input(logger_var.read,input_filename);
  logger_set_input(&logger_var);
  if (PID!=0) { 
    logger_disable(); }   /* disable logger for MPI processes != 0 */
  else {
    FILE* stderrp;
    sprintf(sbuf,">>%s",output_filename);  logger_stdout(sbuf);
    stderrp=freopen(error_filename,"w",stderr);
    error(stderrp==NULL,1,"main [hmc.c]",
	  "Cannot redirect the stderr");
  }
  lprintf("MAIN",0,"Compiled with macros: %s\n",MACROS);
  error(1,1,"MAIN","Need to set MEASURE_FORCE macro for program to measure forces\n"); 
  (void) CI_svnrevision;
}
#else
  FILE* list;
  filename_t fpars;

  read_cmdline(argc,argv);
  
  /* setup process communications */
  setup_process(&argc,&argv);
  
  /* read global variables file */
  read_input(glb_var.read,input_filename);
  
  setup_replicas();
  
  /* logger setup */
  read_input(logger_var.read,input_filename);
  logger_set_input(&logger_var);
  if (PID!=0) { logger_disable(); }   /* disable logger for MPI processes != 0 */
  else {
    FILE* stderrp;
    sprintf(sbuf,">>%s",output_filename);  logger_stdout(sbuf);
    stderrp=freopen(error_filename,"w",stderr);
    error(stderrp==NULL,1,"main [hmc.c]",
	  "Cannot redirect the stderr");
  }
  
  lprintf("MAIN",0,"Compiled with macros: %s\n",MACROS);
  lprintf("MAIN",0,"[RepID: %d][world_size: %d]\n[MPI_ID: %d][MPI_size: %d]\n\n",RID,WORLD_SIZE,MPI_PID,MPI_WORLD_SIZE);

  //  lprintf("MAIN",0,"Logger lelvel: %d\n",logger_getlevel(0));
  
  /* setup lattice geometry */
  if (geometry_init() == 1) { finalize_process(); return 0; }
  geometry_mpi_eo();
  /* test_geometry_mpi_eo(); */ 
  /* setup random numbers */
  read_input(rlx_var.read,input_filename);
  //slower(rlx_var.rlxd_start); //convert start variable to lowercase
  if(strcmp(rlx_var.rlxd_start,"continue")==0 && rlx_var.rlxd_state[0]!='\0') {
    /*load saved state*/
    lprintf("MAIN",0,"Loading rlxd state from file [%s]\n",rlx_var.rlxd_state);
    read_ranlxd_state(rlx_var.rlxd_state);
  } else {
    lprintf("MAIN",0,"RLXD [%d,%d]\n",rlx_var.rlxd_level,rlx_var.rlxd_seed+MPI_PID);
    rlxd_init(rlx_var.rlxd_level,rlx_var.rlxd_seed+MPI_PID); /* use unique MPI_PID to shift seeds */
  }

#ifdef GAUGE_SUN
  lprintf("MAIN",0,"Gauge group: SU(%d)\n",NG);
#elif GAUGE_SON
  lprintf("MAIN",0,"Gauge group: SO(%d)\n",NG);
#else
  lprintf("MAIN",0,"Default gauge group: SU(%d)\n",NG);
#endif
  lprintf("MAIN",0,"Fermion representation: " REPR_NAME " [dim=%d]\n",NF);

  /* Init Monte Carlo */
  init_mc(&flow, input_filename);
  parse_cnfg_filename(cnfg_filename,&fpars);
  list=NULL;
  if(strcmp(list_filename,"")!=0) {
    error((list=fopen(list_filename,"r"))==NULL,1,"main [mk_mesons.c]" ,
	"Failed to open list file\n");
  }

  lprintf("MAIN",0,"Compiled with macros: %s\n",MACROS); 
  lprintf("MAIN",0,"PId =  %d [world_size: %d]\n\n",PID,WORLD_SIZE); 
  lprintf("MAIN",0,"input file [%s]\n",input_filename); 
  lprintf("MAIN",0,"output file [%s]\n",output_filename); 
  if (list_filename[0]!=0) lprintf("MAIN",0,"list file [%s]\n",list_filename); 
  else lprintf("MAIN",0,"cnfg file [%s]\n",cnfg_filename); 

  parse_cnfg_filename(cnfg_filename,&fpars);

  GLB_T=fpars.t; GLB_X=fpars.x; GLB_Y=fpars.y; GLB_Z=fpars.z;
  error(fpars.type==UNKNOWN_CNFG,1,"measure_spectrum.c","Bad name for a configuration file");
  error(fpars.nc!=NG,1,"measure_spectrum.c","Bad NG");

  int i=0;
  while(1){
    struct timeval start, end, etime; /* for timing */
    double times[num_mon()];
    if (force_ave==NULL){
      force_ave = (double*) malloc(num_mon()*sizeof(double));
      force_max = (double*) malloc(num_mon()*sizeof(double));
      n_inv_iter = (int*) malloc(num_mon()*sizeof(int));
    }
    for (int k=0;k<num_mon();k++){
      force_ave[k]=0.0;
      force_max[k]=0.0;
      n_inv_iter[k]=0;
    }


    if(list!=NULL)
      if(fscanf(list,"%s",cnfg_filename)==0 || feof(list)) break;

    i++; 

    lprintf("MAIN",0,"Configuration from %s\n", cnfg_filename);
    read_gauge_field(cnfg_filename);
    represent_gauge_field();

    lprintf("TEST",0,"<p> %1.6f\n",avr_plaquette());
    full_plaquette();

    

    gettimeofday(&start,0);

    corret_pf_dist_hmc();
    
    gettimeofday(&end,0);
    timeval_subtract(&etime,&end,&start);
    lprintf("MAIN",0,"Configuration %d: Time to correct pseudofermion dist: %ld sec %ld usec\n",i,etime.tv_sec,etime.tv_usec);

    gettimeofday(&start,0);
    calc_one_force(0);
    gettimeofday(&end,0);
    timeval_subtract(&etime,&end,&start);
    lprintf("MAIN",0,"Time to calculate gauge force: %ld sec %ld usec\n",etime.tv_sec,etime.tv_usec);
    times[0]=etime.tv_sec + (double) (etime.tv_usec)/1.0e6;

    gettimeofday(&start,0);
    calc_one_force(1);
    gettimeofday(&end,0);
    timeval_subtract(&etime,&end,&start);
    lprintf("MAIN",0,"Time to calculate fermion force: %ld sec %ld usec\n",etime.tv_sec,etime.tv_usec);
    times[1]=etime.tv_sec + (double) (etime.tv_usec)/1.0e6;
    for (int k=2;k<num_mon();k++){
      gettimeofday(&start,0);
      calc_one_force(k);
      gettimeofday(&end,0);
      timeval_subtract(&etime,&end,&start);
      lprintf("MAIN",0,"Time to calculate HB force %d: %ld sec %ld usec\n",k,etime.tv_sec,etime.tv_usec);
      times[k]=etime.tv_sec + (double) (etime.tv_usec)/1.0e6;
    }

    lprintf("FORCE_SUMMARY",10,"Fermion: is the first monomial defined in the input file and Hasen: are the following ones upto number of monomials\n");
    lprintf("FORCE_SUMMARY",0,"%d ave Gauge: %1.6f, Fermion: %1.6f",i,force_ave[0],force_ave[1]);

    for (int k=2;k<num_mon();++k){
      lprintf("FORCE_SUMMARY",0,", Hasen %d: %1.6f",k-2,force_ave[k]);
    }
    lprintf("FORCE_SUMMARY",0,"\n");
    lprintf("FORCE_SUMMARY",0,"%d max Gauge: %1.6f, Fermion: %1.6f",i,force_max[0],force_max[1]);
    for (int k=2;k<num_mon();++k){
      lprintf("FORCE_SUMMARY",0,", Hasen %d: %1.6f",k-2,force_max[k]);
    }
    lprintf("FORCE_SUMMARY",0,"\n");
    lprintf("INV_SUMMARY",0,"%d Iterations in fermion: %d ",i,n_inv_iter[0]);
    for (int k=2;k<num_mon();++k){
      lprintf("INV_SUMMARY",0," Hasenbusch %d: %d",k-2,n_inv_iter[k-1]);
    }
    lprintf("INV_SUMMARY",0,"\n");

    lprintf("TIME_SUMMARY",0,"%d Time in gauge: %g fermion: %g",i,times[0],times[1]);
    for (int k=2;k<num_mon();++k){
      lprintf("TIME_SUMMARY",0," Hasenbusch %d: %1.6f",k-2,times[k]);
    }
  }
  
  free(force_ave);
  free(force_max);
  free(n_inv_iter);
  
  
  /* finalize Monte Carlo */
  end_mc();
  
  /* close communications */
  finalize_process();
  
  return 0;
  
}
#endif
