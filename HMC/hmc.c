/****************************************************************************
* Copyright (c) 2008, Claudio Pica                                          *   
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


/* Mesons parameters */
typedef struct _input_mesons {
  char make[256];
  double precision;
  int nhits;
  double mesmass; /* valence mass */
  
  /* for the reading function */
  input_record_t read[5];
  
} input_mesons;

#define init_input_mesons(varname) \
{ \
  .read={\
    {"make mesons", "mes:make = %s", STRING_T, (varname).make},\
    {"inverter precision", "mes:precision = %lf", DOUBLE_T, &(varname).precision},\
    {"number of noisy sources per cnfg", "mes:nhits = %d", INT_T, &(varname).nhits},\
    {"valence mass", "mes:mass = %lf", DOUBLE_T, &(varname).mesmass},\
    {NULL, NULL, INT_T, NULL}\
  }\
}


input_mesons mes_var = init_input_mesons(mes_var);

/* Polyakov-loop parameters */
typedef struct _input_polyakov {
  char make[256];
  
  /* for the reading function */
  input_record_t read[2];
  
} input_polyakov;

#define init_input_polyakov(varname) \
{ \
  .read={\
    {"make polyakov loops", "poly:make = %s", STRING_T, (varname).make},\
    {NULL, NULL, INT_T, NULL}\
  }\
}


input_polyakov poly_var = init_input_polyakov(poly_var);


/* Lowest-eigenvalue parameters */
typedef struct _input_eigval {
  char make[256];
  int nevt; /* search space dimension */
  int nev; /* number of accurate eigenvalues */
  int kmax; /* max degree of polynomial */
  int maxiter; /* max number of subiterations */
  double omega1; /* absolute precision */
  double omega2; /* relative precision */
  double evamass; /* mass to use in Dirac operator */
  
  /* for the reading function */
  input_record_t read[9];
  
} input_eigval;

#define init_input_eigval(varname) \
{ \
  .read={\
    {"make lowest eigenvalues", "eva:make = %s", STRING_T, (varname).make},\
    {"search space dimension", "eva:nevt = %d", INT_T, &(varname).nevt},\
    {"number of accurate eigenvalues", "eva:nev = %d", INT_T, &(varname).nev},\
    {"max degree of polynomial", "eva:kmax = %d", INT_T, &(varname).kmax},\
    {"max number of subiterations", "eva:maxiter = %d", INT_T, &(varname).maxiter},\
    {"absolute precision", "eva:omega1 = %lf", DOUBLE_T, &(varname).omega1},\
    {"relative precision", "eva:omega2 = %lf", DOUBLE_T, &(varname).omega2},\
    {"Dirac op mass", "eva:mass = %lf", DOUBLE_T, &(varname).evamass},\
    {NULL, NULL, INT_T, NULL}\
  }\
}

input_eigval eigval_var = init_input_eigval(eigval_var);



/* flow control variable */
hmc_flow flow=init_hmc_flow(flow);


char input_filename[256] = "input_file";
char output_filename[256] = "out_0";
char error_filename[256] = "err_0";
static void read_cmdline(int argc, char* argv[]) {
  int i, ai=0, ao=0, am=0, requested=1;

  for (i=1;i<argc;i++) {
    if (strcmp(argv[i],"-i")==0) {ai=i+1;requested+=2;}
    else if (strcmp(argv[i],"-o")==0) {ao=i+1;requested+=2;}
    else if (strcmp(argv[i],"-m")==0) {am=i;requested+=1;}
  }

  if (am != 0) {
    print_compiling_info();
    exit(0);
  }

  error(argc!=requested,1,"read_cmdline [hmc.c]",
      "Arguments: [-i <input file>] [-o <output file>] [-m]");

  if (ao!=0) strcpy(output_filename,argv[ao]);
  if (ai!=0) strcpy(input_filename,argv[ai]);
}


static void H2eva(spinor_field *out, spinor_field *in){
  g5Dphi_sq(eigval_var.evamass, out, in);
}



int main(int argc,char *argv[]) {
  int i,acc, rc;
  char sbuf[128];
 
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
  lprintf("MAIN",0,"[RepID: %d][world_size: %d]\n[MPI_ID: %d][MPI_size: %d]\n",RID,WORLD_SIZE,MPI_PID,MPI_WORLD_SIZE);
  lprintf("MAIN",0,"SVN Revision: %d\n", CI_svnrevision);

  //  lprintf("MAIN",0,"Logger lelvel: %d\n",logger_getlevel(0));
  
  /* setup lattice geometry */
  if (geometry_init() == 1) { finalize_process(); return 0; }
  geometry_mpi_eo();
  /* test_geometry_mpi_eo(); */

  /* setup random numbers */
  read_input(rlx_var.read,input_filename);
  lprintf("MAIN",0,"RLXD [%d,%d]\n",rlx_var.rlxd_level,rlx_var.rlxd_seed+MPI_PID);
  rlxd_init(rlx_var.rlxd_level,rlx_var.rlxd_seed+MPI_PID); /* use unique MPI_PID to shift seeds */

  if(strcmp(rlx_var.rlxd_start,"continue")==0 && rlx_var.rlxd_state[0]!='\0')
  {
    /*load saved state*/
    lprintf("MAIN",0,"Loading rlxd state from file [%s]\n",rlx_var.rlxd_state);
    read_ranlxd_state(rlx_var.rlxd_state);
  }

#ifdef GAUGE_SUN
  lprintf("MAIN",0,"Gauge group: SU(%d)\n",NG);
#elif GAUGE_SON
  lprintf("MAIN",0,"Gauge group: SO(%d)\n",NG);
#else
  lprintf("MAIN",0,"Default gauge group: SU(%d)\n",NG);
#endif
  lprintf("MAIN",0,"Fermion representation: " REPR_NAME " [dim=%d]\n",NF);

  /* read input for measures */
  read_input(mes_var.read,input_filename);
  read_input(poly_var.read,input_filename);
  read_input(eigval_var.read,input_filename);
  
  /* Init Monte Carlo */
  init_mc(&flow, input_filename);
  lprintf("MAIN",0,"MVM during HMC initialization: %ld\n",getMVM());
  lprintf("MAIN",0,"Initial plaquette: %1.8e\n",avr_plaquette());
 
  if(strcmp(mes_var.make,"true")==0) {
    init_meson_correlators(0);
    lprintf("MAIN",0,"Measuring Gamma Gamma correlators and PCAC-mass\n");
    lprintf("OBSERVABLES",0,"Inverter precision for mesons = %e\n",mes_var.precision);
    lprintf("OBSERVABLES",0,"Number of noisy sources for mesons per cnfg = %d\n",mes_var.nhits);
  }
  
  if(strcmp(eigval_var.make,"true")==0) {
    lprintf("OBSERVABLES",0,"EVA Search space dimension  (eva:nevt) = %d\n",eigval_var.nevt);
    lprintf("OBSERVABLES",0,"EVA Number of accurate eigenvalues (eva:nev) = %d\n",eigval_var.nev);
    lprintf("OBSERVABLES",0,"EVA Max degree of polynomial (eva:kmax) = %d\n",eigval_var.kmax);
    lprintf("OBSERVABLES",0,"EVA Max number of subiterations (eva:maxiter) = %d\n",eigval_var.maxiter);
    lprintf("OBSERVABLES",0,"EVA Absolute precision  (eva:omega1) = %e\n",eigval_var.omega1);
    lprintf("OBSERVABLES",0,"EVA Relative precision (eva:omega2) = %e\n",eigval_var.omega2);
  }

  
  double *eva_vals=NULL;
  spinor_field *eva_vecs=NULL;
  if(strcmp(eigval_var.make,"true")==0) {
    eva_vals=malloc(sizeof(double)*eigval_var.nevt);
    eva_vecs=alloc_spinor_field_f(eigval_var.nevt,&glattice);
  }
  rc=acc=0;
  for(i=flow.start;i<flow.end;++i) {
    int rr;
    double perc;
    struct timeval start, end, etime; /* //for trajectory timing */
    lprintf("MAIN",0,"Trajectory #%d...\n",i);
    
    gettimeofday(&start,0);
    
#ifdef MEASURE_FORCE
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
#endif
    
    rr=update_ghmc();

    gettimeofday(&end,0);
    timeval_subtract(&etime,&end,&start);
    lprintf("MAIN",0,"Trajectory #%d: generated in [%ld sec %ld usec]\n",i,etime.tv_sec,etime.tv_usec);
    
    if(rr<0) {
      lprintf("MAIN",0,"Error in updating the gauge field!!\n");
      return 1;
    } else if(rr!=0) {
      acc++;
    }
    rc++;
    perc=(acc==0)?0.:(float)(100*acc)/(float)(rc);
    
    lprintf("MAIN",0,"Trajectory #%d: %d/%d (%3.4f%%) MVM (f;d) = %ld ; %ld\n",i,acc,rc,perc,getMVM_flt(),getMVM());

    if((i%flow.save_freq)==0) {
      save_conf(&flow, i);
      if(u_scalar!=NULL){
	      save_scalar_conf(&flow, i);
      }
      /* Only save state if we have a file to save to */
      if(rlx_var.rlxd_state[0]!='\0') {
          lprintf("MAIN",0,"Saving rlxd state to file %s\n",rlx_var.rlxd_state);
          write_ranlxd_state(rlx_var.rlxd_state);
      }
    }
    
#ifdef MEASURE_FORCE
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
    for (int k=1;k<num_mon();++k){
      lprintf("INV_SUMMARY",0," Hasenbusch %d: %d",k-1,n_inv_iter[k]);
    }
    lprintf("INV_SUMMARY",0,"\n");
#endif

    if((i%flow.meas_freq)==0) {
      /* plaquette */
#ifdef WITH_SMEARING
		 lprintf("MAIN",0,"Plaquette: %1.8e, Smeared: %1.8e\n",avr_plaquette(),avr_smeared_plaquette());
#else
		 lprintf("MAIN",0,"Plaquette: %1.8e\n",avr_plaquette());
#endif

      /* Mesons */
      if(strcmp(mes_var.make,"true")==0) {
			measure_spectrum_semwall(1,&mes_var.mesmass,mes_var.nhits,i,mes_var.precision);
      }

      /* Four fermion observables */
      if(four_fermion_active==1) ff_observables();
      
      /* Polyakov loops */
      if(strcmp(poly_var.make,"true")==0) {
        polyakov();
      }
      
      /* Lowest eigenvalues */
      if(strcmp(eigval_var.make,"true")==0) {
        double max;
        max_H(&H2eva, &glattice, &max);
        max*=1.1;
        int status;
        int ie=eva(eigval_var.nev, eigval_var.nevt, 0, eigval_var.kmax, eigval_var.maxiter, max, eigval_var.omega1, eigval_var.omega2, &H2eva, eva_vecs, eva_vals, &status);
        while (ie!=0) { /* if failed restart EVA */
          lprintf("MAIN",0,"Restarting EVA!\n");
          ie=eva(eigval_var.nev, eigval_var.nevt, 2, eigval_var.kmax, eigval_var.maxiter, max, eigval_var.omega1, eigval_var.omega2, &H2eva, eva_vecs, eva_vals, &status);
        }
        
        for (int n=0;n<eigval_var.nev;++n) {
          lprintf("LOWEIG",0,"Eig %d = %1.15e\n",n,eva_vals[n]);
        }
      }
    }
  }

  /* save final configuration */
  if(((--i)%flow.save_freq)!=0) {
    save_conf(&flow, i);
    if(u_scalar!=NULL){
	    save_scalar_conf(&flow, i);
    }
    /* Only save state if we have a file to save to */
    if(rlx_var.rlxd_state[0]!='\0') {
        lprintf("MAIN",0,"Saving rlxd state to file %s\n",rlx_var.rlxd_state);
        write_ranlxd_state(rlx_var.rlxd_state);
    }
  }
  
#ifdef MEASURE_FORCE
  free(force_ave);
  free(force_max);
  free(n_inv_iter);
#endif
  
  /* finalize Monte Carlo */
  end_mc();
  
  /* close communications */
  finalize_process();
	
  
  return 0;
  
}
