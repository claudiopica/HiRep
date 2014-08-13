/***************************************************************************\
 * Copyright (c) 2009, Claudio Pica                                          *   
 * All rights reserved.                                                      * 
 \***************************************************************************/

/*******************************************************************************
*
* Computation of the lowest eigenvalues of H^2
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

#include "cinfo.c"

#if defined(ROTATED_SF) && defined(BASIC_SF)
#error This code does not work with the Schroedinger functional !!!
#endif

#ifdef FERMION_THETA
#error This code does not work with the fermion twisting !!!
#endif

/* Mesons parameters */
typedef struct _input_eigval {
  /* EVA parameters */
  int nevt; /* search space dimension */
  int nev; /* number of accurate eigenvalues */
  int kmax; /* max degree of polynomial */
  int maxiter; /* max number of subiterations */
  double omega1; /* absolute precision */
  double omega2; /* relative precision */
  double mass; /* quenched mass */

  /* for the reading function */
  input_record_t read[8];

} input_eigval;

#define init_input_eigval(varname) \
{ \
  .read={\
    {"search space dimension", "eva:nevt = %d", INT_T, &(varname).nevt},\
    {"number of accurate eigenvalues", "eva:nev = %d", INT_T, &(varname).nev},\
    {"max degree of polynomial", "eva:kmax = %d", INT_T, &(varname).kmax},\
    {"max number of subiterations", "eva:maxiter = %d", INT_T, &(varname).maxiter},\
    {"absolute precision", "eva:omega1 = %lf", DOUBLE_T, &(varname).omega1},\
    {"relative precision", "eva:omega2 = %lf", DOUBLE_T, &(varname).omega2},\
    {"quark quenched mass", "eva:mass = %lf", DOUBLE_T, &(varname).mass}, \
    {NULL, NULL, INT_T, NULL}\
  }\
}


char cnfg_filename[256]="";
char list_filename[256]="";
char input_filename[256] = "input_file";
char output_filename[256] = "eigval.out";
enum { UNKNOWN_CNFG, DYNAMICAL_CNFG, QUENCHED_CNFG };

input_eigval eig_var = init_input_eigval(eig_var);


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

  hm=sscanf(basename,"%*[^_]_%dx%dx%dx%d%*[Nn]c%db%lfn%d",
      &(fn->t),&(fn->x),&(fn->y),&(fn->z),&(fn->nc),&(fn->b),&(fn->n));
  if(hm==7) {
    fn->type=QUENCHED_CNFG;
    return QUENCHED_CNFG;
  }

  hm=sscanf(basename,"%dx%dx%dx%d%*[Nn]c%db%lfn%d",
      &(fn->t),&(fn->x),&(fn->y),&(fn->z),&(fn->nc),&(fn->b),&(fn->n));
  if(hm==7) {
    fn->type=QUENCHED_CNFG;
    return QUENCHED_CNFG;
  }

  fn->type=UNKNOWN_CNFG;
  return UNKNOWN_CNFG;
}


void read_cmdline(int argc, char* argv[]) {
  int i, ai=0, ao=0, ac=0, al=0, am=0;
  FILE *list=NULL;

  for (i=1;i<argc;i++) {
    if (strcmp(argv[i],"-i")==0) ai=i+1;
    else if (strcmp(argv[i],"-o")==0) ao=i+1;
    else if (strcmp(argv[i],"-c")==0) ac=i+1;
    else if (strcmp(argv[i],"-l")==0) al=i+1;
    else if (strcmp(argv[i],"-m")==0) am=i;
  }

  if (am != 0) {
    print_compiling_info();
    exit(0);
  }

  if (ao!=0) strcpy(output_filename,argv[ao]);
  if (ai!=0) strcpy(input_filename,argv[ai]);

  error((ac==0 && al==0) || (ac!=0 && al!=0),1,"parse_cmdline [eigval.c]",
      "Syntax: eigval { -c <config file> | -l <list file> } [-i <input file>] [-o <output file>] [-m]");

  if(ac != 0) {
    strcpy(cnfg_filename,argv[ac]);
    strcpy(list_filename,"");
  } else if(al != 0) {
    strcpy(list_filename,argv[al]);
    error((list=fopen(list_filename,"r"))==NULL,1,"parse_cmdline [eigval.c]" ,
	"Failed to open list file\n");
    error(fscanf(list,"%s",cnfg_filename)==0,1,"parse_cmdline [eigval.c]" ,
	"Empty list file\n");
    fclose(list);
  }


}

double hevamass=0.;
void H2EVA(spinor_field *out, spinor_field *in){
  g5Dphi_sq(hevamass, out, in);
}
void HEVA(spinor_field *out, spinor_field *in){
  g5Dphi(hevamass, out, in);
}

int main(int argc,char *argv[]) {
  int i,n;
  char tmp[256];
  FILE* list;
  filename_t fpars;

  /* setup process id and communications */
  read_cmdline(argc, argv);
  setup_process(&argc,&argv);

  read_input(glb_var.read,input_filename);
  setup_replicas();

  /* logger setup */
  /* disable logger for MPI processes != 0 */
  logger_setlevel(0,50);
  if (PID!=0) { logger_disable(); }
  if (PID==0) { 
    sprintf(tmp,">%s",output_filename); logger_stdout(tmp);
    sprintf(tmp,"err_%d",PID); freopen(tmp,"w",stderr);
  }

  lprintf("MAIN",0,"Compiled with macros: %s\n",MACROS); 
  lprintf("MAIN",0,"PId =  %d [world_size: %d]\n\n",PID,WORLD_SIZE); 
  lprintf("MAIN",0,"input file [%s]\n",input_filename); 
  lprintf("MAIN",0,"output file [%s]\n",output_filename); 
  if (list_filename!=NULL) lprintf("MAIN",0,"list file [%s]\n",list_filename); 
  else lprintf("MAIN",0,"cnfg file [%s]\n",cnfg_filename); 


  /* read & broadcast parameters */
  parse_cnfg_filename(cnfg_filename,&fpars);

  eig_var.mass=-20.;
  read_input(eig_var.read,input_filename);
  GLB_T=fpars.t; GLB_X=fpars.x; GLB_Y=fpars.y; GLB_Z=fpars.z;
  error(fpars.type==UNKNOWN_CNFG,1,"eigval.c","Bad name for a configuration file");
  error(fpars.nc!=NG,1,"eigval.c","Bad NG");

  lprintf("MAIN",0,"Gauge group: SU(%d)\n",NG);
  lprintf("MAIN",0,"Fermion representation: " REPR_NAME " [dim=%d]\n",NF);

  if(fpars.type==DYNAMICAL_CNFG)
    hevamass = fpars.m;
  if(fpars.type==QUENCHED_CNFG || eig_var.mass>-10.)
    hevamass = eig_var.mass;


  /* setup communication geometry */
  if (geometry_init() == 1) {
    finalize_process();
    return 0;
  }

  /* setup lattice geometry */
  geometry_mpi_eo();
  /* test_geometry_mpi_eo(); */
    
    /* setup random numbers */
    read_input(rlx_var.read,input_filename);
    lprintf("MAIN",0,"RLXD [%d,%d]\n",rlx_var.rlxd_level,rlx_var.rlxd_seed+MPI_PID);
    rlxd_init(rlx_var.rlxd_level,rlx_var.rlxd_seed+MPI_PID); /* use unique MPI_PID to shift seeds */


  init_BCs(NULL);

  /* alloc global gauge fields */
  u_gauge=alloc_gfield(&glattice);
#ifdef ALLOCATE_REPR_GAUGE_FIELD
  u_gauge_f=alloc_gfield_f(&glattice);
#endif


  lprintf("MAIN",0,"EVA Parameters:\n");
  lprintf("MAIN",0,"search space dimension  (eva:nevt) = %d\n",eig_var.nevt);
  lprintf("MAIN",0,"number of accurate eigenvalues (eva:nev) = %d\n",eig_var.nev);
  lprintf("MAIN",0,"max degree of polynomial (eva:kmax) = %d\n",eig_var.kmax);
  lprintf("MAIN",0,"max number of subiterations (eva:maxiter) = %d\n",eig_var.maxiter);
  lprintf("MAIN",0,"absolute precision  (eva:omega1) = %e\n",eig_var.omega1);
  lprintf("MAIN",0,"relative precision (eva:omega2) = %e\n",eig_var.omega2);
  lprintf("MAIN",0,"mass = %f\n",hevamass);


  list=NULL;
  if(strcmp(list_filename,"")!=0) {
    error((list=fopen(list_filename,"r"))==NULL,1,"main [eigval.c]" ,
	"Failed to open list file\n");
  }

  
  /* EVA parameters */
  double max, mupp;
  double *eva_val;
  int status,ie;
  spinor_field *eva_vec, *eva_ws;

  eva_val=malloc(sizeof(double)*eig_var.nevt);
  eva_vec=alloc_spinor_field_f(eig_var.nevt+1,&glattice);
  eva_ws=eva_vec+eig_var.nevt;

  mupp=fabs(hevamass+4)+4;
  mupp*=mupp;
  /* END of EVA parameters */

  i=0;
  while(1) {

    if(list!=NULL)
      if(fscanf(list,"%s",cnfg_filename)==0 || feof(list)) break;

    i++;

    lprintf("MAIN",0,"Configuration from %s\n", cnfg_filename);
    /* NESSUN CHECK SULLA CONSISTENZA CON I PARAMETRI DEFINITI !!! */
    read_gauge_field(cnfg_filename);
    represent_gauge_field();

    lprintf("TEST",0,"<p> %1.6f\n",avr_plaquette());


    int MVM=0; /* counter for matrix-vector multiplications */

    max_H(&H2EVA, &glattice, &max);
    lprintf("MAIN",0,"MAXCHECK: cnfg=%e  uppbound=%e diff=%e %s\n",max,mupp,mupp-max,(mupp-max)<0?"[FAILED]":"[OK]");
    max*=1.1;

    ie=eva(eig_var.nev,eig_var.nevt,0,eig_var.kmax,eig_var.maxiter,max,eig_var.omega1,eig_var.omega2,&H2EVA,eva_vec,eva_val,&status);
    MVM+=status;
    while (ie!=0) { /* if failed restart EVA */
      lprintf("MAIN",0,"Restarting EVA!\n");
      ie=eva(eig_var.nev,eig_var.nevt,2,eig_var.kmax,eig_var.maxiter,max,eig_var.omega1,eig_var.omega2,&H2EVA,eva_vec,eva_val,&status);
      MVM+=status;
    }

    lprintf("MAIN",0,"EVA MVM = %d\n",MVM);
    for (n=0;n<eig_var.nev;++n) {
      HEVA(&eva_ws[0],&eva_vec[n]);
      lprintf("RESULT",0,"Eig %d = %.15e %.15e\n",n,eva_val[n],
        spinor_field_prod_re_f(&eva_ws[0],&eva_vec[n])/spinor_field_sqnorm_f(&eva_vec[n]));
    }
    
    if(list==NULL) break;
  }

  if(list!=NULL) fclose(list);

  free_BCs();

  free(eva_val);
  free_spinor_field_f(eva_vec);

  free_gfield(u_gauge);
#ifdef ALLOCATE_REPR_GAUGE_FIELD
  free_gfield_f(u_gauge_f);
#endif

  finalize_process();

  return 0;
}

