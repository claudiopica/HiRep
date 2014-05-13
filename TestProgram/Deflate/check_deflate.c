/*******************************************************************************
*
* Checks of propagator, spinmatrix and the sequential sources
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
#include "communications.h"
#include "gaugefix.h"
#include "spectrum.h"
#include "gamma_spinor.h"
#include "spin_matrix.h"
#include "propagator.h"

#include "cinfo.c"

#if defined(ROTATED_SF) && defined(BASIC_SF)
#error This code does not work with the Schroedinger functional !!!
#endif

#ifdef FERMION_THETA
#error This code does not work with the fermion twisting !!!
#endif

typedef enum {semwall_src,point_src,gfwall_src} source_type_t;

/* Mesons parameters */
typedef struct _input_mesons {
  char mstring[256];
  double precision;
  /* for the reading function */
  input_record_t read[4];
} input_mesons;

#define init_input_mesons(varname) \
{ \
  .read={\
    {"quark quenched masses", "mes:masses = %s", STRING_T, (varname).mstring}, \
    {"inverter precision", "mes:precision = %lf", DOUBLE_T, &(varname).precision},\
    {NULL, NULL, INT_T, NULL}				\
   }							\
}

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
  double lbnd; /* max wanted eval */
  /* for the reading function */
  input_record_t read[9];

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
    {"quark quenched mass", "eva:lbnd = %lf", DOUBLE_T, &(varname).lbnd}, \
    {NULL, NULL, INT_T, NULL}\
  }\
}

char cnfg_filename[256]="";
char list_filename[256]="";
char input_filename[256] = "input_file";
char output_filename[256] = "mesons.out";
enum { UNKNOWN_CNFG, DYNAMICAL_CNFG, QUENCHED_CNFG };

input_mesons mes_var = init_input_mesons(mes_var);
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

  error((ac==0 && al==0) || (ac!=0 && al!=0),1,"parse_cmdline [mk_mesons.c]",
      "Syntax: mk_mesons { -c <config file> | -l <list file> } [-i <input file>] [-o <output file>] [-m]");

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




int main(int argc,char *argv[]) {
  int i,k;
  char tmp[256], *cptr;
  FILE* list;
  filename_t fpars;
  int nm;
  double m[256];
  
  /* setup process id and communications */
  read_cmdline(argc, argv);
  setup_process(&argc,&argv);

  read_input(glb_var.read,input_filename);
  setup_replicas();


  /* logger setup */
  /* disable logger for MPI processes != 0 */
  logger_setlevel(0,10);
  if (PID!=0) { logger_disable(); }
  if (PID==0) { 
    sprintf(tmp,">%s",output_filename); logger_stdout(tmp);
    sprintf(tmp,"err_%d",PID); 
    if (!freopen(tmp,"w",stderr)) lprintf("MAIN",0,"Error out not open\n");
  }

  lprintf("MAIN",0,"Compiled with macros: %s\n",MACROS); 
  lprintf("MAIN",0,"PId =  %d [world_size: %d]\n\n",PID,WORLD_SIZE); 
  lprintf("MAIN",0,"input file [%s]\n",input_filename); 
  lprintf("MAIN",0,"output file [%s]\n",output_filename); 
  if (list_filename!=NULL) lprintf("MAIN",0,"list file [%s]\n",list_filename); 
  else lprintf("MAIN",0,"cnfg file [%s]\n",cnfg_filename); 


  /* read & broadcast parameters */
  parse_cnfg_filename(cnfg_filename,&fpars);


  read_input(mes_var.read,input_filename);
  GLB_T=fpars.t; GLB_X=fpars.x; GLB_Y=fpars.y; GLB_Z=fpars.z;
  error(fpars.type==UNKNOWN_CNFG,1,"check_propagator.c","Bad name for a configuration file");
  error(fpars.nc!=NG,1,"check_propagator.c","Bad NG");
  


  read_input(rlx_var.read,input_filename);
  lprintf("MAIN",0,"RLXD [%d,%d]\n",rlx_var.rlxd_level,rlx_var.rlxd_seed);
  rlxd_init(rlx_var.rlxd_level,rlx_var.rlxd_seed+MPI_PID);
  srand(rlx_var.rlxd_seed+MPI_PID);

#ifdef GAUGE_SON
  lprintf("MAIN",0,"Gauge group: SO(%d)\n",NG);
#else
  lprintf("MAIN",0,"Gauge group: SU(%d)\n",NG);
#endif
  lprintf("MAIN",0,"Fermion representation: " REPR_NAME " [dim=%d]\n",NF);

  nm=0;
  if(fpars.type==DYNAMICAL_CNFG) {
    nm=1;
    m[0] = fpars.m;
  } else if(fpars.type==QUENCHED_CNFG) {
    strcpy(tmp,mes_var.mstring);
    cptr = strtok(tmp, ";");
    nm=0;
    while(cptr != NULL) {
      m[nm]=atof(cptr);
      nm++;
      cptr = strtok(NULL, ";");
    }            
  }

  read_input(eig_var.read,input_filename);

  /* setup communication geometry */
  if (geometry_init() == 1) {
    finalize_process();
    return 0;
  }

  /* setup lattice geometry */
  geometry_mpi_eo();
  /* test_geometry_mpi_eo(); */

  init_BCs(NULL);

  /* alloc global gauge fields */
  u_gauge=alloc_gfield(&glattice);
#ifdef ALLOCATE_REPR_GAUGE_FIELD
  u_gauge_f=alloc_gfield_f(&glattice);
#endif

  lprintf("MAIN",0,"Inverter precision = %e\n",mes_var.precision);
  for(k=0;k<nm;k++)
    lprintf("MAIN",0,"Mass[%d] = %f\n",k,m[k]);
  

  list=NULL;
  if(strcmp(list_filename,"")!=0) {
    error((list=fopen(list_filename,"r"))==NULL,1,"main [mk_mesons.c]" ,
	"Failed to open list file\n");
  }

  lprintf("MAIN",0,"EVA Parameters:\n");
  lprintf("MAIN",0,"search space dimension  (eva:nevt) = %d\n",eig_var.nevt);
  lprintf("MAIN",0,"number of accurate eigenvalues (eva:nev) = %d\n",eig_var.nev);
  lprintf("MAIN",0,"max degree of polynomial (eva:kmax) = %d\n",eig_var.kmax);
  lprintf("MAIN",0,"max number of subiterations (eva:maxiter) = %d\n",eig_var.maxiter);
  lprintf("MAIN",0,"absolute precision  (eva:omega1) = %e\n",eig_var.omega1);
  lprintf("MAIN",0,"relative precision (eva:omega2) = %e\n",eig_var.omega2);




  i=0;

  spinor_field* source;
  spinor_field* prop_1;
  spinor_field* prop_2;

  source = alloc_spinor_field_f(4,&glattice);
  prop_1 = alloc_spinor_field_f(4,&glattice);
  prop_2 = alloc_spinor_field_f(4,&glattice);

  while(1) {
    struct timeval start, end, etime, e2time;
    struct timeval eig_start, eig_end, eig_time;


    if(list!=NULL)
      if(fscanf(list,"%s",cnfg_filename)==0 || feof(list)) break;

    i++;

    lprintf("MAIN",0,"Configuration from %s\n", cnfg_filename);
    read_gauge_field(cnfg_filename);
    //unit_gauge(u_gauge);
    represent_gauge_field();

    lprintf("TEST",0,"<p> %1.6f\n",avr_plaquette());

    full_plaquette();



  int k, beta;
  init_propagator_eo(1,m, mes_var.precision);//1 for number of masses 



  //calc_propagator_eo(prop_1,source,1);//4 for spin components

  lprintf("SOLVE",10,"\n");


    gettimeofday(&eig_start,0);
    eig_init(eig_var.nev,eig_var.nev+10,eig_var.kmax,eig_var.maxiter,eig_var.lbnd,eig_var.omega1,eig_var.omega2);
  //eig_init(eig_var.nev,eig_var.nevt,eig_var.kmax,eig_var.maxiter,eig_var.omega1,eig_var.omega2);
    gettimeofday(&eig_end,0);
    timeval_subtract(&eig_time,&eig_end,&eig_start);

    lprintf("MAIN",0,"Eigenvectors found in [%lf sec ]\n",(double)eig_time.tv_sec + (double)eig_time.tv_usec*1e-6);

    gettimeofday(&start,0);
    create_point_source(source,0,0);
    calc_deflated_propagator(prop_2, source, 1, 0);
    gettimeofday(&end,0);
    timeval_subtract(&etime,&end,&start);
    lprintf("MAIN",0,"Sfull = %g\n",(double)etime.tv_sec + (double)etime.tv_usec*1e-6 );

    gettimeofday(&start,0);
    create_point_source(source,0,0);
    calc_deflated_propagator(prop_2, source, 1, eig_var.nev);
    gettimeofday(&end,0);
    timeval_subtract(&e2time,&end,&start);
    lprintf("MAIN",0,"Sn = %g\n", (double)e2time.tv_sec + (double)e2time.tv_usec*1e-6 );

    timeval_subtract(&end,&etime,&e2time);
    lprintf("MAIN",0,"saved [%ld sec %ld usec]\n",end.tv_sec,end.tv_usec);


    lprintf("MAIN",0,"m %g n %d prec %g Saved %g\n", 
	m[0], eig_var.nev, eig_var.omega1, 
( 4*NF*GLB_T * ( (double)etime.tv_sec + (double)etime.tv_usec*1e-6 ) ) / 
( 4*NF*GLB_T * ( (double)e2time.tv_sec + (double)e2time.tv_usec*1e-6 ) + (double)eig_time.tv_sec + (double)eig_time.tv_usec*1e-6 ) );

  free_propagator_eo(); 



    if(list==NULL) break;
  }


  free_spinor_field_f(source);  
  free_spinor_field_f(prop_1);  
  free_spinor_field_f(prop_2);  

  if(list!=NULL) fclose(list);

  free_BCs();

  free_gfield(u_gauge);
#ifdef ALLOCATE_REPR_GAUGE_FIELD
  free_gfield_f(u_gauge_f);
#endif

  finalize_process();

  return 0;
}

