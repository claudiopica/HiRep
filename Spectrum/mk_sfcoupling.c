/*******************************************************************************
*
* Computation of the SF coupling
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

#include "cinfo.c"


/* we need the beta for normalization */
typedef struct _input_sfc {
  /* rhmc parameters */
  double beta;

  /* for the reading function */
  input_record_t read[2];
  
} input_sfc;

#define init_input_sfc(varname) \
{ \
  .read={\
    {"beta", "beta = %lf", DOUBLE_T, &(varname).beta},\
    {NULL, NULL, INT_T, NULL}\
  }\
}

input_sfc sfc_var = init_input_sfc( sfc_var );

/* BC variables */
typedef struct _input_bcpar {
  /* rhmc parameters */
  double theta[4];
  double SF_ct;
  double SF_ds;

  /* for the reading function */
  input_record_t read[7];
  
} input_bcpar;

#define init_input_bcpar(varname) \
{ \
  .read={\
    {"theta_T", "theta_T = %lf", DOUBLE_T, &(varname).theta[0]},\
    {"theta_X", "theta_X = %lf", DOUBLE_T, &(varname).theta[1]},\
    {"theta_Y", "theta_Y = %lf", DOUBLE_T, &(varname).theta[2]},\
    {"theta_Z", "theta_Z = %lf", DOUBLE_T, &(varname).theta[3]},\
    {"SF_ds", "SF_ds = %lf", DOUBLE_T, &(varname).SF_ds},\
    {"SF_ct", "SF_ct = %lf", DOUBLE_T, &(varname).SF_ct}, \
    {NULL, NULL, INT_T, NULL}\
  }\
}

input_bcpar bcpar_var = init_input_bcpar( bcpar_var );

char cnfg_filename[256]="";
char list_filename[256]="";
char input_filename[256] = "input_file";
char output_filename[256] = "sfcoupling.out";
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
  int i;
  char tmp[256];
  FILE* list;
  filename_t fpars;
  int nm;
  double m[256];

  /* setup process id and communications */
  read_cmdline(argc, argv);
  setup_process(&argc,&argv);

  /* logger setup */
  /* disable logger for MPI processes != 0 */
  logger_setlevel(0,30);
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


  read_input(glb_var.read,input_filename);
  read_input(sfc_var.read,input_filename);
  read_input(bcpar_var.read,input_filename);
  GLB_T=fpars.t; GLB_X=fpars.x; GLB_Y=fpars.y; GLB_Z=fpars.z;
  error(fpars.type==UNKNOWN_CNFG,1,"mk_sfcoupling.c","Bad name for a configuration file");
  error(fpars.nc!=NG,1,"mk_sfcoupling.c","Bad NG");


  lprintf("MAIN",0,"RLXD [%d,%d]\n",glb_var.rlxd_level,glb_var.rlxd_seed);
  rlxd_init(glb_var.rlxd_level,glb_var.rlxd_seed+PID);
  srand(glb_var.rlxd_seed+PID);

  lprintf("MAIN",0,"Gauge group: SU(%d)\n",NG);
  lprintf("MAIN",0,"Fermion representation: " REPR_NAME " [dim=%d]\n",NF);

  nm=0;
  if(fpars.type==DYNAMICAL_CNFG) {
    nm=1;
    m[0] = fpars.m;
  } else if(fpars.type==QUENCHED_CNFG) {
   lprintf("MAIN",0,"Quenched configurations NOT suported");
   finalize_process();
   return 1;
  }


  /* setup communication geometry */
  if (geometry_init() == 1) {
    finalize_process();
    return 1;
  }

  /* setup lattice geometry */
  geometry_mpi_eo();
  /* test_geometry_mpi_eo(); */

 /* initialize boundary conditions */
  BCs_pars_t BCs_pars = {
    .fermion_twisting_theta = {0.,0.,0.,0.},
    .gauge_boundary_improvement_cs = 1.,
    .gauge_boundary_improvement_ct = 1.,
    .chiSF_boundary_improvement_ds = 1.,
    .SF_BCs = 0
  };
#ifdef FERMION_THETA
  BCs_pars.fermion_twisting_theta[0] = bcpar_var.theta[0];
  BCs_pars.fermion_twisting_theta[1] = bcpar_var.theta[1];
  BCs_pars.fermion_twisting_theta[2] = bcpar_var.theta[2];
  BCs_pars.fermion_twisting_theta[3] = bcpar_var.theta[3];
#endif
#ifdef ROTATED_SF
  BCs_pars.gauge_boundary_improvement_ct = bcpar_var.SF_ct;
  BCs_pars.chiSF_boundary_improvement_ds = bcpar_var.SF_ds;
  BCs_pars.SF_BCs = 1;
#endif
#ifdef BASIC_SF
  BCs_pars.SF_BCs = 1;
#endif
  init_BCs(&BCs_pars);

  /* alloc global gauge fields */
  u_gauge=alloc_gfield(&glattice);
#ifdef ALLOCATE_REPR_GAUGE_FIELD
  u_gauge_f=alloc_gfield_f(&glattice);
#endif

  lprintf("MAIN",0,
          "beta = %.8f\n"
#ifdef ROTATED_SF
	  "rotatedSF ds = %.8f\n"
          "rotatedSF ct = %.8f\n"
#endif /* ROTATED_SF */
          ,sfc_var.beta
          ,bcpar_var.SF_ds
          ,bcpar_var.SF_ct);


  list=NULL;
  if(strcmp(list_filename,"")!=0) {
    error((list=fopen(list_filename,"r"))==NULL,1,"main [mk_sfcoupling.c]" ,
	"Failed to open list file\n");
  }


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

          double gsf=SF_action(sfc_var.beta);
          lprintf("SF_action",10,"gsf = %.10e\n",gsf);

	  if(list==NULL) break;
  }

  if(list!=NULL) fclose(list);

  finalize_process();

  /*
     free_spinor_field_f(pta_qprop[0]);
     free(pta_qprop);
     free(tricorr);
   */

  free_BCs();

  free_gfield(u_gauge);
#ifdef ALLOCATE_REPR_GAUGE_FIELD
  free_gfield_f(u_gauge_f);
#endif

  /* close communications */
  finalize_process();

  return 0;
}

