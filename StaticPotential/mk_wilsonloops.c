/*******************************************************************************
*
* Computation of the Wilson loops for the static potential
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
#error This code does not work with the Schroedinger functional
#endif

#ifdef BC_XYZ_TWISTED
#error This code does not work with the twisted BCs
#endif

#ifdef BC_T_OPEN
#error This code does not work with the open BCs
#endif

/* HYP smearing parameters */
typedef struct _input_HYP {
/*  int nsteps;*/
  double weight[3];

  /* for the reading function */
  input_record_t read[4];

} input_HYP;

#define init_input_HYP(varname) \
{ \
  .read={\
    {"HYP smearing weight[0]", "HYP:weight0 = %lf", DOUBLE_T, &((varname).weight[0])},\
    {"HYP smearing weight[1]", "HYP:weight1 = %lf", DOUBLE_T, &((varname).weight[1])},\
    {"HYP smearing weight[2]", "HYP:weight2 = %lf", DOUBLE_T, &((varname).weight[2])},\
    {NULL, NULL, INT_T, NULL}\
  }\
}

typedef struct _input_wilson {
  int c[10][3];
  int nsteps[10];

  /* for the reading function */
  input_record_t read[41];

} input_WL;

#define init_input_WL(varname) \
{ \
  .read={\
    {"WL load path[0] delta(x)", "WL[0]:delta.x = %d", INT_T, &((varname).c[0][0])},\
    {"WL load path[0] delta(y)", "WL[0]:delta.y = %d", INT_T, &((varname).c[0][1])},\
    {"WL load path[0] delta(z)", "WL[0]:delta.z = %d", INT_T, &((varname).c[0][2])},\
    {"WL load path[0] nsteps", "WL[0]:nsteps = %d", INT_T, &((varname).nsteps[0])},\
    {"WL load path[1] delta(x)", "WL[1]:delta.x = %d", INT_T, &((varname).c[1][0])},\
    {"WL load path[1] delta(y)", "WL[1]:delta.y = %d", INT_T, &((varname).c[1][1])},\
    {"WL load path[1] delta(z)", "WL[1]:delta.z = %d", INT_T, &((varname).c[1][2])},\
    {"WL load path[1] nsteps", "WL[1]:nsteps = %d", INT_T, &((varname).nsteps[1])},\
    {"WL load path[2] delta(x)", "WL[2]:delta.x = %d", INT_T, &((varname).c[2][0])},\
    {"WL load path[2] delta(y)", "WL[2]:delta.y = %d", INT_T, &((varname).c[2][1])},\
    {"WL load path[2] delta(z)", "WL[2]:delta.z = %d", INT_T, &((varname).c[2][2])},\
    {"WL load path[2] nsteps", "WL[2]:nsteps = %d", INT_T, &((varname).nsteps[2])},\
    {"WL load path[3] delta(x)", "WL[3]:delta.x = %d", INT_T, &((varname).c[3][0])},\
    {"WL load path[3] delta(y)", "WL[3]:delta.y = %d", INT_T, &((varname).c[3][1])},\
    {"WL load path[3] delta(z)", "WL[3]:delta.z = %d", INT_T, &((varname).c[3][2])},\
    {"WL load path[3] nsteps", "WL[3]:nsteps = %d", INT_T, &((varname).nsteps[3])},\
    {"WL load path[4] delta(x)", "WL[4]:delta.x = %d", INT_T, &((varname).c[4][0])},\
    {"WL load path[4] delta(y)", "WL[4]:delta.y = %d", INT_T, &((varname).c[4][1])},\
    {"WL load path[4] delta(z)", "WL[4]:delta.z = %d", INT_T, &((varname).c[4][2])},\
    {"WL load path[4] nsteps", "WL[4]:nsteps = %d", INT_T, &((varname).nsteps[4])},\
    {"WL load path[5] delta(x)", "WL[5]:delta.x = %d", INT_T, &((varname).c[5][0])},\
    {"WL load path[5] delta(y)", "WL[5]:delta.y = %d", INT_T, &((varname).c[5][1])},\
    {"WL load path[5] delta(z)", "WL[5]:delta.z = %d", INT_T, &((varname).c[5][2])},\
    {"WL load path[5] nsteps", "WL[5]:nsteps = %d", INT_T, &((varname).nsteps[5])},\
    {"WL load path[6] delta(x)", "WL[6]:delta.x = %d", INT_T, &((varname).c[6][0])},\
    {"WL load path[6] delta(y)", "WL[6]:delta.y = %d", INT_T, &((varname).c[6][1])},\
    {"WL load path[6] delta(z)", "WL[6]:delta.z = %d", INT_T, &((varname).c[6][2])},\
    {"WL load path[6] nsteps", "WL[6]:nsteps = %d", INT_T, &((varname).nsteps[6])},\
    {"WL load path[7] delta(x)", "WL[7]:delta.x = %d", INT_T, &((varname).c[7][0])},\
    {"WL load path[7] delta(y)", "WL[7]:delta.y = %d", INT_T, &((varname).c[7][1])},\
    {"WL load path[7] delta(z)", "WL[7]:delta.z = %d", INT_T, &((varname).c[7][2])},\
    {"WL load path[7] nsteps", "WL[7]:nsteps = %d", INT_T, &((varname).nsteps[7])},\
    {"WL load path[8] delta(x)", "WL[8]:delta.x = %d", INT_T, &((varname).c[8][0])},\
    {"WL load path[8] delta(y)", "WL[8]:delta.y = %d", INT_T, &((varname).c[8][1])},\
    {"WL load path[8] delta(z)", "WL[8]:delta.z = %d", INT_T, &((varname).c[8][2])},\
    {"WL load path[8] nsteps", "WL[8]:nsteps = %d", INT_T, &((varname).nsteps[8])},\
    {"WL load path[9] delta(x)", "WL[9]:delta.x = %d", INT_T, &((varname).c[9][0])},\
    {"WL load path[9] delta(y)", "WL[9]:delta.y = %d", INT_T, &((varname).c[9][1])},\
    {"WL load path[9] delta(z)", "WL[9]:delta.z = %d", INT_T, &((varname).c[9][2])},\
    {"WL load path[9] nsteps", "WL[9]:nsteps = %d", INT_T, &((varname).nsteps[9])},\
    {NULL, NULL, INT_T, NULL}\
  }\
}


input_HYP HYP_var = init_input_HYP(HYP_var);
input_WL WL_var = init_input_WL(WL_var);

char cnfg_filename[256]="";
char list_filename[256]="";
char input_filename[256] = "input_file";
char output_filename[256] = "wilson.out";
enum { UNKNOWN_CNFG=0, DYNAMICAL_CNFG, QUENCHED_CNFG };


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

/*#ifdef REPR_FUNDAMENTAL*/
/*#define repr_name "FUN"*/
/*#elif defined REPR_SYMMETRIC*/
/*#define repr_name "SYM"*/
/*#elif defined REPR_ANTISYMMETRIC*/
/*#define repr_name "ASY"*/
/*#elif defined REPR_ADJOINT*/
/*#define repr_name "ADJ"*/
/*#endif*/
  hm=sscanf(basename,"%*[^_]_%dx%dx%dx%d%*[Nn]c%dr%*[FSA]%*[UYSD]%*[NMYJ]%*[Nn]f%db%lfm%lfn%d",
      &(fn->t),&(fn->x),&(fn->y),&(fn->z),&(fn->nc),&(fn->nf),&(fn->b),&(fn->m),&(fn->n));
  if(hm==9) {
    fn->m=-fn->m; /* invert sign of mass */
    fn->type=DYNAMICAL_CNFG;
    return DYNAMICAL_CNFG;
  }
/*#undef repr_name*/

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

  error((ac==0 && al==0) || (ac!=0 && al!=0),1,"parse_cmdline [mk_wilsonloops.c]",
      "Syntax: mk_wilsonloops { -c <config file> | -l <list file> } [-i <input file>] [-o <output file>] [-m]");

  if(ac != 0) {
    strcpy(cnfg_filename,argv[ac]);
    strcpy(list_filename,"");
  } else if(al != 0) {
    strcpy(list_filename,argv[al]);
    error((list=fopen(list_filename,"r"))==NULL,1,"parse_cmdline [mk_wilsonloops.c]" ,
	"Failed to open list file\n");
    error(fscanf(list,"%s",cnfg_filename)==0,1,"parse_cmdline [mk_wilsonloops.c]" ,
	"Empty list file\n");
    fclose(list);
  }


}


int main(int argc,char *argv[]) {
  int i;
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
  logger_setlevel(0,70);
  if (PID!=0) { logger_disable(); }
  if (PID==0) { 
    sprintf(tmp,">%s",output_filename); logger_stdout(tmp);
    sprintf(tmp,"err_%d",PID); freopen(tmp,"w",stderr);
  }

  lprintf("MAIN",0,"Compiled with macros: %s\n",MACROS); 
  lprintf("MAIN",0,"PId =  %d [world_size: %d]\n\n",PID,WORLD_SIZE); 
  lprintf("MAIN",0,"input file [%s]\n",input_filename); 
  lprintf("MAIN",0,"output file [%s]\n",output_filename); 
  if (strcmp(list_filename,"")!=0) lprintf("MAIN",0,"list file [%s]\n",list_filename); 
  else lprintf("MAIN",0,"cnfg file [%s]\n",cnfg_filename); 


  /* read & broadcast parameters */
  parse_cnfg_filename(cnfg_filename,&fpars);
  
  HYP_var.weight[0]=HYP_var.weight[1]=HYP_var.weight[2]=0.;
  for(i=0;i<10;i++) { WL_var.c[i][0]=WL_var.c[i][1]=WL_var.c[i][2]=WL_var.nsteps[i]=0; }
  read_input(HYP_var.read,input_filename);
  read_input(WL_var.read,input_filename);
  GLB_T=fpars.t; GLB_X=fpars.x; GLB_Y=fpars.y; GLB_Z=fpars.z;
 
  error(fpars.type==UNKNOWN_CNFG,1,"mk_wilsonloops.c","Bad name for a configuration file");
  error(fpars.nc!=NG,1,"mk_wilsonloops.c","Bad NG");


  lprintf("MAIN",0,"Gauge group: SU(%d)\n",NG);
  lprintf("MAIN",0,"Fermion representation: " REPR_NAME " [dim=%d]\n",NF);

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

  lprintf("MAIN",0,"HYP smearing weights: %f %f %f\n",HYP_var.weight[0],HYP_var.weight[1],HYP_var.weight[2]);
/*  lprintf("MAIN",0,"HYP smearing number of steps: %d\n",HYP_var.nsteps);*/

  /* alloc global gauge fields */
  u_gauge=alloc_gfield(&glattice);
  
  list=NULL;
  if(strcmp(list_filename,"")!=0) {
    error((list=fopen(list_filename,"r"))==NULL,1,"main [mk_mesons.c]" ,
	"Failed to open list file\n");
  }

  WL_initialize();
  for(i=0;i<10;i++){
    if(WL_var.c[i][0]*WL_var.c[i][0]+WL_var.c[i][1]*WL_var.c[i][1]+WL_var.c[i][2]*WL_var.c[i][2]!=0 && WL_var.nsteps[i]!=0)
      WL_load_path(WL_var.c[i], WL_var.nsteps[i]);
  }
  
  i=0;
  while(1) {

    if(list!=NULL)
      if(fscanf(list,"%s",cnfg_filename)==0 || feof(list)) break;

    i++;

    lprintf("MAIN",0,"Configuration from %s\n", cnfg_filename);
    /* NESSUN CHECK SULLA CONSISTENZA CON I PARAMETRI DEFINITI !!! */
    read_gauge_field(cnfg_filename);

    lprintf("TEST",0,"<p> %1.6f\n",avr_plaquette());

    full_plaquette();

    WL_wilsonloops(HYP_var.weight);
   
    if(list==NULL) break;
  }

  if(list!=NULL) fclose(list);

  free_BCs();
 
  free_gfield(u_gauge);

  finalize_process();
  
  return 0;
}

