/*******************************************************************************
*
* Computation of the mesonic spectrum
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
#include "spectrum.h"
#include "clover_tools.h"
#include "gaugefix.h"
#include "cinfo.c"




#if defined(ROTATED_SF) && defined(BASIC_SF)
#error This code does not work with the Schroedinger functional !!!
#endif

#ifdef FERMION_THETA
#error This code does not work with the fermion twisting !!!
#endif


char cnfg_filename[256]="";
char list_filename[256]="";
char input_filename[256] = "input_file";
char output_filename[256] = "mesons.out";
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

static void prepend_name(char* result, char* str1, char* prefix){
    // Find the last slash
    char * loc = strchr(str1,'/');
    int n;
    if (loc == NULL){ // No slashes at all
        n = -1;
    }
    else // At least 1 slash
    {
        do
        {
            n = (int)(loc - str1);
            loc = strchr(loc + 1, '/');
        } while (loc != NULL);
    }
    // Copy characters after the slash into tail variable
    char tail[50]; 
    strcpy(tail,str1+n+1);

    // Copy characters up to and including the slash
    strncpy(result,str1,n+1);
    result[n+1]='\0';

    //Add prefix and tail
    strcat(result,prefix);
    strcat(result,tail);
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
  lprintf("MAIN",0,"cnfg file [%s]\n",cnfg_filename); 

  GLB_T=fpars.t; GLB_X=fpars.x; GLB_Y=fpars.y; GLB_Z=fpars.z;
  error(fpars.type==UNKNOWN_CNFG,1,"measure_spectrum.c","Bad name for a configuration file");
  error(fpars.nc!=NG,1,"measure_spectrum.c","Bad NG");

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



  /* setup communication geometry */
  if (geometry_init() == 1) {
    finalize_process();
    return 0;
  }

  /* setup lattice geometry */
  geometry_mpi_eo();
  /* test_geometry_mpi_eo(); */

  init_BCs(NULL);

  /* alloc global fields */
  u_gauge=alloc_gfield(&glattice);
  u_scalar=alloc_scalar_field(&glattice);
#ifdef ALLOCATE_REPR_GAUGE_FIELD
  u_gauge_f=alloc_gfield_f(&glattice);
#endif

 
  list=NULL;
  if(strcmp(list_filename,"")!=0) {
    error((list=fopen(list_filename,"r"))==NULL,1,"main [mk_mesons.c]" ,
	"Failed to open list file\n");
  }


  i=0;

  while(++i) {
    struct timeval start, end, etime;

    if(list!=NULL)
      if(fscanf(list,"%s",cnfg_filename)==0 || feof(list)) break;


    lprintf("MAIN",0,"Configuration from %s\n", cnfg_filename);
    read_gauge_field(cnfg_filename);
    
    char scalar_cnfg_filename[256];
    prepend_name(scalar_cnfg_filename,cnfg_filename,"scalar_");
    lprintf("DEBUG",0,"%s %s \n", cnfg_filename, scalar_cnfg_filename);
    lprintf("MAIN",0,"Configuration from %s\n", scalar_cnfg_filename);
    read_scalar_field(scalar_cnfg_filename);
    
    represent_gauge_field();

    gettimeofday(&start,0);

    lprintf("TEST",0,"Plaquette before gauge fixing = %1.6f\n",avr_plaquette());
    lprintf("TEST",0,"SUS before gauge fixing = %1.6f %1.6f \n",average_SUS().re,average_SUS().im);

    for(int cont=0; cont<NG; cont++){
	    lprintf("TEST",0,"Average scalar[%d], before gauge fixing = %1.6f %1.6f \n", cont, average_S().c[cont].re, average_S().c[cont].im);
    }

    double act = scalar_gaugefix(10, //= 0, 1, 2, 3 for Coulomb gauge else Landau
	1.8,	//overrelax
	10000,	//maxit
	1e-10, //tolerance
	u_gauge, //gauge
	u_scalar //scalar
	);
    lprintf("TEST",0,"action  %1.6f\n",act);

    lprintf("TEST",0,"Plaquette after gauge fixing = %1.6f\n",avr_plaquette());
    lprintf("TEST",0,"SUS after gauge fixing = %1.6f %1.6f \n",average_SUS().re,average_SUS().im);

    for(int cont=0; cont<NG; cont++){
	    lprintf("TEST",0,"Average scalar[%d], after gauge fixing = %1.6f %1.6f \n", cont, average_S().c[cont].re, average_S().c[cont].im);
    }

    gettimeofday(&end,0);
    timeval_subtract(&etime,&end,&start);
    lprintf("MAIN",0,"Configuration #%d: analysed in [%ld sec %ld usec]\n",i,etime.tv_sec,etime.tv_usec);
    if(list==NULL) break;
  }

  if(list!=NULL) fclose(list);
  free_meson_observables();
  free_BCs();


  free_gfield(u_gauge);
  free_scalar_field(u_scalar);
#ifdef ALLOCATE_REPR_GAUGE_FIELD
  free_gfield_f(u_gauge_f);
#endif

  finalize_process();

  return 0;
}

