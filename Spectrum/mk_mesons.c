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

#include "cinfo.c"

#if defined(ROTATED_SF) && defined(BASIC_SF)
#error This code does not work with the Schroedinger functional !!!
#endif

#ifdef FERMION_THETA
#error This code does not work with the fermion twisting !!!
#endif



/* Mesons parameters */
typedef struct _input_mesons {
  char mstring[256];
  double precision;
  int nhits;

  /* for the reading function */
  input_record_t read[4];

} input_mesons;

#define init_input_mesons(varname) \
{ \
  .read={\
    {"quark quenched masses", "mes:masses = %s", STRING_T, (varname).mstring},\
    {"inverter precision", "mes:precision = %lf", DOUBLE_T, &(varname).precision},\
    {"number of inversions per cnfg", "mes:nhits = %d", INT_T, &(varname).nhits},\
    {NULL, NULL, INT_T, NULL}\
  }\
}


char cnfg_filename[256]="";
char list_filename[256]="";
char input_filename[256] = "input_file";
char output_filename[256] = "mesons.out";
enum { UNKNOWN_CNFG, DYNAMICAL_CNFG, QUENCHED_CNFG };

input_mesons mes_var = init_input_mesons(mes_var);

void inline_mk_mesons(double *m, int nm, double prec) {
    int k, n, g0[4];
    spinor_field **pta_qprop=0;
    double* tricorr;

    tricorr=(double*)malloc(GLB_T*sizeof(double));
    pta_qprop=(spinor_field**)malloc(sizeof(spinor_field*)*nm);
    pta_qprop[0]=alloc_spinor_field_f(4*NF*nm,&glattice);

    for(k=0;k<nm;++k)
	    pta_qprop[k]=pta_qprop[0]+4*NF*k;

    g0[0]=rand()%GLB_T; g0[1]=rand()%GLB_X; g0[2]=rand()%GLB_Y; g0[3]=rand()%GLB_Z;
    if((g0[0]+g0[1]+g0[2]+g0[3])%2!=0)
	    g0[3]=(g0[3]+1)%GLB_Z;

    bcast_int(g0,4);

    lprintf("MAIN",0,"PTA meson source in (%d,%d,%d,%d)\n",g0[0],g0[1],g0[2],g0[3]);

    pta_qprop_QMR_eo(g0, pta_qprop, nm, m, prec);

    for (k=0;k<nm;++k){

#define CORR(name) \
	    name##_correlator(tricorr, g0[0], pta_qprop[k]);\
	    lprintf("MAIN",0,"conf #0 mass=%2.6f TRIPLET " #name "= ",m[k]);\
	    for(n=0;n<GLB_T;++n) {\
		    lprintf("MAIN",0,"%e ",tricorr[n]);\
	    }\
	    lprintf("MAIN",0,"\n");\
	    fflush(stdout)

	    CORR(id);
	    CORR(g5);
	    CORR(g0);
	    CORR(g0g5);
	    CORR(g1);
	    CORR(g2);
	    CORR(g3);
	    CORR(g0g1);
	    CORR(g0g2);
	    CORR(g0g3);
	    CORR(g5g1);
	    CORR(g5g2);
	    CORR(g5g3);
	    CORR(g0g5g1);
	    CORR(g0g5g2);
	    CORR(g0g5g3);
	    CORR(g5_g0g5_re);
	    CORR(g5_g0g5_im);

    }
#undef CORR

  free_spinor_field_f(pta_qprop[0]);
  free(pta_qprop);
  free(tricorr);

}

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
  read_input(rlx_var.read,input_filename);
  setup_replicas();

  /* logger setup */
  /* disable logger for MPI processes != 0 */
  logger_setlevel(0,30);
  if (PID!=0) { logger_disable(); }
  if (PID==0) { 
    sprintf(tmp,">%s",output_filename); logger_stdout(tmp);
    sprintf(tmp,"err_%d",PID); 
    if(!freopen(tmp,"w",stderr)){
      error(1,1,"mk_mesons.c","Cannot open something!\n");
    }
  }

  lprintf("MAIN",0,"Compiled with macros: %s\n",MACROS); 
  lprintf("MAIN",0,"PId =  %d [world_size: %d]\n\n",PID,WORLD_SIZE); 
  lprintf("MAIN",0,"input file [%s]\n",input_filename); 
  lprintf("MAIN",0,"output file [%s]\n",output_filename); 
  if (list_filename!=NULL) lprintf("MAIN",0,"list file [%s]\n",list_filename); 
  else lprintf("MAIN",0,"cnfg file [%s]\n",cnfg_filename); 


  /* read & broadcast parameters */
  parse_cnfg_filename(cnfg_filename,&fpars);

/*
 * x Agostino: Serve veramente??
 * Claudio

#define remove_parameter(NAME,PAR) \
  { \
    for(i=0;(PAR).read[i].name!=NULL;i++) { \
      if(strcmp((PAR).read[i].name,#NAME)==0) { \
	(PAR).read[i].descr=NULL; \
	break; \
      } \
    } \
  }

  if(fpars.type==DYNAMICAL_CNFG || fpars.type==QUENCHED_CNFG) {
    remove_parameter(GLB_T,glb_var);
    remove_parameter(GLB_X,glb_var);
    remove_parameter(GLB_Y,glb_var);
    remove_parameter(GLB_Z,glb_var);
  }
  if(fpars.type==DYNAMICAL_CNFG) remove_parameter(quark quenched masses,mes_var);
#undef remove_parameter
*/

  read_input(mes_var.read,input_filename);
  GLB_T=fpars.t; GLB_X=fpars.x; GLB_Y=fpars.y; GLB_Z=fpars.z;
  error(fpars.type==UNKNOWN_CNFG,1,"mk_mesons.c","Bad name for a configuration file");
  error(fpars.nc!=NG,1,"mk_mesons.c","Bad NG");


  lprintf("MAIN",0,"RLXD [%d,%d]\n",rlx_var.rlxd_level,rlx_var.rlxd_seed);
  rlxd_init(rlx_var.rlxd_level,rlx_var.rlxd_seed+PID);
  srand(rlx_var.rlxd_seed+PID);

  lprintf("MAIN",0,"Gauge group: SU(%d)\n",NG);
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

/*
  tricorr=(double*)malloc(GLB_T*sizeof(double));
  pta_qprop=(spinor_field**)malloc(sizeof(spinor_field*)*nm);
  pta_qprop[0]=alloc_spinor_field_f(4*NF*nm,&glattice);
  for(k=0;k<nm;++k)
    pta_qprop[k]=pta_qprop[0]+4*NF*k;
*/

  list=NULL;
  if(strcmp(list_filename,"")!=0) {
    error((list=fopen(list_filename,"r"))==NULL,1,"main [mk_mesons.c]" ,
	"Failed to open list file\n");
  }


  i=0;
  while(1) {
	  int nn;

	  if(list!=NULL)
		  if(fscanf(list,"%s",cnfg_filename)==0 || feof(list)) break;

	  i++;

	  lprintf("MAIN",0,"Configuration from %s\n", cnfg_filename);
	  /* NESSUN CHECK SULLA CONSISTENZA CON I PARAMETRI DEFINITI !!! */
	  read_gauge_field(cnfg_filename);
	  represent_gauge_field();

	  lprintf("TEST",0,"<p> %1.6f\n",avr_plaquette());

	  full_plaquette();

	  for (nn=0;nn<mes_var.nhits;++nn) {
		  lprintf("MAIN",0,"Configuration #%d hit %d\n",i,nn);
		  inline_mk_mesons(m,nm,mes_var.precision);
	  }

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

