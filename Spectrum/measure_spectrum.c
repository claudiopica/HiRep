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
#include "cinfo.c"




#if defined(ROTATED_SF) && defined(BASIC_SF)
#error This code does not work with the Schroedinger functional !!!
#endif

#ifdef FERMION_THETA
#error This code does not work with the fermion twisting !!!
#endif

//typedef enum {semwall_src,point_src} source_type_t;
/* Mesons parameters */
typedef struct _input_mesons {
	char mstring[1024];
	int use_input_mass;
	double precision;
	int meas_mixed;
	int nhits_2pt;
	int nhits_disc;
	int def_semwall;
	int def_point;
	int def_baryon;
	int def_glueball;
	int ext_semwall;
	int ext_point;
	int fixed_semwall;
	int fixed_point;
	int fixed_gfwall;
	int discon_semwall;
	int discon_gfwall;
	int discon_volume;
	int def_gfwall;
	int dt;
	int n_mom;
	int dilution;
	int background_field;
	int nEz;
	double Q;
	double csw;
	double rho_s;
	double rho_t;

	//Currently only implemented for ff
	int nhits_hopping;  //Multiplies the number of hits in the fast part of the hopping parameter expansion
	int degree_hopping;  // The degree of the hopping parameter expasion

	/* for the reading function */
	input_record_t read[31];
} input_mesons;

#define init_input_mesons(varname) \
{ \
  .read={\
    {"quark quenched masses", "mes:masses = %s", STRING_T, (varname).mstring}, \
    {"use input mass", "mes:use_input_mass = %d",INT_T, &(varname).use_input_mass},\
    {"inverter precision", "mes:precision = %lf", DOUBLE_T, &(varname).precision},\
    {"measure mixed correlators", "mes:meas_mixed = %d",INT_T,&(varname).meas_mixed},\
    {"number of noisy sources per cnfg for 2pt fn", "mes:nhits_2pt = %d", INT_T, &(varname).nhits_2pt}, \
    {"number of noisy sources per cnfg for disconnected", "mes:nhits_disc = %d", INT_T, &(varname).nhits_disc}, \
    {"enable default semwall", "mes:def_semwall = %d",INT_T, &(varname).def_semwall},	\
    {"enable default point", "mes:def_point = %d",INT_T, &(varname).def_point},		\
    {"enable default gfwall", "mes:def_gfwall = %d",INT_T, &(varname).def_gfwall},	\
    {"enable extended semwall", "mes:ext_semwall = %d",INT_T, &(varname).ext_semwall},	\
    {"enable extended point", "mes:ext_point = %d",INT_T, &(varname).ext_point},		\
    {"enable Dirichlet semwall", "mes:dirichlet_semwall = %d",INT_T, &(varname).fixed_semwall},	\
    {"enable Dirichlet point", "mes:dirichlet_point = %d",INT_T, &(varname).fixed_point},	\
    {"enable Dirichlet gfwall", "mes:dirichlet_gfwall = %d",INT_T, &(varname).fixed_gfwall},	\
    {"enable discon semwall", "mes:discon_semwall = %d",INT_T, &(varname).discon_semwall},	\
    {"enable discon gfwall", "mes:discon_gfwall = %d",INT_T, &(varname).discon_gfwall},	\
    {"enable discon volume", "mes:discon_volume = %d",INT_T, &(varname).discon_volume},	\
    {"volume source dilution", "mes:dilution = %d",INT_T, &(varname).dilution},	\
    {"Distance of t_initial from Dirichlet boundary", "mes:dirichlet_dt = %d", INT_T, &(varname).dt},\
    {"maximum component of momentum", "mes:momentum = %d", INT_T, &(varname).n_mom}, \
    {"enable baryon", "mes:def_baryon = %d",INT_T, &(varname).def_baryon},		\
    {"enable glueball", "mes:def_glueball = %d",INT_T, &(varname).def_glueball},		\
    {"enable background electric field", "mes:background_field = %d",INT_T, &(varname).background_field},	\
    {"electric charge", "mes:Q = %lf",DOUBLE_T, &(varname).Q},	\
    {"electric field nEz", "mes:nEz = %d",INT_T, &(varname).nEz},	\
    {"csw coefficient", "mes:csw = %lg",DOUBLE_T, &(varname).csw},	\
    {"smearing space", "mes:rho_s = %lg",DOUBLE_T, &(varname).rho_s},	\
    {"smearing time", "mes:rho_t = %lg",DOUBLE_T, &(varname).rho_t},	\
    {"hopping expansion degree", "mes:degree_hopping = %d",INT_T, &(varname).degree_hopping}, \
    {"hopping expansion hits", "mes:nhits_hopping = %d",INT_T, &(varname).nhits_hopping}, \
    {NULL, NULL, INT_T, NULL}				\
   }							\
}

/* Turn four fermion on/off */
typedef struct _input_ff {
  char make[256]; 
  /* for the reading function */
  input_record_t read[2];
  
} input_ff;

#define init_input_ff(varname) \
{ \
  .read={\
    {"Include four fermion interactions", "ff:on = %s", STRING_T, (varname).make},\
    {NULL, NULL, INT_T, NULL}\
  }\
}



char cnfg_filename[256]="";
char list_filename[256]="";
char input_filename[256] = "input_file";
char output_filename[256] = "mesons.out";
enum { UNKNOWN_CNFG, DYNAMICAL_CNFG, QUENCHED_CNFG };

input_mesons mes_var = init_input_mesons(mes_var);
input_ff ff_var = init_input_ff(ff_var);



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
  int i,k,tau;
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
  read_input(ff_var.read,input_filename);
  read_input(mes_var.read,input_filename);

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

  if(mes_var.use_input_mass) {
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
#ifdef WITH_CLOVER
  clover_init(mes_var.csw);
#endif
#ifdef WITH_SMEARING
	init_smearing(mes_var.rho_s, mes_var.rho_t);
#endif

  lprintf("MAIN",0,"Inverter precision = %e\n",mes_var.precision);
  lprintf("MAIN",0,"Mass[%d] = %f",0,m[0]);
  for(k=1;k<nm;k++)
    lprintf("MAIN",0,", Mass[%d] = %f",k,m[k]);
  lprintf("MAIN",0,"\n",k,m[k]);
  lprintf("MAIN",0,"Number of noisy sources per cnfg = %d. Does not affect point sources\n",mes_var.nhits_2pt);
  if (mes_var.def_semwall){
    lprintf("MAIN",0,"Spin Explicit Method (SEM) wall sources\n");    
  }
  if (mes_var.def_point){
    lprintf("MAIN",0,"Point sources\n");    
  }
 	if (mes_var.def_baryon){
    lprintf("MAIN",0,"Baryon masses\n");    
  }
 	if (mes_var.def_glueball){
    lprintf("MAIN",0,"Glueball masses\n");    
  }
  if (mes_var.def_gfwall){
    lprintf("MAIN",0,"Gauge Fixed Wall Source\n");    
  }
  if (mes_var.ext_semwall){
    lprintf("MAIN",0,"Spin Explicit Method (SEM) wall sources on extended lattice\n");    
  }
  if (mes_var.ext_point){
    lprintf("MAIN",0,"Point sources on extended lattice\n");    
  }
  if (mes_var.fixed_semwall){
    lprintf("MAIN",0,"Spin Explicit Method (SEM) wall sources with Dirichlet boundary conditions\n"); 
    lprintf("MAIN",0,"Distance between tau and the boundary: %d\n",mes_var.dt); 
  }
  if (mes_var.fixed_point){
    lprintf("MAIN",0,"Point sources with Dirichlet boundary conditions\n");    
    lprintf("MAIN",0,"Distance between tau and the boundary: %d\n",mes_var.dt); 
  }
  if (mes_var.fixed_gfwall){
    lprintf("MAIN",0,"Gauge Fixed Wall Source with Dirichlet boundary conditions\n"); 
    lprintf("MAIN",0,"Distance between tau and the boundary: %d\n",mes_var.dt); 
  }
  if (mes_var.n_mom>1){
    lprintf("MAIN",0,"Number of maximum monentum component %d\n",mes_var.n_mom-1);
    if (mes_var.def_semwall || mes_var.ext_semwall || mes_var.fixed_semwall){
      lprintf("MAIN",0,"WARGING: wall sources measure only with zero momenta\n");
    }
  }
  if (mes_var.n_mom==0){ mes_var.n_mom++;}
 	if (mes_var.background_field){
    lprintf("MAIN",0,"Electric background field in the z direction with charge Q=%1.6f , with E =  %d x 2 pi /(Q*T*L) \n",mes_var.Q,mes_var.nEz);    
  }
 
  list=NULL;
  if(strcmp(list_filename,"")!=0) {
    error((list=fopen(list_filename,"r"))==NULL,1,"main [mk_mesons.c]" ,
	"Failed to open list file\n");
  }

  if (mes_var.meas_mixed){
    init_meson_correlators(1);
    lprintf("MAIN",0,"Measuring all 256 correlators %d\n");
  }
  else{
    lprintf("MAIN",0,"Measuring Gamma Gamma correlators and PCAC-mass\n");
    init_meson_correlators(0);
  }

  if(strcmp(ff_var.make,"true")==0){
    four_fermion_active = 1;
    ff_sigma = alloc_sfield(1,&glattice);
    ff_pi = alloc_sfield(1,&glattice); 
  }

  init_discon_correlators();
  init_cvc_correlators();
  if(four_fermion_active==1) init_triplet_discon_correlators();
  i=0;

  while(++i) {
    struct timeval start, end, etime;

    if(list!=NULL)
      if(fscanf(list,"%s",cnfg_filename)==0 || feof(list)) break;


    lprintf("MAIN",0,"Configuration from %s\n", cnfg_filename);
    read_gauge_field(cnfg_filename);
    represent_gauge_field();

    lprintf("TEST",0,"<p> %1.6f\n",avr_plaquette());

		// if non zero background field : apply abelian field and boundary correction. Then measure all plaquettes.
		if (mes_var.background_field){
//		apply_background_field_zdir(u_gauge,mes_var.Q,mes_var.nEz);
	      measure_diquark_semwall_background(nm,m,mes_var.nhits_2pt,i,mes_var.precision,mes_var.Q,mes_var.nEz);
		}
    full_plaquette();
    gettimeofday(&start,0);

    if(four_fermion_active==1) ff_observables(); 

    tau=0;
    if (four_fermion_active==0) {

     if (mes_var.def_semwall){
       measure_spectrum_semwall(nm,m,mes_var.nhits_2pt,i,mes_var.precision);
     }
     if (mes_var.def_point){
       measure_spectrum_pt(tau,nm,m,mes_var.n_mom,mes_var.nhits_2pt,i,mes_var.precision);
     }
     if (mes_var.def_baryon){
       measure_baryons(m,i,mes_var.precision);
     }
     if (mes_var.def_glueball){
       //measure_glueballs(); //This does not seem to exist 
     }
     if (mes_var.def_gfwall){
       measure_spectrum_gfwall(nm,m,i,mes_var.precision);
     }
     if (mes_var.ext_semwall){
      measure_spectrum_semwall_ext(nm,m,mes_var.nhits_2pt,i,mes_var.precision);
     }
     if (mes_var.ext_point){
       measure_spectrum_pt_ext(tau,nm,m,mes_var.n_mom,mes_var.nhits_2pt,i,mes_var.precision);
     }
     if (mes_var.fixed_semwall){
       measure_spectrum_semwall_fixedbc(mes_var.dt,nm,m,mes_var.nhits_2pt,i,mes_var.precision);
     }
     if (mes_var.fixed_point){
       measure_spectrum_pt_fixedbc(tau,mes_var.dt,nm,m,mes_var.n_mom,mes_var.nhits_2pt,i,mes_var.precision);
     }
     if (mes_var.fixed_gfwall){
       measure_spectrum_gfwall_fixedbc(mes_var.dt,nm,m,i,mes_var.precision);
     }
     if (mes_var.discon_semwall){
        measure_spectrum_discon_semwall(nm,m,mes_var.nhits_disc,i,mes_var.precision); 
     }
     if (mes_var.discon_gfwall){
       measure_spectrum_discon_gfwall(nm,m,i,mes_var.precision);
     }
     if (mes_var.discon_volume){
       measure_spectrum_discon_volume(nm,m,i,mes_var.precision,mes_var.dilution);
     }

    } else {
     //With four fermion interactions
     if (mes_var.def_semwall){
       measure_spectrum_ff_semwall(nm,m,mes_var.nhits_2pt,i,mes_var.precision);
     }
     if (mes_var.def_point){
       measure_spectrum_ff_pt(tau,nm,m,mes_var.n_mom,mes_var.nhits_2pt,i,mes_var.precision);
     }
     if (mes_var.ext_semwall){
       measure_spectrum_semwall_ff_ext(nm,m,mes_var.nhits_2pt,i,mes_var.precision);
     }
     if (mes_var.discon_semwall){
       measure_spectrum_discon_ff_semwall(nm,m,mes_var.nhits_disc,mes_var.degree_hopping,mes_var.nhits_hopping,i,mes_var.precision);
     }

    }

    gettimeofday(&end,0);
    timeval_subtract(&etime,&end,&start);
    lprintf("MAIN",0,"Configuration #%d: analysed in [%ld sec %ld usec]\n",i,etime.tv_sec,etime.tv_usec);
    if(list==NULL) break;
  }

  if(list!=NULL) fclose(list);
  free_meson_observables();
  free_BCs();

  if(four_fermion_active==1) free_triplet_discon_observables();

  free_gfield(u_gauge);
#ifdef ALLOCATE_REPR_GAUGE_FIELD
  free_gfield_f(u_gauge_f);
#endif

  finalize_process();

  return 0;
}

