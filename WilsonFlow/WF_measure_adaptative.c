/*******************************************************************************
*
* Computation of the observable E(t) evolved with the Wilson flow
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
#include "random.h"
#include "wilsonflow.h"

#include "cinfo.c"


typedef struct _input_WF {
  double tmax;
  int nmeas;
  double eps;
	double delta;
  int def_glueball;

  /* for the reading function */
  input_record_t read[6];

} input_WF;

#define init_input_WF(varname) \
{ \
  .read={\
    {"WF max integration time", "WF:tmax = %lf", DOUBLE_T, &((varname).tmax)},\
    {"WF number of measures", "WF:nmeas = %d", DOUBLE_T, &((varname).nmeas)},\
    {"WF initial epsilon", "WF:eps = %lf", DOUBLE_T, &((varname).eps)},\
    {"WF delta", "WF:delta = %lf", DOUBLE_T, &((varname).delta)},\
		{"enable glueball", "WF:def_glueball = %d",INT_T, &(varname).def_glueball},                \
    {NULL, NULL,0,NULL}\
  }\
}

input_WF WF_var = init_input_WF(WF_var);


#ifdef ROTATED_SF

typedef struct _input_SF {
  double ct;
  double beta;

  /* for the reading function */
  input_record_t read[3];

} input_SF;

#define init_input_SF(varname) \
{ \
  .read={\
    {"SF ct", "SF:ct = %lf", DOUBLE_T, &((varname).ct)},\
    {"SF beta", "SF:beta = %lf", DOUBLE_T, &((varname).beta)},\
    {NULL, NULL, 0, NULL}\
  }\
}

input_SF SF_var = init_input_SF(SF_var);

#endif /* ROTATED_SF */


char cnfg_filename[256]="";
char list_filename[256]="";
char input_filename[256] = "input_file";
char output_filename[256] = "wilsonflow.out";
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

  error((ac==0 && al==0) || (ac!=0 && al!=0),1,"parse_cmdline [WF_measure.c]",
      "Syntax: mk_wilsonloops { -c <config file> | -l <list file> } [-i <input file>] [-o <output file>] [-m]");

  if(ac != 0) {
    strcpy(cnfg_filename,argv[ac]);
    strcpy(list_filename,"");
  } else if(al != 0) {
    strcpy(list_filename,argv[al]);
    error((list=fopen(list_filename,"r"))==NULL,1,"parse_cmdline [WF_measure.c]" ,
	"Failed to open list file\n");
    error(fscanf(list,"%s",cnfg_filename)==0,1,"parse_cmdline [WF_measure.c]" ,
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
  read_input(rlx_var.read,input_filename);
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

  read_input(WF_var.read,input_filename);
#ifdef ROTATED_SF
  read_input(SF_var.read,input_filename);
#endif
  GLB_T=fpars.t; GLB_X=fpars.x; GLB_Y=fpars.y; GLB_Z=fpars.z;

  error(fpars.type==UNKNOWN_CNFG,1,"WF_measure.c","Bad name for a configuration file");
  error(fpars.nc!=NG,1,"WF_measure.c","Bad NG");

  lprintf("MAIN",0,"RLXD [%d,%d]\n",rlx_var.rlxd_level,rlx_var.rlxd_seed);
  rlxd_init(rlx_var.rlxd_level,rlx_var.rlxd_seed+PID);
  srand(rlx_var.rlxd_seed+PID);

#ifdef GAUGE_SUN
  lprintf("MAIN",0,"Gauge group: SU(%d)\n",NG);
#elif GAUGE_SON
  lprintf("MAIN",0,"Gauge group: SO(%d)\n",NG);
#else
  lprintf("MAIN",0,"Default gauge group: SU(%d)\n",NG);
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
  
  BCs_pars_t BCs_pars = {
    .fermion_twisting_theta = {0.,0.,0.,0.},
    .gauge_boundary_improvement_cs = 1.,
#ifdef ROTATED_SF
    .gauge_boundary_improvement_ct = SF_var.ct,
#else
    .gauge_boundary_improvement_ct = 1.,
#endif
    .chiSF_boundary_improvement_ds = 1.,
    .SF_BCs = 1
  };

  init_BCs(&BCs_pars);

  /* alloc global gauge fields */
  u_gauge=alloc_gfield(&glattice);

		double dt = (double)WF_var.tmax/(double)WF_var.nmeas;
  
  lprintf("MAIN",0,"WF tmax: %e\n",WF_var.tmax);
  lprintf("MAIN",0,"WF number of measures: %d\n",WF_var.nmeas);
  lprintf("MAIN",0,"WF initial epsilon: %e\n",WF_var.eps);
  lprintf("MAIN",0,"WF delta: %e\n",WF_var.delta);
  lprintf("MAIN",0,"WF measurement interval dt : %e\n",dt);


#ifdef ROTATED_SF
  lprintf("MAIN",0,"SF beta=%e\n",SF_var.beta);      
  lprintf("MAIN",0,"SF ct=%e\n",SF_var.ct);
#endif

  
  list=NULL;
  if(strcmp(list_filename,"")!=0) {
    error((list=fopen(list_filename,"r"))==NULL,1,"main [WF_measure.c]" ,
	"Failed to open list file\n");
  }

  WF_initialize();

#ifndef ROTATED_SF
  double E, Esym, TC;
#else
  int j;
  double E[2*GLB_T];
  double Esym[2*GLB_T];
  double Eavg[2];
  double Esymavg[2];
#endif

  struct timeval start, end, etime;

	
  i=0;
  while(1) {

    if(list!=NULL)
      if(fscanf(list,"%s",cnfg_filename)==0 || feof(list)) break;

    i++;

    lprintf("MAIN",0,"Configuration from %s\n", cnfg_filename);
    /* NESSUN CHECK SULLA CONSISTENZA CON I PARAMETRI DEFINITI !!! */

		gettimeofday(&start,0);
    read_gauge_field(cnfg_filename);

    apply_BCs_on_fundamental_gauge_field();

    full_plaquette();

    int k, n;
    double epsilon=WF_var.eps;
    double t=0.;

//  measurement at 0 wilson flow time
#ifndef ROTATED_SF

      E=WF_E(u_gauge);
      Esym=WF_Esym(u_gauge);
      TC=WF_topo(u_gauge);
			lprintf("WILSONFLOW",0,"WF (ncnfg,t,E,t2*E,Esym,t2*Esym,TC) = %d %e %e %e %e %e %e\n",i,t,E,t*t*E,Esym,t*t*Esym,TC);
//			if (WF_var.def_glueball){
//							measure_glueballs();
//			}
#else

			WF_E_T(E,u_gauge);
			WF_Esym_T(Esym,u_gauge);
			Eavg[0]=Eavg[1]=Esymavg[0]=Esymavg[1]=0.0;
			for(j=1;j<GLB_T-1;j++){
							lprintf("WILSONFLOW",0,"WF (ncnfg,T,t,Etime,Espace,Esymtime,Esymspace) = %d %d %e %e %e %e %e\n",i,j,t,E[2*j],E[2*j+1],Esym[2*j],Esym[2*j+1]);
							Eavg[0] += E[2*j];
							Eavg[1] += E[2*j+1];
							Esymavg[0] += Esym[2*j];
							Esymavg[1] += Esym[2*j+1];
			}

			Eavg[0] /= GLB_T-2;
			Eavg[1] /= GLB_T-3;
			Esymavg[0] /= GLB_T-2;
			Esymavg[1] /= GLB_T-3;

			lprintf("WILSONFLOW",0,"WF avg (ncnfg,t,Etime,Espace,Esymtime,Esymspace,Pltime,Plspace) = %d %e %e %e %e %e %e %e\n",i,t,Eavg[0],Eavg[1],Esymavg[0],Esymavg[1],(NG-Eavg[0]),(NG-Eavg[1]));
			lprintf("WILSONFLOW",0,"SF dS/deta= %e\n", SF_action(SF_var.beta));

#endif
		k=1;	
		double epsilon_new=0;
		while (t < WF_var.tmax)
		{	
          if (t+epsilon > (double)k*dt) 	epsilon = (double)k*dt - t; 

					epsilon_new=WilsonFlow3_adaptative(u_gauge,epsilon,WF_var.delta);
	      
          if ( fabs(epsilon_new+1.) > 1e-7) t=t+epsilon;
				
					if ( fabs(t - (double)k*dt ) < 1e-7 ) {
						k=k+1;
#ifndef ROTATED_SF

		     	 E=WF_E(u_gauge);
    		 	 Esym=WF_Esym(u_gauge);
      		 TC=WF_topo(u_gauge);
					 lprintf("WILSONFLOW",0,"WF (ncnfg,t,E,t2*E,Esym,t2*Esym,TC) = %d %e %e %e %e %e %e\n",i,t,E,t*t*E,Esym,t*t*Esym,TC);
//					 if (WF_var.def_glueball){
//							measure_glueballs();
//					 }
#else

					 WF_E_T(E,u_gauge);
					 WF_Esym_T(Esym,u_gauge);
					 Eavg[0]=Eavg[1]=Esymavg[0]=Esymavg[1]=0.0;
					 for(j=1;j<GLB_T-1;j++){
							lprintf("WILSONFLOW",0,"WF (ncnfg,T,t,Etime,Espace,Esymtime,Esymspace) = %d %d %e %e %e %e %e\n",i,j,t,E[2*j],E[2*j+1],Esym[2*j],Esym[2*j+1]);
							Eavg[0] += E[2*j];
							Eavg[1] += E[2*j+1];
							Esymavg[0] += Esym[2*j];
							Esymavg[1] += Esym[2*j+1];
					 }

					Eavg[0] /= GLB_T-2;
					Eavg[1] /= GLB_T-3;
					Esymavg[0] /= GLB_T-2;
					Esymavg[1] /= GLB_T-3;

					lprintf("WILSONFLOW",0,"WF avg (ncnfg,t,Etime,Espace,Esymtime,Esymspace,Pltime,Plspace) = %d %e %e %e %e %e %e %e\n",i,t,Eavg[0],Eavg[1],Esymavg[0],Esymavg[1],(NG-Eavg[0]),(NG-Eavg[1]));
					lprintf("WILSONFLOW",0,"SF dS/deta= %e\n", SF_action(SF_var.beta));

#endif

					}
					if (fabs(epsilon_new + 1.) > 1e-7) epsilon=epsilon_new;	
					if (fabs(epsilon_new +1.) < 1e-7 ) epsilon=epsilon/2;	
    
 }
	gettimeofday(&end,0);
	timeval_subtract(&etime,&end,&start);
	lprintf("TIMING",0,"Wilson Flow evolution and measurements for configuration  [%s] done [%ld sec %ld usec]\n",cnfg_filename,etime.tv_sec,etime.tv_usec);

	

		if(list==NULL) break;
	} // end loop configurations

	if(list!=NULL) fclose(list);

	WF_free();

	free_BCs();

	free_gfield(u_gauge);

	finalize_process();

	return 0;
}

