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


char cnfg_filename[256];
char list_filename[256];
char input_filename[256] = "input_file";
char output_filename[256] = "mesons.out";
enum { UNKNOWN_CNFG, DYNAMICAL_CNFG, QUENCHED_CNFG };

input_glb glb_ip = init_input_glb(glb_ip);
input_mesons mes_ip = init_input_mesons(mes_ip);


typedef struct {
  char string[256];
  int t, x, y, z;
  int nc, nf;
  float b, k;
  int n;
  int type;
} filename_t;


int parse_cnfg_filename(char* filename, filename_t* fn) {
  int hm;
  char *tmp = NULL;
  char *basename;
  
  basename = filename;
  tmp = strchr(filename, '/');
  while(tmp != NULL) {
    basename = tmp+1;
    tmp = strchr(tmp+1, '/');
  }            

  hm=sscanf(basename,"%dx%dx%dx%dNc%dNf%db%fk%fn%d",
            &(fn->t),&(fn->x),&(fn->y),&(fn->z),&(fn->nc),&(fn->nf),&(fn->b),&(fn->k),&(fn->n));
  if(hm==9) {
    fn->type=DYNAMICAL_CNFG;
    return DYNAMICAL_CNFG;
  }

  hm=sscanf(basename,"%dx%dx%dx%dNc%db%fn%d",
            &(fn->t),&(fn->x),&(fn->y),&(fn->z),&(fn->nc),&(fn->b),&(fn->n));
  if(hm==7) {
    fn->type=QUENCHED_CNFG;
    return QUENCHED_CNFG;
  }
	
	fn->type=UNKNOWN_CNFG;
	return UNKNOWN_CNFG;
}


void read_cmdline(int argc, char* argv[]) {
  int i, ai=0, ao=0, ac=0, al=0;
  FILE *list=NULL;
  
  for (i=1;i<argc;i++) {
    if (strcmp(argv[i],"-i")==0) ai=i+1;
    else if (strcmp(argv[i],"-o")==0) ao=i+1;
    else if (strcmp(argv[i],"-c")==0) ac=i+1;
    else if (strcmp(argv[i],"-l")==0) al=i+1;
  }

  if (ao!=0) strcpy(output_filename,argv[ao]);
  if (ai!=0) strcpy(input_filename,argv[ai]);

  error((ac==0 && al==0) || (ac!=0 && al!=0),1,"parse_cmdline [mk_mesons.c]",
        "Syntax: mk_mesons { -c <config file> | -l <list file> } [-i <input file>] [-o <output file>]");

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

  lprintf("MAIN",0,"input file %s\n",input_filename); 
  lprintf("MAIN",0,"output file %s\n",output_filename); 
  lprintf("MAIN",0,"cnfg file %s\n",cnfg_filename); 
  lprintf("MAIN",0,"list file %s\n",list_filename); 

}


int main(int argc,char *argv[]) {
	int i,k,n;
  char tmp[256], *cptr;
  FILE* list;
	spinor_field **pta_qprop=0;
 	double* tricorr;
  filename_t fpars;
  int nm;
  double m[256];
  spinor_field *test;

  /* setup process id and communications */
  setup_process(&argc,&argv);

  /* logger setup */
  /* disable logger for MPI processes != 0 */
  if (PID!=0) { logger_disable(); }
  logger_setlevel(0,40);
  if (PID==0) sprintf(tmp,">%s",output_filename); logger_stdout(tmp);
  sprintf(tmp,"err_%d",PID); freopen(tmp,"w",stderr);

  lprintf("MAIN",0,"PId =  %d [world_size: %d]\n\n",PID,WORLD_SIZE); 


  /* read & broadcast parameters */
  read_cmdline(argc, argv);
  parse_cnfg_filename(cnfg_filename,&fpars);
  GLB_T=fpars.t; GLB_X=fpars.x; GLB_Y=fpars.y; GLB_Z=fpars.z;
  error(fpars.type==UNKNOWN_CNFG,1,"mk_mesons.c","Bad name for a configuration file");
  error(fpars.nc!=NG,1,"mk_mesons.c","Bad NG");


#define remove_parameter(NAME,PAR) \
{ \
  for(i=0;(PAR).read[i].name!=NULL;i++) { \
    if(strcmp((PAR).read[i].name,#NAME)==0) { \
      (PAR).read[i].descr=NULL; \
      break; \
    } \
  } \
}
  
  remove_parameter(GLB_T,glb_ip);
  remove_parameter(GLB_X,glb_ip);
  remove_parameter(GLB_Y,glb_ip);
  remove_parameter(GLB_Z,glb_ip);
  if(fpars.type==DYNAMICAL_CNFG) remove_parameter(quark quenched masses,mes_ip);

#undef remove_parameter

  read_input(glb_ip.read,input_filename);
  read_input(mes_ip.read,input_filename);

  nm=0;
  if(fpars.type==DYNAMICAL_CNFG) {
    nm=1;
	  m[0] = 0.5/fpars.k - 4.0;
	} else if(fpars.type==QUENCHED_CNFG) {
    strcpy(tmp,mes_ip.mstring);
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

  lprintf("MAIN",0,"Gauge group: SU(%d)\n",NG);
  lprintf("MAIN",0,"Fermion representation: " REPR_NAME " [dim=%d]\n",NF);
  lprintf("MAIN",0,"global size is %dx%dx%dx%d\n",GLB_T,GLB_X,GLB_Y,GLB_Z);
  lprintf("MAIN",0,"proc grid is %dx%dx%dx%d\n",NP_T,NP_X,NP_Y,NP_Z);

  /* setup lattice geometry */
  geometry_mpi_eo();
  /* test_geometry_mpi_eo(); */

  lprintf("MAIN",0,"local size is %dx%dx%dx%d\n",T,X,Y,Z);
  lprintf("MAIN",0,"extended local size is %dx%dx%dx%d\n",T_EXT,X_EXT,Y_EXT,Z_EXT);

  lprintf("MAIN",0,"RLXD [%d,%d]\n",glb_ip.rlxd_level,glb_ip.rlxd_seed);
  rlxd_init(glb_ip.rlxd_level,glb_ip.rlxd_seed+PID);
 
  /* alloc global gauge fields */
  u_gauge=alloc_gfield(&glattice);
#ifndef REPR_FUNDAMENTAL
  u_gauge_f=alloc_gfield_f(&glattice);
#endif

  for(k=0;k<nm;k++)
    lprintf("MAIN",0,"Mass[%d] = %f\n",k,m[k]);

  tricorr=(double*)malloc(GLB_T*sizeof(double));
	pta_qprop=(spinor_field**)malloc(sizeof(spinor_field*)*nm);
	pta_qprop[0]=alloc_spinor_field_f(4*NF*nm,&glattice);
	for(k=0;k<nm;++k)
		pta_qprop[k]=pta_qprop[0]+4*NF*k;
	
  list=NULL;
  if(strcmp(list_filename,"")!=0) {
    error((list=fopen(list_filename,"r"))==NULL,1,"main [mk_mesons.c]" ,
        "Failed to open list file\n");
  }

  test=alloc_spinor_field_f(1,&glattice);

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

	  pta_qprop_QMR_eo(pta_qprop, nm, m, 1e-9);
	
	  for (k=0;k<nm;++k){

	    lprintf("MAIN",0,"conf #%d mass=%2.6f \n",i,m[k]);

#define CORR(name) \
	name##_correlator(tricorr, pta_qprop[k]);\
	lprintf("MAIN",0,"conf #%d mass=%2.6f TRIPLET " #name "= ",i,m[k]);\
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
	  
	  if(list==NULL) break;
	}

  if(list!=NULL) fclose(list);

	free_spinor_field(pta_qprop[0]);
	free(pta_qprop);
	free(tricorr);

	free_gfield(u_gauge);
#ifndef REPR_FUNDAMENTAL
	free_gfield_f(u_gauge_f);
#endif

	return 0;
}
