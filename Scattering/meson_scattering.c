/*******************************************************************************
 *
 * Computation of the Pseudo scalar scattering length
 * VD 2014
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
#include "gamma_spinor.h"
#include "spin_matrix.h"
#include "scattering.h"

#define PI 3.141592653589793238462643383279502884197


#include "cinfo.c"

#if defined(ROTATED_SF) && defined(BASIC_SF)
#error This code does not work with the Schroedinger functional !!!
#endif

#ifdef FERMION_THETA
#error This code does not work with the fermion twisting !!!
#endif
/* Mesons parameters */
typedef struct _input_scatt {
	char mstring[256];
	double precision;
	int nhits;
	int tsrc;

	/* for the reading function */
	input_record_t read[5];

} input_scatt;

#define init_input_scatt(varname) \
{ \
	.read={\
		{"quark quenched masses", "mes:masses = %s", STRING_T, (varname).mstring},\
		{"inverter precision", "mes:precision = %lf", DOUBLE_T, &(varname).precision},\
		{"number of inversions per cnfg", "mes:nhits = %d", INT_T, &(varname).nhits},\
		{"Source time:", "mes:tsrc = %d", INT_T, &(varname).tsrc},\
		{NULL, NULL, INT_T, NULL}\
	}\
}


char cnfg_filename[256]="";
char list_filename[256]="";
char prop_filename[256]="";
char source_filename[256]="";
char input_filename[256] = "input_file";
char output_filename[256] = "meson_scattering.out";
int Nsource;
double M;

enum { UNKNOWN_CNFG, DYNAMICAL_CNFG, QUENCHED_CNFG };

input_scatt mes_var = init_input_scatt(mes_var);

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
	int i, ai=0, ao=0, ac=0, al=0, am=0,ap=0,as=0;
	FILE *list=NULL;

	for (i=1;i<argc;i++) {
		if (strcmp(argv[i],"-i")==0) ai=i+1;
		else if (strcmp(argv[i],"-o")==0) ao=i+1;
		else if (strcmp(argv[i],"-c")==0) ac=i+1;
		else if (strcmp(argv[i],"-p")==0) ap=i+1;
		else if (strcmp(argv[i],"-l")==0) al=i+1;
		else if (strcmp(argv[i],"-s")==0) as=i+1;
		else if (strcmp(argv[i],"-m")==0) am=i;
	}


	if (am != 0) {
		print_compiling_info();
		exit(0);
	}


	if (ao!=0) strcpy(output_filename,argv[ao]);
	if (ai!=0) strcpy(input_filename,argv[ai]);

	error ((ap != 0 && as ==0),1,"parse_cmdline [discc]",
			"Syntax: disc { -c <config file> | -l <list file> } [-i <input file>] [-o <output file>] [-m] -p <propagator name> -s <source_name> ");

	error((ac==0 && al==0) || (ac!=0 && al!=0),1,"parse_cmdline [meson_scattering.c]",
			"Syntax: disc { -c <config file> | -l <list file> } [-i <input file>] [-o <output file>] [-m] -p <propagator name> -s <source_name> ");

	if(ap != 0) {
		strcpy(prop_filename,argv[ap]);
	} 
	if(as != 0) {
		strcpy(source_filename,argv[as]);
	} 

	if(ac != 0) {
		strcpy(cnfg_filename,argv[ac]);
		strcpy(list_filename,"");
	} else if(al != 0) {
		strcpy(list_filename,argv[al]);
		error((list=fopen(list_filename,"r"))==NULL,1,"parse_cmdline [meson_scattering.c]" ,
				"Failed to open list file\n");
		error(fscanf(list,"%s",cnfg_filename)==0,1,"parse_cmdline [meson_scattering.c]" ,
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



	struct timeval start, end, etime;
	gettimeofday(&start,0);
	/* setup process id and communications */
	read_cmdline(argc, argv);
	setup_process(&argc,&argv);

	read_input(glb_var.read,input_filename);
	setup_replicas();

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

	read_input(mes_var.read,input_filename);
	GLB_T=fpars.t; GLB_X=fpars.x; GLB_Y=fpars.y; GLB_Z=fpars.z;

	/* setup lattice geometry */
	if (geometry_init() == 1) { finalize_process(); return 0; }
	 // geometry_mpi_eo();
	/* test_geometry_mpi_eo(); */ 
	// test_geometry_mpi_eo();  
	/* setup random numbers */
	read_input(rlx_var.read,input_filename);
	//slower(rlx_var.rlxd_start); //convert start variable to lowercase
	if(strcmp(rlx_var.rlxd_start,"continue")==0 && rlx_var.rlxd_state[0]!='\0') {
		/*load saved state*/
		lprintf("MAIN",0,"Loading rlxd state from file [%s]\n",rlx_var.rlxd_state);
		read_ranlxd_state(rlx_var.rlxd_state);
	} else {
		lprintf("MAIN",0,"RLXD [%d,%d]\n",rlx_var.rlxd_level,rlx_var.rlxd_seed+MPI_PID);
		rlxd_init(rlx_var.rlxd_level,rlx_var.rlxd_seed+MPI_PID); /* use unique MPI_PID to shift seeds */
	}

	//                                    
	//lprintf("MAIN",0,"RLXD [%d,%d]\n",glb_var.rlxd_level,glb_var.rlxd_seed);
	//rlxd_init(glb_var.rlxd_level,glb_var.rlxd_seed+PID);
	//srand(glb_var.rlxd_seed+PID);

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
			printf(" %3.3e \n",m[nm]);
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
	{
		lprintf("MAIN",0,"Mass[%d] = %f\n",k,m[k]);
		lprintf("CORR",0,"Mass[%d] = %f\n",k,m[k]);
	}
	/* if a propagator and a source are provided , then read them and perform contractions [ debug only ] */
	if(strcmp(prop_filename,"")!=0 || strcmp(source_filename,"")!=0) {



	}
	else {
		list=NULL;
		if(strcmp(list_filename,"")!=0) {
			error((list=fopen(list_filename,"r"))==NULL,1,"main [meson_scattering.c]" ,
					"Failed to open list file\n");
		}


		i=0;
		lprintf("MAIN",0,"tsrc set to %i\n",mes_var.tsrc);
		while(1) {

			if(list!=NULL)
				if(fscanf(list,"%s",cnfg_filename)==0 || feof(list)) break;


			i++;

			lprintf("MAIN",0,"Configuration from %s\n", cnfg_filename);
			/* NESSUN CHECK SULLA CONSISTENZA CON I PARAMETRI DEFINITI !!! */
			read_gauge_field(cnfg_filename);
			represent_gauge_field();

			lprintf("TEST",0,"<p> %1.6f\n",avr_plaquette());
			full_plaquette();

			//								lprintf("CORR",0,"Number of noise vector : nhits = %i \n", mes_var.nhits);
			measure_pion_scattering( m, mes_var.nhits,i,  mes_var.precision,mes_var.tsrc);
			if(list==NULL) break;
		}

		if(list!=NULL) fclose(list);

	}




	finalize_process();



	free_BCs();

	free_gfield(u_gauge);
#ifdef ALLOCATE_REPR_GAUGE_FIELD
	free_gfield_f(u_gauge_f);
#endif

	/* close communications */
	finalize_process();

	gettimeofday(&end,0);
	timeval_subtract(&etime,&end,&start);
	lprintf("TIMING",0,"Inversions and contractions for configuration  [%s] done [%ld sec %ld usec]\n",cnfg_filename,etime.tv_sec,etime.tv_usec);




	return 0;
}
