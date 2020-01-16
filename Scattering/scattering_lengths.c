/*******************************************************************************
 *
 * Computation of the Pseudo scalar scattering lengths
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
#include "setup.h"
#include "clover_tools.h"
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
	char mstring[256],configlist[256],outpath[256];
    double csw;	
	double precision;
	int nhits;

	/* for the reading function */
	input_record_t read[8];

} input_scatt;

#define init_input_scatt(varname) \
{ \
	.read={\
		{"quark quenched masses", "mes:masses = %s", STRING_T, (varname).mstring},\
		{"csw", "mes:csw = %lf", DOUBLE_T, &(varname).csw},\
		{"inverter precision", "mes:precision = %lf", DOUBLE_T, &(varname).precision},\
		{"number of inversions per cnfg", "mes:nhits = %d", INT_T, &(varname).nhits},\
		{"Configuration list:", "mes:configlist = %s", STRING_T, &(varname).configlist},\
		{"outpath:", "mes:outpath = %s", STRING_T, &(varname).outpath},\
		{NULL, NULL, INT_T, NULL}\
	}\
}


char cnfg_filename[256]="";
char list_filename[256]="";
char source_filename[256]="";
char input_filename[256] = "input_file";
char output_filename[256] = "meson_scattering.out";
char output_dir[256] = "./output/";
int Nsource;
double M;

enum { UNKNOWN_CNFG, DYNAMICAL_CNFG, QUENCHED_CNFG };

input_scatt mes_var = init_input_scatt(mes_var);

typedef struct {
	char string[256];
	char configlist[256];
	char outpath[256];
	int t, x, y, z;
	int nc, nf;
	double b, m;
	int n;
	int type;
} filename_t;


//to move




/* Random timeslice not previously chosen */
static int random_tau(){
  static int* slices=NULL;
  if (slices == NULL) slices = (int*) malloc(GLB_T*sizeof(int));
  static int counter = 0;
  int itmp,tau,i;
  double ran;

  if (counter == 0){
    for (i=0;i<GLB_T;++i){
      slices[i]=i;
    }
    counter=GLB_T;
  }
  do{
    ranlxd(&ran,1);
    itmp=(int)(ran*counter);
  } while(itmp==counter);
  counter--;
  tau = slices[itmp];
  slices[itmp]=slices[counter];
  slices[counter]=tau;
  bcast_int(&tau,1);
  return tau;
}

static int gi(int i){
	if(i==1) return _g1;
	if(i==2) return _g2;
	if(i==3) return _g3;

	return -1;
}
static void do_global_sum(meson_observable *mo, double norm)
{
  meson_observable *motmp = mo;
  int i;
  while (motmp != NULL)
  {
    global_sum(motmp->corr_re, motmp->corr_size);
    global_sum(motmp->corr_im, motmp->corr_size);
    for (i = 0; i < motmp->corr_size; i++)
    {
      motmp->corr_re[i] *= norm;
      motmp->corr_im[i] *= norm;
    }
    motmp = motmp->next;
  }
}
void measure_pion_scattering_I2(double* m, int numsources, double precision,char* path,char* cnfg_filename){
	int ts;
	meson_observable *rho1[3][3],*rho2[3][3];
	meson_observable *pi1,*pi2;
	meson_observable *AD;
	meson_observable *BC;
	spinor_field* source_ts1 = alloc_spinor_field_f(4,&glattice);
	spinor_field* source_ts2 = alloc_spinor_field_f(4,&glattice);
	char auxname[256];

	spinor_field* prop_ts1 =  alloc_spinor_field_f(4 ,&glattice);
	spinor_field* prop_ts2=  alloc_spinor_field_f(4 ,&glattice);

	pi1 = (meson_observable*) malloc(sizeof(meson_observable));
	pi2 = (meson_observable*) malloc(sizeof(meson_observable));
	AD = (meson_observable*) malloc(sizeof(meson_observable));
	BC = (meson_observable*) malloc(sizeof(meson_observable));

	for(int i=0; i<3; i++){
		for(int j=0; j<3; j++){
			rho1[i][j] = (meson_observable*) malloc(sizeof(meson_observable));
			rho2[i][j] = (meson_observable*) malloc(sizeof(meson_observable));
		}
	}
	init_mo(pi1,"Pi1",GLB_T);
	init_mo(pi2,"Pi2",GLB_T);
	init_mo(AD,"AD",GLB_T);
	init_mo(BC,"BC",GLB_T);

	pi1->ind1 = _g5;
	pi1->ind2 = _g5;
	pi2->ind1 = _g5;
	pi2->ind2 = _g5;
		
	for(int i=0; i<3; i++){
		for(int j=0; j<3; j++){
				init_mo(rho1[i][j],"rho1",GLB_T);
				init_mo(rho2[i][j],"rho2",GLB_T);
				rho1[i][j]->ind1 = gi(i);
				rho1[i][j]->ind2 = gi(j);
				rho2[i][j]->ind1 = gi(i);
				rho2[i][j]->ind2 = gi(j);
		}
	}

	for (int src=0;src<numsources;++src)
   	{
		reset_mo(pi1);
		reset_mo(pi2);
		reset_mo(AD);
		reset_mo(BC);
		for(int i=0; i<3; i++){
			for(int j=0; j<3; j++){
				reset_mo(rho1[i][j]);
				reset_mo(rho2[i][j]);
		}
	}
	
	init_propagator_eo(1, m, precision);
	ts=random_tau();
	spinor_field_zero_f(source_ts1);
	create_diluted_source_equal_atau(source_ts1, ts);
	calc_propagator(prop_ts1,source_ts1,4);
	spinor_field_zero_f(source_ts2);
	create_diluted_source_equal_atau(source_ts2, ts);
	calc_propagator(prop_ts2,source_ts2,4);
	lprintf("MAIN",0,"Start to perform the contractions ...\n");
	
	// "standard" two points : pi and rho 
	measure_mesons_core(prop_ts1,prop_ts1,source_ts1,pi1,1,ts,1,0,GLB_T);
	measure_mesons_core(prop_ts2,prop_ts2,source_ts2,pi2,1,ts,1,0,GLB_T);
	do_global_sum(pi1,1.0);
	do_global_sum(pi2,1.0);

	for(int i=0; i<3; i++){
			for(int j=0; j<3; j++){
				measure_mesons_core(prop_ts1,prop_ts1,source_ts1,rho1[i][j],1,ts,1,0,GLB_T);
				measure_mesons_core(prop_ts2,prop_ts2,source_ts2,rho2[i][j],1,ts,1,0,GLB_T);
				do_global_sum(rho1[i][j],1.0);
				do_global_sum(rho2[i][j],1.0);

			}
	}

	// contraction 4 particles two-points.
	measure_scattering_AD_core(AD, prop_ts1,prop_ts1,prop_ts2,prop_ts2, ts, 0,0,0,0,0); 
	measure_scattering_BC_core(BC, prop_ts1,prop_ts1,prop_ts2,prop_ts2, ts, 0,0,0,0,0);
	
	lprintf("MAIN",0,"Contraction done\n");
	// Printing.
	io2pt(pi1,1,src,path,"pi1",cnfg_filename);
	io2pt(pi2,1,src,path,"pi2",cnfg_filename);
	for(int i=0; i<3; i++){	
		for(int j=0; j<3; j++){
			sprintf(auxname, "rho1_%d%d",i,j);
			io2pt(rho1[i][j],1,src,path,auxname,cnfg_filename);
			sprintf(auxname, "rho2_%d%d",i,j);
			io2pt(rho2[i][j],1,src,path,auxname,cnfg_filename);
		}
	}
	io4pt(AD,0,src,path,"AD",cnfg_filename);
	io4pt(BC,0,src,path,"BC",cnfg_filename);
	
	
	}

	//free memory
  	free_spinor_field_f(source_ts1);
  	free_spinor_field_f(source_ts2);
	free_spinor_field_f(prop_ts1);
  	free_spinor_field_f(prop_ts2);
	free_mo(pi1);
	free_mo(pi2);
	free_mo(AD);
	free_mo(BC);
	for(int i=0; i<3; i++){
		for(int j=0; j<3; j++){
			free_mo(rho1[i][j]);
			free_mo(rho2[i][j]);
		}
	}

}

// end : to move








int main(int argc,char *argv[]) {
	int i,k;
	FILE* list;
	int nm;
	double m[256];

 	/* setup process communications */
  	setup_process(&argc,&argv);

	setup_gauge_fields();

	read_input(glb_var.read,get_input_filename());
	read_input(mes_var.read,get_input_filename());
    read_input(rlx_var.read,get_input_filename());

	strcpy(list_filename,mes_var.configlist);
	strcpy(output_dir,mes_var.outpath);

	lprintf("MAIN",0,"list_filename = %s %s\n", list_filename,mes_var.configlist);	
	lprintf("MAIN",0,"output_dir : %s \n",mes_var.outpath);	

  	if(strcmp(list_filename,"")!=0) {
    	error((list=fopen(list_filename,"r"))==NULL,1,"main [mk_mesons.c]" ,
		"Failed to open list file\n");
  	}

  	#ifdef WITH_CLOVER
		clover_init(mes_var.csw);//  VD: not nice here.		
  	#endif
	

	nm=0;
	m[0] = -atof(mes_var.mstring); // 	
  	init_propagator_eo(1,m,mes_var.precision);
	

	lprintf("MAIN",0,"Inverter precision = %e\n",mes_var.precision);
	for(k=0;k<nm;k++)
	{
		lprintf("MAIN",0,"Mass[%d] = %f\n",k,m[k]);
		lprintf("CORR",0,"Mass[%d] = %f\n",k,m[k]);
	}
	/* if a propagator and a source are provided , then read them and perform contractions [ debug only ] */
	

	i=0;
	lprintf("CORR",0,"Number of noise vector : nhits = %i \n", mes_var.nhits);
	while(1){
    	struct timeval start, end, etime;
	    gettimeofday(&start,0);
    	if(list!=NULL)
      		if(fscanf(list,"%s",cnfg_filename)==0 || feof(list)) break;

    	lprintf("MAIN",0,"Configuration from %s\n", cnfg_filename);
    	read_gauge_field(cnfg_filename);
    	represent_gauge_field();
		lprintf("TEST",0,"<p> %1.6f\n",avr_plaquette());
		full_plaquette();	
		measure_pion_scattering_I2( m, mes_var.nhits,  mes_var.precision,mes_var.outpath,cnfg_filename);
		gettimeofday(&end,0);
		timeval_subtract(&etime,&end,&start);
	    lprintf("MAIN",0,"Configuration : analysed in [%ld sec %ld usec]\n",etime.tv_sec,etime.tv_usec);

    	if(list==NULL) break;
	}
  lprintf("DEBUG",0,"ALL done, deallocating\n");

  if(list!=NULL) fclose(list);
  
  finalize_process();	
  return 0;

}
