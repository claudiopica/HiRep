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
  int ti;
  int tf;
  int ff_fixed_point;
  int dt;
  /* for the reading function */
  input_record_t read[8];
} input_mesons;

#define init_input_mesons(varname) \
{ \
  .read={\
    {"quark quenched masses", "mes:masses = %s", STRING_T, (varname).mstring}, \
    {"inverter precision", "mes:precision = %lf", DOUBLE_T, &(varname).precision},\
    {"t_initial", "mes:ti = %d", INT_T, &(varname).ti}, \
    {"t_final", "mes:tf = %d", INT_T, &(varname).tf}, \
    {"enable Dirichlet point", "mes:ff_dirichlet_point = %d", INT_T, &(varname).ff_fixed_point}, \
    {"Distance of t_initial from Dirichlet boundary", "mes:dirichlet_dt = %d", INT_T, &(varname).dt}, \
    {NULL, NULL, INT_T, NULL}				\
   }							\
}


char cnfg_filename[256]="";
char list_filename[256]="";
char input_filename[256] = "input_file";
char output_filename[256] = "mesons.out";
enum { UNKNOWN_CNFG, DYNAMICAL_CNFG, QUENCHED_CNFG };

input_mesons mes_var = init_input_mesons(mes_var);


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

//Copied From meson_measurements
static void fix_T_bc(int tau){
  int index;
  int ix,iy,iz;
  suNf *u;
  if (--tau<0) tau+= GLB_T;
  lprintf("meson_measurements",15,"Setting Dirichlet boundary conidtion at global time slice %d, %d\n",tau,T_BORDER);
  if((zerocoord[0]-1<=tau && zerocoord[0]+T>tau) || (zerocoord[0]==0 && tau==GLB_T-1)) { 
    for (ix=0;ix<X_EXT;++ix) for (iy=0;iy<Y_EXT;++iy) for (iz=0;iz<Z_EXT;++iz){
	  if( ( (tau==zerocoord[0]-1) || (zerocoord[0]==0 && tau==GLB_T-1)) && (NP_T>1) ){
	//printf("PID = %d, zc = %d, tau = %d\n", PID, zerocoord[0], tau);
	    index=ipt_ext(0,ix,iy,iz);
	  }
	  else{
	    index=ipt_ext(T_BORDER+tau-zerocoord[0],ix,iy,iz); 
	  }
	  if(index!=-1) {
	    u=pu_gauge_f(index,0);
	    _suNf_zero(*u);
	  }
	}
  }
  lprintf("meson_measurements",50,"Boundaries set!\n");
}

//Why all the curly brackets? This is how mathematica likes to get matrices.
static void print_prop(suNf_propagator S){
      int i,j;
	lprintf("PROP", 10, "{");
      for(i=0;i<4*NF;i++){ 
	lprintf("PROP", 10, "{");
      for(j=0;j<4*NF;j++){ 
	 if(j<4*NF-1){ lprintf("PROP", 10, "%.10f + I*%.10f, ", _PROP_IDX(S,i,j).re, _PROP_IDX(S,i,j).im ); }
	 else{ lprintf("PROP", 10, "%.10f + I*%.10f", _PROP_IDX(S,i,j).re, _PROP_IDX(S,i,j).im ); }
      }
      	if(i<4*NF-1){ lprintf("PROP", 10, "},"); }
	else{ lprintf("PROP", 10, "}"); }
      }
      lprintf("PROP", 10, "};\n");
}

//Check: g5 D^dag(x,0) g5 = D(0,x)
static void check_g5herm(spinor_field *prop1, int t1, spinor_field *prop2, int t2){

  int beta, a, ix1, ix2;

  suNf_propagator sp1,sp2,spdag;
  complex tr;
		lprintf("CK_G5HERM",0,"Only Works in serial!\n");
		ix1 = ipt(t1-zerocoord[0],0,0,0);
		ix2 = ipt(t2-zerocoord[0],0,0,0);
		for (a=0;a<NF;++a){
		    for (beta=0;beta<4;beta++){ 
		      _propagator_assign(sp1, *_FIELD_AT(&prop1[a*4+beta],ix1),a,beta); //S( (2,0,0,0), (0,0,0,0) )
		      _propagator_assign(sp2, *_FIELD_AT(&prop2[a*4+beta],ix2),a,beta); //S( (0,0,0,0), (2,0,0,0) )
		    }
		}
		_propagator_dagger(spdag,sp2);
		_g5_propagator(sp2,spdag); _propagator_g5(spdag,sp2);
		lprintf("CK_G5HERM",0,"Propagator1\n");
		print_prop(sp1);
		lprintf("CK_G5HERM",0,"g5 Propagator2^dagger g5 \n");
		print_prop(spdag);
		_propagator_sub(sp2,sp1,spdag);
		lprintf("CK_G5HERM",0,"Propagator1 - g5 Propagator2^dagger g5 \n");
		print_prop(sp2);
		_propagator_trace(tr,sp2);
		lprintf("CK_G5HERM",0,"Tr[ g5 Propagator1^dag g5 - Propagator2 ] = %g + I%g\n", tr.re, tr.im);
     
}

//source = g5 prop( x, 0 ) delta( x, (tf,0,0,0) )
void create_sequential_source_point(spinor_field *source, int tf, spinor_field* prop){

  int beta, a, ix;

  suNf_propagator sp0,sp1;

   for (beta=0;beta<4*NF;++beta){
     spinor_field_zero_f(&source[beta]);
   }

	ix = ipt(tf-zerocoord[0],0,0,0);
	for (a=0;a<NF;++a){
	    for (beta=0;beta<4;beta++){ 
	      _propagator_assign(sp0, *_FIELD_AT(&prop[a*4+beta],ix),a,beta);
	    }
	}
	_g5_propagator(sp1,sp0); //g5 Prop
	_propagator_transpose(sp0,sp1);

	for (a=0;a<NF;++a){
	    for (beta=0;beta<4;beta++){ 
	      *_FIELD_AT(&source[a*4+beta],ix) = sp0.c[a].c[beta];
	    }
	}

  for (beta=0;beta<4*NF;++beta){
     start_sf_sendrecv(source + beta);
     complete_sf_sendrecv(source + beta);
  }

}
//Check sequential Gamma Seq(0,0) = Gamma g5 S^dag(x,0) g5 g5 S(x,0)
static void check_sequential_point(spinor_field *prop_1, spinor_field *prop_2, spinor_field *prop_seq, int ti){

  lprintf("CK_SEQ",0,"Only Works in serial!\n");
  int ix1 = ipt(0,0,0,0);
  int ix2 = ipt(ti,0,0,0);
  suNf_propagator sp1,sp2,sp3,sptmp1, sptmp2;
  int a, beta;
  complex tr;

	for (a=0;a<NF;++a){
	    for (beta=0;beta<4;beta++){ 
	      _propagator_assign(sp1, *_FIELD_AT(&prop_seq[a*4+beta],ix1),a,beta); //S( (0,0,0,0), (2,0,0,0) ) g5 S( (2,0,0,0), (0,0,0,0) )
	      _propagator_assign(sp2, *_FIELD_AT(&prop_2[a*4+beta],ix1),a,beta);   //S( (0,0,0,0), (2,0,0,0) )
	      _propagator_assign(sp3, *_FIELD_AT(&prop_1[a*4+beta],ix2),a,beta);   //S( (2,0,0,0), (0,0,0,0) )
	    }
	}

	_g5_propagator(sptmp1,sp3);
	_propagator_mul(sptmp2,sp2,sptmp1);

	_propagator_trace(tr,sptmp2);
	lprintf("CK_SEQ",0,"S(0,x) g5 S(x,0) point = %g + I%g\n", tr.re, tr.im);
	print_prop(sptmp2);

	_propagator_trace(tr,sp1);
	lprintf("CK_SEQ",0,"S(0,x) g5 S(x,0) seq   = %g + I%g\n", tr.re, tr.im);
	print_prop(sp1);

	_propagator_sub(sptmp1,sp1,sptmp2);
	_propagator_trace(tr,sptmp1);
	lprintf("CK_SEQ",0,"point - seq   = %g + I%g\n", tr.re, tr.im);
	print_prop(sptmp1);
}

static void check_sequential(spinor_field *prop_seq, spinor_field *prop_1, int tau){

lprintf("CK_SEQ",0,"Only Works in serial!\n");

double Corr[2][GLB_T];
int ix,t,x,y,z,a,beta,tc;
complex tr;
suNf_propagator sp0,sp1,spdag, sptmp;

int i;
for (i=0;i<2;i++){
for (t=0; t<GLB_T; t++) { 
	Corr[i][t] = 0; 
}}


	for (t=0; t<T; t++) {
	    tc = (zerocoord[0]+t+GLB_T-tau)%GLB_T;	 
	    for (x=0; x<X; x++) for (y=0; y<Y; y++) for (z=0; z<Z; z++) { 
		  ix=ipt(t,x,y,z);					
		  
		  for (a=0;a<NF;++a){
		    for (beta=0;beta<4;beta++){ 
		      _propagator_assign(sp0, *_FIELD_AT(&prop_seq[a*4+beta],ix),a,beta);
		      _propagator_assign(sp1, *_FIELD_AT(&prop_1[a*4+beta],ix),a,beta);
		    }
		  }
                  _propagator_dagger(spdag,sp1);

			//Pion
			_propagator_mul(sptmp,sp1,spdag); _propagator_trace(tr, sptmp); 
			Corr[0][tc] += tr.re;

			//g5 Seq(0,0)
			if(t == 0 && x==0 && y==0 && z==0){
				_g5_propagator(sptmp,sp0); 
				_propagator_trace(tr, sptmp);
				Corr[1][tc] += tr.re;
			}
		} //END SPATIAL LOOP
	  } //END TIME LOOP

	lprintf("CORR",0,"Pion: "); for (t=0; t<T; t++){ lprintf("CORR",0,"%g ",Corr[0][t]/(GLB_VOL3) ); } lprintf("CORR",0,"\n");
	lprintf("CORR",0,"SeqPion: "); for (t=0; t<T; t++){ lprintf("CORR",0,"%g ",Corr[1][t]/(GLB_VOL3) ); } lprintf("CORR",0,"\n");

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
  logger_setlevel(0,30);
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

  i=0;

  spinor_field* source;
  spinor_field* prop_1;
  spinor_field* prop_2;
  spinor_field* source_seq;
  spinor_field* prop_seq;
  suNf_field* u_gauge_old=alloc_gfield_f(&glattice);

  source = alloc_spinor_field_f(4,&glattice);
  source_seq = alloc_spinor_field_f(4*NF,&glattice);
  prop_1 = alloc_spinor_field_f(4*NF,&glattice);
  prop_2 = alloc_spinor_field_f(4*NF,&glattice);
  prop_seq = alloc_spinor_field_f(4*NF,&glattice);
  
  while(1) {
    struct timeval start, end, etime;

    if(list!=NULL)
      if(fscanf(list,"%s",cnfg_filename)==0 || feof(list)) break;

    i++;

    lprintf("MAIN",0,"Configuration from %s\n", cnfg_filename);
    read_gauge_field(cnfg_filename);
    //unit_gauge(u_gauge);
    represent_gauge_field();

    lprintf("TEST",0,"<p> %1.6f\n",avr_plaquette());

    full_plaquette();
    gettimeofday(&start,0);

  int k;

  suNf_field_copy(u_gauge_old,u_gauge_f);//Save the gaugefield
  if(mes_var.ff_fixed_point){
    lprintf("MAIN",10,"Applying Dirichlet Boundaries\n");
    fix_T_bc(mes_var.ti-mes_var.dt);//Apply fixed boundaryconditions by zeroing links at time slice tau to direction 0.
  } else {
    lprintf("MAIN",10,"Default Boundaries\n");
  }

  init_propagator_eo(1,m, mes_var.precision);//1 for number of masses 
  for (k=0;k<NF;++k){
    create_point_source(source,0,k);
    calc_propagator(prop_1+4*k,source,4);//4 for spin components
    create_point_source(source,mes_var.ti,k);
    calc_propagator(prop_2+4*k,source,4);//4 for spin components
  }
  check_g5herm(prop_1, mes_var.ti, prop_2, 0);

  create_sequential_source_point(source_seq, mes_var.ti, prop_1);
  calc_propagator(prop_seq,source_seq,4*NF);
  check_sequential_point(prop_1, prop_2, prop_seq, mes_var.ti); 

  create_sequential_source(source_seq,mes_var.tf,prop_1); 
  calc_propagator(prop_seq,source_seq,4*NF);
  check_sequential(prop_seq,prop_1,0);

  suNf_field_copy(u_gauge_f,u_gauge_old);//Restore the gaugefield
  free_propagator_eo(); 

    gettimeofday(&end,0);
    timeval_subtract(&etime,&end,&start);
    lprintf("MAIN",0,"Configuration #%d: analysed in [%ld sec %ld usec]\n",i,etime.tv_sec,etime.tv_usec);
    if(list==NULL) break;
  }

  free_gfield_f(u_gauge_old);
  free_spinor_field_f(source);  
  free_spinor_field_f(source_seq);
  free_spinor_field_f(prop_1);
  free_spinor_field_f(prop_2);
  free_spinor_field_f(prop_seq);

  if(list!=NULL) fclose(list);


  free_BCs();

  free_gfield(u_gauge);
#ifdef ALLOCATE_REPR_GAUGE_FIELD
  free_gfield_f(u_gauge_f);
#endif

  finalize_process();

  return 0;
}

