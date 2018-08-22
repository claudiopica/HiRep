// This code will contain all the contractions necessary for rho to pi pi scattering
#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "global.h"
#include "io.h"
#include "random.h"
#include "error.h"
#include "geometry.h"
#include "memory.h"
#include "statistics.h"
#include "update.h"
#include "scattering.h"
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

#include "cinfo.c"
#include "IOroutines.c"

#define PI 3.141592653589793238462643383279502884197

// Function for initiating meson observable (used to store the correlation function)
void init_mo(meson_observable* mo, char* name, int size)
{
  int i;
  //ind1 and ind2 don't do anything for the moment
  mo->ind1 = _g5;
  mo->ind2 = _g5;
  strcpy(mo->channel_name,name);
  strcpy(mo->channel_type, "Pi Pi scattering");
  mo->sign=1.0;
  mo->corr_size=size;
  mo->corr_re = (double * ) malloc(size * sizeof(double));
  mo->corr_im = (double * ) malloc(size * sizeof(double));
  if (mo->corr_re == NULL || mo->corr_im == NULL)
  {
    fprintf(stderr, "malloc failed in init_mo \n");
    return;
  }
  mo->next=NULL;
  for (i=0; i<size; ++i)
  {
    mo->corr_re[i]=0.0;
    mo->corr_im[i]=0.0;
  }
}

void reset_mo(meson_observable* mo)
{
  int i;
  for (i=0; i< mo->corr_size; ++i)
  {
    mo->corr_re[i]=0.0;
    mo->corr_im[i]=0.0;
  }
}

static void do_global_sum(meson_observable* mo, double norm){
  meson_observable* motmp=mo;
  int i;
  while (motmp!=NULL){
      global_sum(motmp->corr_re,motmp->corr_size);
      global_sum(motmp->corr_im,motmp->corr_size);
      for(i=0; i<motmp->corr_size; i++){
	motmp->corr_re[i] *= norm;
	motmp->corr_im[i] *= norm;
      }
    motmp=motmp->next;
  }
}

void free_mo(meson_observable* mo)
{
  free(mo->corr_re);
  free(mo->corr_im);
  free(mo);
}

#define BASENAME(filename) (strrchr((filename),'/') ? strrchr((filename),'/')+1 : filename )

#define corr_ind(px,py,pz,n_mom,tc,nm,cm) ((px)*(n_mom)*(n_mom)*GLB_T*(nm)+(py)*(n_mom)*GLB_T*(nm)+(pz)*GLB_T*(nm)+ ((cm)*GLB_T) +(tc))
inline void io2pt(meson_observable* mo, int pmax, int sourceno, char* path, char* name)
{
	FILE* file;
	char outfile[256] = {};
	int px,py,pz,t;
	if(PID==0){
		sprintf(outfile,"%s/%s_src_%d_%s", path, name, sourceno, BASENAME(cnfg_filename) );
		file=fopen(outfile,"w+");
		//Factor of 2 to correct for the noise source normalisation
		for(px=0;px<pmax;++px) for(py=0;py<pmax;++py) for(pz=0;pz<pmax;++pz) for(t=0;t<GLB_T;++t) fprintf(file,"%i %i %i %i %3.10e %3.10e \n", px, py, pz, t,2*(mo->corr_re[corr_ind(px,py,pz,pmax,t,1,0)]), 2*(mo->corr_im[corr_ind(px,py,pz,pmax,t,1,0)]));
		fclose(file);
	}
	return;
}

#define INDEX(px,py,pz,n_mom,tc) ((px + n_mom)*(2*n_mom+1)*(2*n_mom+1)*(GLB_T)+(py + n_mom)*(2*n_mom+1)*(GLB_T)+(pz + n_mom)*(GLB_T)+ (tc))
inline void io4pt(meson_observable* mo, int pmax, int sourceno, char* path, char* name)
{
	FILE* file;
	char outfile[256] = {};
	int px,py,pz,t;
	if(PID==0){
		sprintf(outfile,"%s/%s_src_%d_%s", path, name, sourceno, BASENAME(cnfg_filename) );
		file=fopen(outfile,"w+");
		//Factor of 4 to correct for the noise source normalisation
		for(px=-pmax;px<=pmax;++px) for(py=-pmax;py<=pmax;++py) for(pz=-pmax;pz<=pmax;++pz) for(t=0;t<GLB_T;++t) fprintf(file, "%i %i %i %i %3.10e %3.10e \n", px, py, pz, t,4*(mo->corr_re[INDEX(px,py,pz, pmax,t)]), 4*(mo->corr_im[INDEX(px,py,pz,pmax,t)]));
		fclose(file);
	}
	return;
}

int main(int argc,char *argv[])
{
  int src,t;
  int px,py,pz, px2, py2, pz2;
  int tau=0;
  filename_t fpars;
  int nm;
  char tmp[256], *cptr;
  int i,k;
  double m[256];
  FILE* list;

  //Copy I/O from another file
  read_cmdline(argc, argv);
  setup_process(&argc,&argv);

  read_input(glb_var.read,input_filename);
  setup_replicas();

  /* logger setup */
  /* disable logger for MPI processes != 0 */
  //logger_setlevel(0,500);
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

  list=NULL;
  if(strcmp(list_filename,"")!=0) {
    error((list=fopen(list_filename,"r"))==NULL,1,"main [mk_mesons.c]" ,
	"Failed to open list file\n");
  }

  /* read & broadcast parameters */
  parse_cnfg_filename(cnfg_filename,&fpars);

  read_input(mes_var.read,input_filename);
  char *path=mes_var.outdir;
  printf("The momenta are (%s) and (%s) \n", mes_var.p1, mes_var.p2);
  sscanf(mes_var.p1, "(%d,%d,%d)" ,&px,&py,&pz);
  sscanf(mes_var.p2, "(%d,%d,%d)" ,&px2,&py2,&pz2);
  printf("The momenta are (%d %d %d) and (%d %d %d) \n", px, py, pz, px2, py2, pz2);
  int numsources = mes_var.nhits;
  GLB_T=fpars.t; GLB_X=fpars.x; GLB_Y=fpars.y; GLB_Z=fpars.z;

  /* setup lattice geometry */
  if (geometry_init() == 1) { finalize_process(); return 0; }
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

  read_gauge_field(cnfg_filename);
  represent_gauge_field();
  //End of the I/O block

  init_propagator_eo(nm,m,mes_var.precision);

#define OBSERVABLE_LIST X(pi2p) X(pi_p1) X(pi_p2) X(twopt_nomom_g1) X(twopt_nomom_g2) X(twopt_nomom_g3) X(twopt_g1) X(twopt_g2) X(twopt_g3) X(d_nomom) X(d) X(r1) X(r2) X(r3) X(r4) X(t1_g1) X(t1_g2) X(t1_g3) X(t2_g1) X(t2_g2) X(t2_g3) X(twoptp2_g1) X(twoptp2_g2) X(twoptp2_g3) X(dp2) X(r1p2) X(r2p2) X(r3p2) X(r4p2) X(t1p2_g1) X(t1p2_g2) X(t1p2_g3) X(t2p2_g1) X(t2p2_g2) X(t2p2_g3)
#define TMP_OBSERVABLE_LIST X(r1) X(r2) X(r3) X(r4) X(r1p2) X(r2p2) X(r3p2) X(r4p2)
#define X(NAME) meson_observable* mo_##NAME = malloc(sizeof(meson_observable)); init_mo(mo_##NAME,#NAME,27*GLB_T);
OBSERVABLE_LIST
#undef X
#define X(NAME) meson_observable* mo_tmp##NAME = malloc(sizeof(meson_observable)); init_mo(mo_tmp##NAME,#NAME,27*GLB_T);
TMP_OBSERVABLE_LIST
#undef X
mo_twopt_nomom_g1->ind1=_g1;
mo_twopt_nomom_g1->ind2=_g1;
mo_twopt_g1->ind1=_g1;
mo_twopt_g1->ind2=_g1;
mo_t1_g1->ind2=_g1;
mo_t2_g1->ind2=_g1;
mo_twoptp2_g1->ind1=_g1;
mo_twoptp2_g1->ind2=_g1;
mo_t1p2_g1->ind2=_g1;
mo_t2p2_g1->ind2=_g1;

mo_twopt_nomom_g2->ind1=_g2;
mo_twopt_nomom_g2->ind2=_g2;
mo_twopt_g2->ind1=_g2;
mo_twopt_g2->ind2=_g2;
mo_t1_g2->ind2=_g2;
mo_t2_g2->ind2=_g2;
mo_twoptp2_g2->ind1=_g2;
mo_twoptp2_g2->ind2=_g2;
mo_t1p2_g2->ind2=_g2;
mo_t2p2_g2->ind2=_g2;

mo_twopt_nomom_g3->ind1=_g3;
mo_twopt_nomom_g3->ind2=_g3;
mo_twopt_g3->ind1=_g3;
mo_twopt_g3->ind2=_g3;
mo_t1_g3->ind2=_g3;
mo_t2_g3->ind2=_g3;
mo_twoptp2_g3->ind1=_g3;
mo_twoptp2_g3->ind2=_g3;
mo_t1p2_g3->ind2=_g3;
mo_t2p2_g3->ind2=_g3;

#define SOURCE_LIST_Q X(0) X(0_eta) X(p) X(mp) X(p2) X(mp2) 
#define SOURCE_LIST_W X(0_p) X(0_mp) X(p_0) X(mp_0) X(0_p2) X(0_mp2) X(p2_0) X(mp2_0)
#define X(NAME) lprintf("DEBUG",0,"Creating source source_%s \n", #NAME);spinor_field* source_##NAME = alloc_spinor_field_f(4*NF,&glattice); for(i=0;i<4*NF;++i) {spinor_field_zero_f(&source_##NAME[i]);}
    SOURCE_LIST_Q
    SOURCE_LIST_W
    spinor_field* source_0_0 = alloc_spinor_field_f(4*NF,&glattice); spinor_field_zero_f(source_0_0);
#undef X
#define X(NAME) lprintf("DEBUG",0,"Creating propagator Q_%s \n", #NAME); spinor_field* Q_##NAME = alloc_spinor_field_f(4*NF,&glattice);
    SOURCE_LIST_Q
#undef X
#define X(NAME) lprintf("DEBUG",0,"Creating propagator W_%s \n", #NAME);spinor_field* W_##NAME = alloc_spinor_field_f(4*NF,&glattice);
    SOURCE_LIST_W
#undef X
    spinor_field* W_0_0 = alloc_spinor_field_f(4*NF,&glattice);


while(1){
    if(list!=NULL)
      if(fscanf(list,"%s",cnfg_filename)==0 || feof(list)) break;

    lprintf("MAIN",0,"Configuration from %s\n", cnfg_filename);
    read_gauge_field(cnfg_filename);
    represent_gauge_field();

    for (src=0;src<numsources;++src)
    {
	    // Clear all meson observables
#define X(NAME) reset_mo(mo_##NAME);
	    OBSERVABLE_LIST
#undef X


    	// Sources for non-sequential props
  	create_diluted_source_equal_atau(source_0, tau);
  	create_diluted_source_equal_atau(source_0_eta, tau);
    	add_momentum(source_p, source_0, px, py, pz);
    	add_momentum(source_mp, source_0, -px, -py, -pz);
    	add_momentum(source_p2, source_0, px2, py2, pz2);
    	add_momentum(source_mp2, source_0, -px2, -py2, -pz2);
	lprintf("DEBUG",0,"Created sources \n");
	// Non-sequential prop inversions
#define X(NAME) lprintf("DEBUG",0,"Calculating propagator Q_%s \n", #NAME); calc_propagator(Q_##NAME, source_##NAME,4);
	SOURCE_LIST_Q
#undef X
	lprintf("DEBUG",0,"Created propagators \n");
	// Sequential sources
	create_sequential_source(source_0_0,tau,Q_0);
	create_sequential_source(source_p_0,tau,Q_p);
	create_sequential_source(source_mp_0,tau,Q_mp);
	add_momentum(source_0_p,source_0_0,px,py,pz);
	add_momentum(source_0_mp,source_0_0,-px,-py,-pz);

	create_sequential_source(source_p2_0,tau,Q_p2);
	create_sequential_source(source_mp2_0,tau,Q_mp2);
	add_momentum(source_0_p2,source_0_0,px2,py2,pz2);
	add_momentum(source_0_mp2,source_0_0,-px2,-py2,-pz2);

	lprintf("DEBUG",0,"Created seq sources \n");

	//Sequential propagators
#define X(NAME) calc_propagator(W_##NAME, source_##NAME,4);
	SOURCE_LIST_W
#undef X
	lprintf("DEBUG",0,"Created seq propagators \n");

	// pipi-> pipi direct (2 traces) and pipi-> rho contractions go here
	
	//2-point pi -> pi

	measure_mesons_core(Q_0, Q_0, source_0, mo_pi2p, 1, tau, 2, 0, GLB_T);
	do_global_sum(mo_pi2p,1.0);
	measure_mesons_core(Q_0, Q_p, source_0, mo_pi_p1, 1, tau, 2, 0, GLB_T);
	do_global_sum(mo_pi_p1,1.0);
	measure_mesons_core(Q_0, Q_p2, source_0, mo_pi_p2, 1, tau, 2, 0, GLB_T);
	do_global_sum(mo_pi_p2,1.0);

	// 2-point rho->rho
	measure_mesons_core(Q_0, Q_0, source_0, mo_twopt_nomom_g1, 1, tau, 2, 0, GLB_T);
	do_global_sum(mo_twopt_nomom_g1,1.0);
	measure_mesons_core(Q_0, Q_p, source_0, mo_twopt_g1, 1, tau, 2, 0, GLB_T);
	do_global_sum(mo_twopt_g1,1.0);
	measure_mesons_core(Q_0, Q_p2, source_0, mo_twoptp2_g1, 1, tau, 2, 0, GLB_T);
	do_global_sum(mo_twoptp2_g1,1.0);

	measure_mesons_core(Q_0, Q_0, source_0, mo_twopt_nomom_g2, 1, tau, 2, 0, GLB_T);
	do_global_sum(mo_twopt_nomom_g2,1.0);
	measure_mesons_core(Q_0, Q_p, source_0, mo_twopt_g2, 1, tau, 2, 0, GLB_T);
	do_global_sum(mo_twopt_g2,1.0);
	measure_mesons_core(Q_0, Q_p2, source_0, mo_twoptp2_g2, 1, tau, 2, 0, GLB_T);
	do_global_sum(mo_twoptp2_g2,1.0);

	measure_mesons_core(Q_0, Q_0, source_0, mo_twopt_nomom_g3, 1, tau, 2, 0, GLB_T);
	do_global_sum(mo_twopt_nomom_g3,1.0);
	measure_mesons_core(Q_0, Q_p, source_0, mo_twopt_g3, 1, tau, 2, 0, GLB_T);
	do_global_sum(mo_twopt_g3,1.0);
	measure_mesons_core(Q_0, Q_p2, source_0, mo_twoptp2_g3, 1, tau, 2, 0, GLB_T);
	do_global_sum(mo_twoptp2_g3,1.0);
	
	// direct 1
	measure_scattering_AD_core(mo_d_nomom, Q_0, Q_0, Q_0_eta, Q_0_eta, tau, 0, 1, 0, 0, 0 );
	measure_scattering_AD_core(mo_d, Q_p, Q_0, Q_0_eta, Q_0_eta, tau, 0, 1, px, py, pz );
	measure_scattering_AD_core(mo_dp2, Q_p2, Q_0, Q_0_eta, Q_0_eta, tau, 0, 1, px2, py2, pz2 );
	//Triangle pipi->rho
	measure_mesons_core(W_0_mp, Q_0, source_0, mo_t1_g1, 1, tau, 2, 0, GLB_T);
	do_global_sum(mo_t1_g1,1.0);
	measure_mesons_core(Q_0, W_0_p, source_0, mo_t2_g1, 1, tau, 2, 0, GLB_T);
	do_global_sum(mo_t2_g1,1.0);
	measure_mesons_core(W_0_mp2, Q_0, source_0, mo_t1p2_g1, 1, tau, 2, 0, GLB_T);
	do_global_sum(mo_t1p2_g1,1.0);
	measure_mesons_core(Q_0, W_0_p2, source_0, mo_t2p2_g1, 1, tau, 2, 0, GLB_T);
	do_global_sum(mo_t2p2_g1,1.0);

	measure_mesons_core(W_0_mp, Q_0, source_0, mo_t1_g2, 1, tau, 2, 0, GLB_T);
	do_global_sum(mo_t1_g2,1.0);
	measure_mesons_core(Q_0, W_0_p, source_0, mo_t2_g2, 1, tau, 2, 0, GLB_T);
	do_global_sum(mo_t2_g2,1.0);
	measure_mesons_core(W_0_mp2, Q_0, source_0, mo_t1p2_g2, 1, tau, 2, 0, GLB_T);
	do_global_sum(mo_t1p2_g2,1.0);
	measure_mesons_core(Q_0, W_0_p2, source_0, mo_t2p2_g2, 1, tau, 2, 0, GLB_T);
	do_global_sum(mo_t2p2_g2,1.0);

	measure_mesons_core(W_0_mp, Q_0, source_0, mo_t1_g3, 1, tau, 2, 0, GLB_T);
	do_global_sum(mo_t1_g3,1.0);
	measure_mesons_core(Q_0, W_0_p, source_0, mo_t2_g3, 1, tau, 2, 0, GLB_T);
	do_global_sum(mo_t2_g3,1.0);
	measure_mesons_core(W_0_mp2, Q_0, source_0, mo_t1p2_g3, 1, tau, 2, 0, GLB_T);
	do_global_sum(mo_t1p2_g3,1.0);
	measure_mesons_core(Q_0, W_0_p2, source_0, mo_t2p2_g3, 1, tau, 2, 0, GLB_T);
	do_global_sum(mo_t2p2_g3,1.0);

	//File IO
	io2pt(mo_pi2p, 2, src, path, "pi");
	io2pt(mo_pi_p1, 2, src, path, "pi_p1");
	io2pt(mo_pi_p2, 2, src, path, "pi_p2");
	io4pt(mo_d_nomom, 1, src, path, "d_p0");
	io4pt(mo_d, 1, src, path, "d");
	io4pt(mo_dp2, 1, src, path, "dp2");

	io2pt(mo_twopt_nomom_g1, 2, src, path, "rho_p0_g1");
	io2pt(mo_twopt_g1, 2, src, path, "rho_g1");
	io2pt(mo_t1_g1, 2, src, path, "t1_g1");
	io2pt(mo_t2_g1, 2, src, path, "t2_g1");
	io2pt(mo_twoptp2_g1, 2, src, path, "rhop2_g1");
	io2pt(mo_t1p2_g1, 2, src, path, "t1p2_g1");
	io2pt(mo_t2p2_g1, 2, src, path, "t2p2_g1");

	io2pt(mo_twopt_nomom_g2, 2, src, path, "rho_p0_g2");
	io2pt(mo_twopt_g2, 2, src, path, "rho_g2");
	io2pt(mo_t1_g2, 2, src, path, "t1_g2");
	io2pt(mo_t2_g2, 2, src, path, "t2_g2");
	io2pt(mo_twoptp2_g2, 2, src, path, "rhop2_g2");
	io2pt(mo_t1p2_g2, 2, src, path, "t1p2_g2");
	io2pt(mo_t2p2_g2, 2, src, path, "t2p2_g2");

	io2pt(mo_twopt_nomom_g3, 2, src, path, "rho_p0_g3");
	io2pt(mo_twopt_g3, 2, src, path, "rho_g3");
	io2pt(mo_t1_g3, 2, src, path, "t1_g3");
	io2pt(mo_t2_g3, 2, src, path, "t2_g3");
	io2pt(mo_twoptp2_g3, 2, src, path, "rhop2_g3");
	io2pt(mo_t1p2_g3, 2, src, path, "t1p2_g3");
	io2pt(mo_t2p2_g3, 2, src, path, "t2p2_g3");

	//Will have to change the range; this is too big
	for (t=0; t<GLB_T; ++t)
	{
#define X(NAME) reset_mo(mo_tmp##NAME);
		TMP_OBSERVABLE_LIST
#undef X
		//Invert the backward propagator
		create_sequential_source(source_0_0,t,Q_0);
		calc_propagator(W_0_0, source_0_0, 1);
		//Rectangle contractions and rho->pipi contractions go here
		// (Even though rho->pipi are not strictly necessary)
		measure_mesons_core(W_mp_0, W_0_0, source_0, mo_tmpr1, 1, tau, 2, 0, GLB_T);
		measure_mesons_core(W_0_0, W_p_0, source_0, mo_tmpr2, 1, tau, 2, 0, GLB_T);
		measure_mesons_core(W_0_0, W_0_p, source_0, mo_tmpr3, 1, tau, 2, 0, GLB_T);
		measure_mesons_core(W_0_mp, W_0_0, source_0, mo_tmpr4, 1, tau, 2, 0, GLB_T);

		measure_mesons_core(W_mp2_0, W_0_0, source_0, mo_tmpr1p2, 1, tau, 2, 0, GLB_T);
		measure_mesons_core(W_0_0, W_p2_0, source_0, mo_tmpr2p2, 1, tau, 2, 0, GLB_T);
		measure_mesons_core(W_0_0, W_0_p2, source_0, mo_tmpr3p2, 1, tau, 2, 0, GLB_T);
		measure_mesons_core(W_0_mp2, W_0_0, source_0, mo_tmpr4p2, 1, tau, 2, 0, GLB_T);

		//Pick only the relevant contributions
		for(int px_=0;px_<2;++px_)for(int py_=0;py_<2;++py_)for(int pz_=0;pz_<2;++pz_){
#define X(NAME) mo_##NAME->corr_re[corr_ind(px_,py_,pz_,2,t,1,0)] = mo_tmp##NAME->corr_re[corr_ind(px_,py_,pz_,2,t,1,0)]; mo_##NAME->corr_im[corr_ind(px_,py_,pz_,2,t,1,0)] = mo_tmp##NAME->corr_im[corr_ind(px_,py_,pz_,2,t,1,0)];
			TMP_OBSERVABLE_LIST
#undef X
		}
	}

#define X(NAME) do_global_sum(mo_##NAME, 1.0); io2pt(mo_##NAME,2,src,path, #NAME);
	TMP_OBSERVABLE_LIST
#undef X
  }

    if(list==NULL) break;
}
  lprintf("DEBUG",0,"ALL done, deallocating\n");

#define X(NAME) free_mo(mo_##NAME);
  OBSERVABLE_LIST
#undef X
#define X(NAME) free_mo(mo_tmp##NAME);
	  TMP_OBSERVABLE_LIST
#undef X
  if(list!=NULL) fclose(list);
  finalize_process();
  free_BCs();
  free_gfield(u_gauge);
  free_propagator_eo();

  return 0;
}
