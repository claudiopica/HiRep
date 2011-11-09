/***************************************************************************\
 * Copyright (c) 2008, Claudio Pica                                          *   
 * All rights reserved.                                                      * 
 \***************************************************************************/

/*******************************************************************************
 *
 * Main RHMC program
 *
 *******************************************************************************/

#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "io.h"
#include "ranlux.h"
#include "geometry.h"
#include "update.h"
#include "global.h"
#include "observables.h"
#include "dirac.h"
#include "logger.h"
#include "rhmc_utils.h"
#include "memory.h"
#include "communications.h"
#include "observables.h"
#include "utils.h"

#include "cinfo.c"

#if defined(ROTATED_SF) && defined(UPDATE_EO)
#error ROTATED_SF DOES NOT WORK WITH E/O PRECONDITIONING
#endif


/* Mesons parameters */
typedef struct _input_mesons {
  int domes;
  double precision;
  int nhits;

  /* for the reading function */
  input_record_t read[4];

} input_mesons;

#define init_input_mesons(varname) \
{ \
  .domes=0,\
  .read={\
    {"do mesons", "mes:domes = %d", INT_T, &(varname).domes},\
    {"inverter precision", "mes:precision = %lf", DOUBLE_T, &(varname).precision},\
    {"number of inversions per cnfg", "mes:nhits = %d", INT_T, &(varname).nhits},\
    {NULL, NULL, 0, NULL}\
  }\
}

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

  free_spinor_field(pta_qprop[0]);
  free(pta_qprop);
  free(tricorr);

}

/* flow control variable */
rhmc_flow flow=init_rhmc_flow(flow);


void read_cmdline(int argc, char* argv[]) {
  int i, am=0;

  for (i=1;i<argc;i++) {
    if (strcmp(argv[i],"-m")==0) am=i;
  }

  if (am != 0) {
    print_compiling_info();
    exit(0);
  }

}

int main(int argc,char *argv[])
{
  int i, acc, rc;
  char sbuf[128];
  double mass;

  /* setup process id and communications */
  setup_process(&argc,&argv);

  read_cmdline(argc,argv);

  /* logger setup */
    logger_setlevel(0,10);
  /* disable logger for MPI processes != 0 */
  if (PID!=0) { logger_disable(); }

  if (PID==0) {
    sprintf(sbuf,">>out_%d",PID);  logger_stdout(sbuf); 
    /* sprintf(sbuf,"err_%d",PID); freopen(sbuf,"w",stderr);  */
  }

  lprintf("MAIN",0,"Compiled with macros: %s\n",MACROS); 
  lprintf("MAIN",0,"PId =  %d [world_size: %d]\n\n",PID,WORLD_SIZE); 

  /* read input file */
  read_input(glb_var.read,"input_file");
  read_input(mes_var.read,"input_file");
  
  if(glb_var.rlxd_state[0]!='\0')
  {
  	/*load saved state*/
	lprintf("MAIN",0,"Loading rlxd state from file %s\n",glb_var.rlxd_state);
	read_ranlxd_state(glb_var.rlxd_state);
  }
  else
  {
  lprintf("MAIN",0,"RLXD [%d,%d]\n",glb_var.rlxd_level,glb_var.rlxd_seed+PID);
  rlxd_init(glb_var.rlxd_level,glb_var.rlxd_seed+PID);
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
  lprintf("MAIN",0,"Fermion boundary conditions: %.2f,%.2f,%.2f,%.2f\n",bc[0],bc[1],bc[2],bc[3]);

  /* setup lattice geometry */
  geometry_mpi_eo();
  /* test_geometry_mpi_eo(); */
  lprintf("MAIN",0,"Geometry buffers: %d\n",glattice.nbuffers);

  /* Init Monte Carlo */
  init_mc(&flow, "input_file");
  lprintf("MAIN",0,"MVM during RHMC initialzation: %ld\n",getMVM());
  lprintf("MAIN",0,"Initial plaquette: %1.8e\n",avr_plaquette());
#ifdef BASIC_SF
  lprintf("MAIN",0,"Initial SF_action: %1.8e\n",SF_action((&flow)->rhmc_v->rhmc_p.beta));
#ifndef NDEBUG
  lprintf("MAIN",0,"Initial SF_test_gauge_bcs: %1.8e\n",SF_test_gauge_bcs());
#endif /*NDEBUG*/
#endif /* BASIC_SF */

  mass=flow.rhmc_v->rhmc_p.mass;

  rc=acc=0;
  for(i=flow.start;i<flow.end;++i) {
    int rr;
    double perc;
    struct timeval start, end, etime; /* //for trajectory timing */
    lprintf("MAIN",0,"Trajectory #%d...\n",i);
    
    gettimeofday(&start,0);
    
    rr=update_rhmc();

    gettimeofday(&end,0);
    timeval_subtract(&etime,&end,&start);
    lprintf("MAIN",0,"Trajectory #%d: generated in [%ld sec %ld usec]\n",i,etime.tv_sec,etime.tv_usec);
      
    if(rr<0) {
      lprintf("MAIN",0,"Error in updating the gauge field!!\n");
      return 1;
    } else if(rr!=0) {
      acc++;
    }
    rc++;
    perc=(acc==0)?0.:(float)(100*acc)/(float)(rc);

    lprintf("MAIN",0,"Trajectory #%d: %d/%d (%3.4f%%) MVM = %ld\n",i,acc,rc,perc,getMVM());

    if((i%flow.save_freq)==0) {
      save_conf(&flow, i);
    }

#ifdef BASIC_SF
    lprintf("MAIN",0,"SF action: %1.8e\n",SF_action((&flow)->rhmc_v->rhmc_p.beta));
#ifndef NDEBUG
      lprintf("MAIN",0,"SF_test_gauge_bcs: %1.8e\n",SF_test_gauge_bcs());
#endif /*NDEBUG*/
#endif /* BASIC_SF */

    if((i%flow.meas_freq)==0) {
      /* plaquette */
      lprintf("MAIN",0,"Plaquette: %1.8e\n",avr_plaquette());
#ifdef BASIC_SF

    gettimeofday(&start,0);
	    
    lprintf("MAIN",0,"PCAC mass: %1.8e\n",SF_PCAC_wall_mass((&flow)->rhmc_v->rhmc_p.mass));
	    
    gettimeofday(&end,0);
    timeval_subtract(&etime,&end,&start);
    lprintf("MAIN",0,"SF Propagators generated in [%ld sec %ld usec]\n",etime.tv_sec,etime.tv_usec);

#endif /* BASIC_SF */

      /* Mesons */
      if (mes_var.domes) {
        int nn;
        for (nn=0;nn<mes_var.nhits;++nn) {
          inline_mk_mesons(&mass,1,mes_var.precision);
        }
      }
    }
  }
  /* save final configuration */
  if(((--i)%flow.save_freq)!=0) {
    save_conf(&flow, i);
  }

 /* Only save state if we have a file to save to */
  if(glb_var.rlxd_state[0]!='\0') {
      lprintf("MAIN",0,"Saving rlxd state to file %s\n",glb_var.rlxd_state);
      write_ranlxd_state(glb_var.rlxd_state);
  }
  
  /* finalize Monte Carlo */
  end_mc();

  /* close communications */
  finalize_process();

  return 0;

}
