/***************************************************************************\
 * Copyright (c) 2008, Agostino Patella, Claudio Pica                        *   
 * All rights reserved.                                                      * 
 \***************************************************************************/

#include "global.h"
#include "suN.h"
#include "utils.h"
#include "update.h"
#include "memory.h"
#include "random.h"
#include "representation.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "logger.h"
#include "communications.h"

/* State quantities for RHMC */
suNg_av_field *momenta=NULL;
puregauge_par _update_par={0};
integrator_par *integrator = NULL;
action_par puregauge_action_par;
/* END of State */


static short int init=0;


static suNg_field *u_gauge_old=NULL;
static scalar_field *la=NULL; /* local action field for Metropolis test */


void init_puregauge(puregauge_par *par){
  
	if (init) {
	  /* already initialized */
	  lprintf("GAUGE",0,"WARNING: GAUGE already initialized!\nWARNNG: Ignoring call to init_puregauge.\n");
		return;
	}
	
  lprintf("GAUGE",0,"Initializing...\n");
  
  /* copy input parameters into the internal variable and make some tests */
  _update_par=*par;
  /* no tests yet... */
  
  lprintf("GAUGE",10,
	  "beta = %.8f\n"
	  "Metropolis test precision = %.8e\n"
	  "MD precision = %.8e\n"
	  "MD trajectory length = %.8f\n"
	  "MD steps = %d\n"
	  ,_update_par.beta
	  ,_update_par.MT_prec
	  ,_update_par.MD_prec
	  ,_update_par.tlen
	  ,_update_par.gsteps
	  );
	
	/* allocate space for the backup copy of gfield */
	if(u_gauge_old==NULL) u_gauge_old=alloc_gfield(&glattice);
	suNg_field_copy(u_gauge_old,u_gauge);
	
	/* allocate momenta */
	if(momenta==NULL) momenta = alloc_avfield(&glattice);
	
	/* allocate memory for the local action */
	if(la==NULL) la=alloc_sfield(1, &glattice);

  /* integrator */
  integrator = (integrator_par*)malloc(sizeof(integrator_par)*2);

  integrator[0].level = 0;
  integrator[0].nsteps = _update_par.gsteps;
  integrator[0].force = &force0;
  integrator[0].force_par = (void*)malloc(sizeof(double));
  *((double*)integrator[0].force_par) = _update_par.beta;
  integrator[0].integrator = &O2MN_multistep;
  integrator[0].next = &integrator[1];

  integrator[1].level = 1;
  integrator[1].nsteps = 1;
  integrator[1].force = NULL;
  integrator[1].force_par = NULL;
  integrator[1].integrator = &gauge_integrator;
  integrator[1].next = NULL;
	
  puregauge_action_par.beta = _update_par.beta;
  puregauge_action_par.n_pf = 0;
	init = 1;
	
	lprintf("GAUGE",0,"Initialization done.\n");
	
}

void free_puregauge(){
  
  if (!init) {
    /* not initialized */
    lprintf("GAUGE",0,"WARNING: GAUGE not initialized!\nWARNNG: Ignoring call to free_puregauge.\n");
		return;
	}
	
	/* free momenta */
	if(u_gauge_old!=NULL) free_gfield(u_gauge_old); u_gauge_old=NULL;
	if(momenta!=NULL) free_avfield(momenta); momenta=NULL;
	if(la!=NULL) free_sfield(la); la=NULL;
  
  int i=0;
  while(1) {
    if(integrator[i].force_par != NULL) free(integrator[i].force_par);
    if(integrator[i].next == NULL) break;
    i++;
  }
  free(integrator);
  
	init = 0;
	
	lprintf("GAUGE",0,"Memory deallocated.\n");
	
}

int update_puregauge(){
    
    #ifdef TIMING
    struct timeval start1, end1;
    struct timeval start2, end2;
    struct timeval etime;
    
    #ifdef TIMING_WITH_BARRIERS
    MPI_Barrier(GLB_COMM);
    #endif
    gettimeofday(&start1,0);
    #endif
    
    double deltaH;
    _DECLARE_INT_ITERATOR(ix);
    
    if (!init) {
        /* not initialized */
        lprintf("GAUGE",0,"WARNING: GAUGE not initialized!\nWARNNG: Ignoring call to update_puregauge.\n");
		return -1;
    }
    
    /* generate new momenta and pseudofermions */
    lprintf("GAUGE",30,"Generating gaussian momenta...\n");
    gaussian_momenta(momenta);

    /* compute starting action */
    lprintf("GAUGE",30,"Computing action density...\n");
    local_hmc_action(NEW, &puregauge_action_par, la, momenta, NULL, NULL);
    
    /* integrate molecular dynamics */
    lprintf("GAUGE",30,"MD integration...\n");
    #ifdef TIMING
    #ifdef TIMING_WITH_BARRIERS
    MPI_Barrier(GLB_COMM);
    #endif
    gettimeofday(&start2,0);
    #endif
    (*(integrator[0].integrator))(momenta,_update_par.tlen,&integrator[0]);
    #ifdef TIMING
    #ifdef TIMING_WITH_BARRIERS
    MPI_Barrier(GLB_COMM);
    #endif
    gettimeofday(&end2,0);
    timeval_subtract(&etime,&end2,&start2);
    lprintf("TIMING",0,"integrator in update_puregauge %.6f s\n",1.*etime.tv_sec+1.e-6*etime.tv_usec);
    #endif
    
    /* project gauge field */
    project_gauge_field();
        
    lprintf("GAUGE",30,"Computing new action density...\n");
    
    /* compute new action */
    local_hmc_action(DELTA, &puregauge_action_par, la, momenta, NULL, NULL);
    
    /* Metropolis test */
    deltaH=0.;
    _MASTER_FOR(la->type,ix) {
        deltaH+=*_FIELD_AT(la,ix);
    }
    global_sum(&deltaH, 1);
    lprintf("GAUGE",10,"[DeltaS = %1.8e][exp(-DS) = %1.8e]\n",deltaH,exp(-deltaH));
    
    if(deltaH<0.) {
        suNg_field_copy(u_gauge_old,u_gauge);
    } else {
        double r;
        if (PID==0) { ranlxd(&r,1); if(r<exp(-deltaH)) r=1.0; else r=-1.0;}  /* make test on on PID 0 */
        bcast(&r,1);
        if(r>0.) {
            suNg_field_copy(u_gauge_old,u_gauge);
        } else {
            lprintf("GAUGE",10,"Configuration rejected.\n");
            suNg_field_copy(u_gauge,u_gauge_old);
            start_gf_sendrecv(u_gauge); /* this may not be needed if we always guarantee that we copy also the buffers */
            complete_gf_sendrecv(u_gauge);
            
            return 0;
        }
    }
    
    lprintf("GAUGE",10,"Configuration accepted.\n");

    #ifdef TIMING
    #ifdef TIMING_WITH_BARRIERS
    MPI_Barrier(GLB_COMM);
    #endif
    gettimeofday(&end1,0);
    timeval_subtract(&etime,&end1,&start1);
    lprintf("TIMING",0,"update_puregauge %.6f s\n",1.*etime.tv_sec+1.e-6*etime.tv_usec);
    #endif
    
    return 1;
}



