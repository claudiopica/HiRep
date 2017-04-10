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
#include "dirac.h"
#include "representation.h"
#include "linear_algebra.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "logger.h"
#include "communications.h"

/* State quantities for HMC */
static suNg_field *u_gauge_old=NULL;
static suNg_scalar_field *u_scalar_old=NULL;
static scalar_field *ff_sigma_old=NULL;
static scalar_field *ff_pi_old=NULL;

static scalar_field *la=NULL; /* local action field for Metropolis test */

static ghmc_par update_par;
static int init=0;


void init_ghmc(ghmc_par *par){
  
  if (init) {
    /* already initialized */
    lprintf("GHMC",0,"WARNING: GHMC already initialized!\nWARNING: Ignoring call to init_ghmc.\n");
    return;
  }
	
  lprintf("GHMC",0,"Initializing...\n");
	  
  /* allocate space for the backup copy of gfield */
  if(u_gauge_old==NULL) u_gauge_old=alloc_gfield(&glattice);
  suNg_field_copy(u_gauge_old,u_gauge);
 
  /* allocate space for the backup copy of the scalar field */
  if(u_scalar!=NULL){
	  if(u_scalar_old==NULL) u_scalar_old=alloc_scalar_field(&glattice);
	  suNg_scalar_field_copy(u_scalar_old,u_scalar);
  }

  /* allocate space for backup copy of four fermion fields */
  if( four_fermion_active ) {
    if(ff_sigma_old==NULL) ff_sigma_old=alloc_sfield(1,&glattice);
    if(ff_pi_old==NULL) ff_pi_old=alloc_sfield(1,&glattice); 
    scalar_field_copy( ff_sigma_old, ff_sigma );
    scalar_field_copy( ff_pi_old, ff_pi );
  }
  

  /* allocate momenta */
  if(suN_momenta==NULL) suN_momenta = alloc_avfield(&glattice);
  if(u_scalar!=NULL){
	  if(scalar_momenta==NULL) scalar_momenta = alloc_scalar_field(&glattice);
  }
  
  /* allocate pseudofermions */
  /* we allocate one more pseudofermion for the computation
   * of the final action density
   */

  /* allocate memory for the local action */
  /* NOTE: should this be moved into local_action.c ? */
  if(la==NULL) la=alloc_sfield(1, &glattice);
  
  /* represent gauge field */
  represent_gauge_field();
  
  /* copy update parameters */
  update_par = *par;
  
//#ifdef ROTATED_SF
//  hmc_action_par.SF_ct = _update_par.SF_ct;
//#endif
  init = 1;
  
  lprintf("HMC",0,"Initialization done.\n");
}

void free_ghmc()
{
	if (!init)
	{
		/* not initialized */
		lprintf("HMC",0,"WARNING: HMC not initialized!\nWARNNG: Ignoring call to free_hmc.\n");
		return;
	}

	/* free momenta */
	if(u_gauge_old!=NULL) { free_gfield(u_gauge_old); u_gauge_old=NULL; }
	if(u_scalar_old!=NULL){ free_scalar_field(u_scalar_old); u_scalar_old=NULL;}
	if(suN_momenta!=NULL) { free_avfield(suN_momenta); suN_momenta=NULL; }
	if(scalar_momenta!=NULL){ free_scalar_field(scalar_momenta); scalar_momenta=NULL;}
	if(la!=NULL) { free_sfield(la); la=NULL; }
  
	/*Free integrator */
	integrator_par *ip = update_par.integrator;
	while(ip != NULL)
	{
		update_par.integrator=ip->next;
		free(ip->mon_list);
		free(ip);
		ip = update_par.integrator;
	}
	update_par.integrator = NULL;

	//free_force_hmc();
	init = 0;
	lprintf("HMC",0,"Memory deallocated.\n");
}


int update_ghmc()
{
  double deltaH;

  if(!init)
    {
      /* not initialized */
      lprintf("HMC",0,"WARNING: GHMC not initialized!\nWARNNG: Ignoring call to update_ghmc.\n");
      return -1;
    }

  /* generate new momenta */
  lprintf("HMC",30,"Generating gaussian momenta and pseudofermions...\n");
  gaussian_momenta(suN_momenta);
  if(u_scalar!=NULL){
	  gaussian_scalar_momenta(scalar_momenta);
  }

  /* generate new pseudofermions */
  for (int i=0;i<num_mon();++i) {
    const monomial *m = mon_n(i);
    m->gaussian_pf(m);
  }

  /* compute starting action */
  lprintf("HMC",30,"Computing action density...\n");
  local_hmc_action(NEW, la, suN_momenta, scalar_momenta);

  /* correct pseudofermion distribution */
  for (int i=0;i<num_mon();++i) {
    const monomial *m = mon_n(i);
    m->correct_pf(m);
  }
  

  /* integrate molecular dynamics */
  lprintf("HMC",30,"MD integration...\n");
  update_par.integrator->integrator(update_par.tlen,update_par.integrator);


  /* project and represent gauge field */
  project_gauge_field();
  represent_gauge_field();

  /* compute new action */
  lprintf("HMC",30,"Computing new action density...\n");
  for (int i=0;i<num_mon();++i) {
    const monomial *m = mon_n(i);
    m->correct_la_pf(m);
  }
  local_hmc_action(DELTA, la, suN_momenta, scalar_momenta);

  /* Metropolis test */
  deltaH = 0.0;
  _MASTER_FOR_SUM(la->type,i,deltaH) {
    deltaH += *_FIELD_AT(la,i);
  }

  global_sum(&deltaH, 1);
  lprintf("HMC",10,"[DeltaS = %1.8e][exp(-DS) = %1.8e]\n",deltaH,exp(-deltaH));

  if(deltaH < 0) {
    suNg_field_copy(u_gauge_old,u_gauge);
    if(u_scalar!=NULL){
	    suNg_scalar_field_copy(u_scalar_old,u_scalar);
    }
    if( four_fermion_active ) {
      scalar_field_copy( ff_sigma_old, ff_sigma );
      scalar_field_copy( ff_pi_old, ff_pi );
    }
  } else {
    double r;
    if(PID == 0) {
      ranlxd(&r,1);
      if(r < exp(-deltaH)) {
        r = 1.0;
      } else {
        r = -1.0;
      }
    }

    bcast(&r, 1);

    if(r > 0) {
      suNg_field_copy(u_gauge_old,u_gauge);
      if(u_scalar!=NULL){
	      suNg_scalar_field_copy(u_scalar_old,u_scalar);
      }
      if( four_fermion_active ) {
        scalar_field_copy( ff_sigma_old, ff_sigma );
        scalar_field_copy( ff_pi_old, ff_pi );
      }
    } else {
      lprintf("HMC",10,"Configuration rejected.\n");
      suNg_field_copy(u_gauge,u_gauge_old);
      if(u_scalar!=NULL){
	      suNg_scalar_field_copy(u_scalar,u_scalar_old);
      }
      if( four_fermion_active ) {
        scalar_field_copy( ff_sigma, ff_sigma_old );
        scalar_field_copy( ff_pi, ff_pi_old );
      }
      start_gf_sendrecv(u_gauge); /* this may not be needed if we always guarantee that we copy also the buffers */
      if(u_scalar!=NULL){
	      start_sc_sendrecv(u_scalar); /* this may not be needed if we always guarantee that we copy also the buffers */
      }
      represent_gauge_field();
      return 0;
    }
  }

  lprintf("HMC",10,"Configuration accepted.\n");
  return 1;
}


#ifdef MEASURE_FORCEHMC
/*Functions to check forces */
void corret_pf_dist_hmc(){
    /* init monomials */
  for (int i=0;i<num_mon();++i) {
    const monomial *m = mon_n(i);
    m->init_traj(m);
  }
  
  /* generate new momenta */
  lprintf("HMC",30,"Generating gaussian momenta and pseudofermions...\n");
  gaussian_momenta(suN_momenta);

  /* generate new pseudofermions */
  for (int i=0;i<num_mon();++i) {
    const monomial *m = mon_n(i);
    m->gaussian_pf(m);
  }

  /* compute starting action */
  lprintf("HMC",30,"Computing action density...\n");
  local_hmc_action(NEW, la, suN_momenta, scalar_momenta);

  /* correct pseudofermion distribution */
  for (int i=0;i<num_mon();++i) {
    const monomial *m = mon_n(i);
    m->correct_pf(m);
  }
}

void calc_one_force(int n_force){
  integrator_par *ip=update_par.integrator;
  for (;;){
    error(ip==NULL,1,"calc_one_force","Error in force index\n");
    for(int n = 0; n < ip->nmon; n++){
      const monomial *m=ip->mon_list[n];
      if (m->data.id==n_force){
        m->force_f(1,m->force_par);
        return;
      }
    }
    ip=ip->next;
  }
}
#endif
