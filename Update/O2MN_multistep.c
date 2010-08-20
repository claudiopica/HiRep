/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *   
* All rights reserved.                                                      * 
\***************************************************************************/

#include "global.h"
#include "suN.h"
#include "utils.h"
#include "update.h"
#include "representation.h"
#include "logger.h"
#include "communications.h"
#include <assert.h>

#define _PROJ_BIT (1<<4) /* project gauge field every 2^_PROJ_BIT changes */
#define _proj_gfield(c) if((c)&_PROJ_BIT) \
                          {(c)=0;project_gauge_field();} \
                        else \
                          {++(c); start_gf_sendrecv(u_gauge);} \
                        complete_gf_sendrecv(u_gauge); (void*)0

const double lambda = 0.1931833275037836; /* Omelyan et al */ 
                                         /* lamdba = 1./6. Sexton Weingarten */ 

static unsigned int count=0; /*used to project gfield */

static void O2MN_ms_gauge(suNg_av_field *momenta, double dt, unsigned int nsteps){
  _DECLARE_INT_ITERATOR(i);

  int n;

  assert(nsteps>0);

  lprintf("MD_INT_GAUGE",30,"GAUGE_INT starting. gauge_steps = %d\n",nsteps);

  /* Update of momenta */
  Force0(lambda*dt,momenta);

  for(n=1;n<nsteps;++n) {
    /* update of the gauge field */
    _MASTER_FOR(&glattice,i) {
      ExpX(dt*0.5,_4FIELD_AT(momenta,i,0), _4FIELD_AT(u_gauge,i,0));
      ExpX(dt*0.5,_4FIELD_AT(momenta,i,1), _4FIELD_AT(u_gauge,i,1));
      ExpX(dt*0.5,_4FIELD_AT(momenta,i,2), _4FIELD_AT(u_gauge,i,2));
      ExpX(dt*0.5,_4FIELD_AT(momenta,i,3), _4FIELD_AT(u_gauge,i,3));
    }
    _proj_gfield(count);
    //represent_gauge_field(); /* not needed in intermediate steps */
    
    /* Update of momenta */
    Force0((1.-2.*lambda)*dt,momenta);

    /* update of the gauge field */
    _MASTER_FOR(&glattice,i) {
      ExpX(dt*0.5,_4FIELD_AT(momenta,i,0), _4FIELD_AT(u_gauge,i,0));
      ExpX(dt*0.5,_4FIELD_AT(momenta,i,1), _4FIELD_AT(u_gauge,i,1));
      ExpX(dt*0.5,_4FIELD_AT(momenta,i,2), _4FIELD_AT(u_gauge,i,2));
      ExpX(dt*0.5,_4FIELD_AT(momenta,i,3), _4FIELD_AT(u_gauge,i,3));
    }
    _proj_gfield(count);
    //represent_gauge_field(); /* not needed in intermediate steps */

    /* Update of momenta */
    Force0(2.*lambda*dt,momenta);

    lprintf("MD_INT_GAUGE",30,"GAUGE_INT step: %d/%d\n",n,nsteps);
  }
   
  /* update of the gauge field */
  _MASTER_FOR(&glattice,i) {
    ExpX(dt*0.5,_4FIELD_AT(momenta,i,0), _4FIELD_AT(u_gauge,i,0));
    ExpX(dt*0.5,_4FIELD_AT(momenta,i,1), _4FIELD_AT(u_gauge,i,1));
    ExpX(dt*0.5,_4FIELD_AT(momenta,i,2), _4FIELD_AT(u_gauge,i,2));
    ExpX(dt*0.5,_4FIELD_AT(momenta,i,3), _4FIELD_AT(u_gauge,i,3));
  }
  _proj_gfield(count);
  //represent_gauge_field(); /* not needed in intermediate steps */
    
  /* Update of momenta */
  Force0((1.-2.*lambda)*dt,momenta);
    
  /* update of the gauge field */
  _MASTER_FOR(&glattice,i) {
    ExpX(dt*0.5,_4FIELD_AT(momenta,i,0), _4FIELD_AT(u_gauge,i,0));
    ExpX(dt*0.5,_4FIELD_AT(momenta,i,1), _4FIELD_AT(u_gauge,i,1));
    ExpX(dt*0.5,_4FIELD_AT(momenta,i,2), _4FIELD_AT(u_gauge,i,2));
    ExpX(dt*0.5,_4FIELD_AT(momenta,i,3), _4FIELD_AT(u_gauge,i,3));
  }
  _proj_gfield(count);
  represent_gauge_field();
    
  /* Update of momenta */
  Force0(lambda*dt,momenta);

  lprintf("MD_INT_GAUGE",30,"GAUGE_INT completed.\n");
}

void O2MN_multistep(suNg_av_field *momenta, int_par *traj_par){

  /* check input types */
#ifndef CHECK_SPINOR_MATCHING
   _TWO_SPINORS_MATCHING(u_gauge,momenta);
#endif

  count = 0; /* reset counter */

  if (traj_par->nsteps>0) {
    int n;
    double dt=traj_par->tlen/((double)traj_par->nsteps);
    double gdt=dt/(double)(2*traj_par->gsteps);
		  
    lprintf("MD_INT",10,"Starting new MD trajectory with O2MN_multistep.\n");
    lprintf("MD_INT",20,"MD parameters: len=%1.4f steps=%d gauge_steps=%d => dt=%1.4f gdt=%1.4f\n",
	    traj_par->tlen,traj_par->nsteps,traj_par->gsteps,dt,gdt);

    /* Update of momenta */
    Force_rhmc_f(lambda*dt,momenta);
    
    for(n=1;n<traj_par->nsteps;++n) {
      
      /* Update gfield */
      O2MN_ms_gauge(momenta, gdt, traj_par->gsteps);
      
      /* Update of momenta */
      Force_rhmc_f((1.-2.*lambda)*dt,momenta);
      
      /* Update gfield */
      O2MN_ms_gauge(momenta, gdt, traj_par->gsteps);
      
      /* Update of momenta */
      Force_rhmc_f(2.*lambda*dt,momenta);

      lprintf("MD_INT",10,"MD step: %d/%d\n",n,traj_par->nsteps);
      
    }
    
    /* Update gfield */
    O2MN_ms_gauge(momenta, gdt, traj_par->gsteps);
    
    /* Update of momenta */
    Force_rhmc_f((1.-2.*lambda)*dt,momenta);
    
    /* Update gfield */
    O2MN_ms_gauge(momenta, gdt, traj_par->gsteps);
    
    /* Update of momenta */
    Force_rhmc_f(lambda*dt,momenta);

    lprintf("MD_INT",10,"MD trajectory completed\n");
  }
}


#undef _PROJ_BIT
#undef _proj_gfield
