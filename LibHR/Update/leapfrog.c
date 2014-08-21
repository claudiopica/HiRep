/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *   
* All rights reserved.                                                      * 
\***************************************************************************/
#if (1==0) //Old code. leapfrog Moved to file O2MN_multistep.c

#include "global.h"
#include "suN.h"
#include "utils.h"
#include "update.h"
#include "representation.h"
#include "logger.h"
#include "communications.h"

#define _PROJ_BIT (1<<4) /* project gauge field every 2^_PROJ_BIT changes */
#define _proj_gfield(c) if((c)&_PROJ_BIT){(c)=0;project_gauge_field();} else {++(c); start_gf_sendrecv(u_gauge);}

void leapfrog(suNg_av_field *momenta, int_par *traj_par){
  unsigned int count=0;
  _DECLARE_INT_ITERATOR(i);

  /* check input types */
#ifndef CHECK_SPINOR_MATCHING
  _TWO_SPINORS_MATCHING(u_gauge,momenta);
#endif

  if (traj_par->nsteps>0) {
    int n;
    double dt=traj_par->tlen/((double)traj_par->nsteps);

    lprintf("MD_INT",10,"Starting new MD trajectory using leapfrog.\n");
    lprintf("MD_INT",20,"MD parameters: len=%1.4f steps=%d => dt=%1.4f\n",
        traj_par->tlen,traj_par->nsteps,dt);

    /* half step for the gauge field */
    _MASTER_FOR(&glattice,i) {
      ExpX(dt*0.5,_4FIELD_AT(momenta,i,0), _4FIELD_AT(u_gauge,i,0));
      ExpX(dt*0.5,_4FIELD_AT(momenta,i,1), _4FIELD_AT(u_gauge,i,1));
      ExpX(dt*0.5,_4FIELD_AT(momenta,i,2), _4FIELD_AT(u_gauge,i,2));
      ExpX(dt*0.5,_4FIELD_AT(momenta,i,3), _4FIELD_AT(u_gauge,i,3));
    }
    _proj_gfield(count);
    represent_gauge_field();

    for(n=1;n<traj_par->nsteps;++n) {
      /* Update of momenta */
      Force(dt,momenta);

      /* update of the gauge field */
      _MASTER_FOR(&glattice,i) {
        ExpX(dt,_4FIELD_AT(momenta,i,0), _4FIELD_AT(u_gauge,i,0));
        ExpX(dt,_4FIELD_AT(momenta,i,1), _4FIELD_AT(u_gauge,i,1));
        ExpX(dt,_4FIELD_AT(momenta,i,2), _4FIELD_AT(u_gauge,i,2));
        ExpX(dt,_4FIELD_AT(momenta,i,3), _4FIELD_AT(u_gauge,i,3));
      }
      _proj_gfield(count);
      represent_gauge_field();

      lprintf("MD_INT",10,"MD step: %d/%d\n",n,traj_par->nsteps);
    }

    /* Update of momenta */
    Force(dt,momenta);

    /* half step for the gauge field */
    _MASTER_FOR(&glattice,i) {
      ExpX(dt*0.5,_4FIELD_AT(momenta,i,0), _4FIELD_AT(u_gauge,i,0));
      ExpX(dt*0.5,_4FIELD_AT(momenta,i,1), _4FIELD_AT(u_gauge,i,1));
      ExpX(dt*0.5,_4FIELD_AT(momenta,i,2), _4FIELD_AT(u_gauge,i,2));
      ExpX(dt*0.5,_4FIELD_AT(momenta,i,3), _4FIELD_AT(u_gauge,i,3));
    }
    _proj_gfield(count);
    represent_gauge_field();
  }

  lprintf("MD_INT",10,"MD trajectory completed.\n");

}

#undef _PROJ_BIT
#undef _proj_leapfrog
#endif
