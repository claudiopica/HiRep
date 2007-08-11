#include "global.h"
#include "suN.h"
#include "utils.h"
#include "update.h"
#include "representation.h"
#include "logger.h"

#define _PROJ_BIT (1<<6) /* project gauge field every 2^_PROJ_BIT changes */
#define _proj_leapfrog(c) if((c)&_PROJ_BIT){(c)=0;project_gauge_field();} else ++c

void leapfrog(suNg_algebra_vector *momenta, int_par *traj_par){
  static unsigned int count=0;

  if (traj_par->nsteps>0) {
    int i, n;
    double dt=traj_par->tlen/((double)traj_par->nsteps);

		lprintf("MD_INT",10,"Starting new MD trajectory using leapfrog.\n");
		lprintf("MD_INT",20,"MD parameters: len=%1.4f steps=%d => dt=%1.4f\n",
				traj_par->tlen,traj_par->nsteps,dt);

    /* half step for the gauge field */
    for(i=0;i<4*VOLUME;++i){
      ExpX(dt/2.,momenta+i, u_gauge+i);
    }
    _proj_leapfrog(count);
    represent_gauge_field();

    for(n=1;n<traj_par->nsteps;++n) {
      /* Update of momenta */
      Force(dt,momenta);

      /* update of the gauge field */
      for(i=0;i<4*VOLUME;++i){
				ExpX(dt,momenta+i, u_gauge+i);
      }
      _proj_leapfrog(count);
      represent_gauge_field();

			lprintf("MD_INT",10,"MD step: %d/%d\n",n,traj_par->nsteps);
    }
   
    /* Update of momenta */
    Force(dt,momenta);


    /* half step for the gauge field */
    for(i=0;i<4*VOLUME;++i){
      ExpX(dt/2.,momenta+i, u_gauge+i);
    }
    _proj_leapfrog(count);
    represent_gauge_field();
  }

	lprintf("MD_INT",10,"MD trajectory completed.\n");

}

#undef _PROJ_BIT
#undef _proj_leapfrog
