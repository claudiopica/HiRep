#include "global.h"
#include "suN.h"
#include "utils.h"
#include "update.h"
#include "representation.h"


/* extern suNg_algebra_vector *momenta; */


#define _PROJ_BIT (1<<6) /* project gauge field every 2^_PROJ_BIT changes */
#define _proj_leapfrog(c) if((c)&_PROJ_BIT){(c)=0;project_gauge_field();} else ++c


void leapfrog(suNg_algebra_vector *momenta, float tlen, unsigned int nsteps){
  static unsigned int count=0;

  if (nsteps>0) {
    int i, n;
    float dt=tlen/((float)nsteps);

    /* half step for the gauge field */
    for(i=0;i<4*VOLUME;++i){
      ExpX(dt/2.,momenta+i, u_gauge+i);
    }
    _proj_leapfrog(count);
    represent_gauge_field();

    for(n=1;n<nsteps;++n) {
      /* Update of momenta */
      Force(dt,momenta);

      /* update of the gauge field */
      for(i=0;i<4*VOLUME;++i){
				ExpX(dt,momenta+i, u_gauge+i);
      }
      _proj_leapfrog(count);
      represent_gauge_field();
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

}

#undef _PROJ_BIT
#undef _proj_leapfrog
