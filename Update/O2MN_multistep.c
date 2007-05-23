#include "global.h"
#include "suN.h"
#include "utils.h"
#include "update.h"
#include "representation.h"

#include <assert.h>

#define _PROJ_BIT (1<<6) /* project gauge field every 2^_PROJ_BIT changes */
#define _proj_gfield(c) if((c)&_PROJ_BIT){(c)=0;project_gauge_field();} else ++(c)

const float lambda = 0.1931833275037836; /* Omelyan et al */ 
                                         /* lamdba = 1./6. Sexton Weingarten */ 

static void O2MN_ms_gauge(suNg_algebra_vector *momenta, float dt, unsigned int nsteps){
  static unsigned int count=0; /*used to project gfield */

    int i, n;

    assert(nsteps>0);

    /* Update of momenta */
    Force0(lambda*dt,momenta);

    for(n=1;n<nsteps;++n) {
      /* update of the gauge field */
      for(i=0;i<4*VOLUME;++i){
				ExpX(0.5*dt,momenta+i, u_gauge+i);
      }
      _proj_gfield(count);
      represent_gauge_field();

      /* Update of momenta */
      Force0((1.-2.*lambda)*dt,momenta);

      /* update of the gauge field */
      for(i=0;i<4*VOLUME;++i){
				ExpX(0.5*dt,momenta+i, u_gauge+i);
      }
      _proj_gfield(count);
      represent_gauge_field();

      /* Update of momenta */
      Force0(2.*lambda*dt,momenta);
    }
   
    /* update of the gauge field */
    for(i=0;i<4*VOLUME;++i){
      ExpX(0.5*dt,momenta+i, u_gauge+i);
    }
    _proj_gfield(count);
    represent_gauge_field();
    
    /* Update of momenta */
    Force0((1.-2.*lambda)*dt,momenta);
    
    /* update of the gauge field */
    for(i=0;i<4*VOLUME;++i){
      ExpX(0.5*dt,momenta+i, u_gauge+i);
    }
    _proj_gfield(count);
    represent_gauge_field();
    
    /* Update of momenta */
    Force0(lambda*dt,momenta);
}

void O2MN_multistep(suNg_algebra_vector *momenta, float tlen, unsigned int nsteps, unsigned int gsteps){

  if (nsteps>0) {
    int n;
    float dt=tlen/((float)nsteps);
    float gdt=dt/(double)(2*gsteps);
		  
    /* Update of momenta */
    Force_rhmc_f(lambda*dt,momenta);
    
    for(n=1;n<nsteps;++n) {
      
      /* Update gfield */
      O2MN_ms_gauge(momenta, gdt, gsteps);
      
      /* Update of momenta */
      Force_rhmc_f((1.-2.*lambda)*dt,momenta);
      
      /* Update gfield */
      O2MN_ms_gauge(momenta, gdt, gsteps);
      
      /* Update of momenta */
      Force_rhmc_f(2.*lambda*dt,momenta);
      
    }
    
    /* Update gfield */
    O2MN_ms_gauge(momenta, gdt, gsteps);
    
    /* Update of momenta */
    Force_rhmc_f((1.-2.*lambda)*dt,momenta);
    
    /* Update gfield */
    O2MN_ms_gauge(momenta, gdt, gsteps);
    
    /* Update of momenta */
    Force_rhmc_f(lambda*dt,momenta);
  }
}


#undef _PROJ_BIT
#undef _proj_gfield
