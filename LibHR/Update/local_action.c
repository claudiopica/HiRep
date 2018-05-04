/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *   
* All rights reserved.                                                      * 
\***************************************************************************/

#include "global.h"
#include "observables.h"
#include "update.h"
#include "suN.h"
#include "linear_algebra.h"
#include "error.h"
#include "representation.h"
#include "logger.h"

/*
 * compute the local action at every site for the HMC
 * H = | momenta |^2 + S_g + < phi1, phi2>
 */
void local_hmc_action(local_action_type type,
                      scalar_field *loc_action,
                      suNg_av_field *momenta,
                      suNg_scalar_field *momenta_s
                      ){

  
  /* check input types */
  _TWO_SPINORS_MATCHING(u_gauge,loc_action); /* check that action is defined on the global lattice */
  _TWO_SPINORS_MATCHING(loc_action,momenta);


  switch(type) {
  case NEW:
    _MASTER_FOR(&glattice,i) {
      *_FIELD_AT(loc_action,i)=0.;
    }
    break;
  case DELTA:
    _MASTER_FOR(&glattice,i) {
      *_FIELD_AT(loc_action,i)= -*_FIELD_AT(loc_action,i);
    }
    break;
  default:
    error(1,1,"local_hmc_action","Invalid type");
  }

#ifdef DEBUGACTION
   double total_action = 0.0;
    _MASTER_FOR_SUM(&glattice,j,total_action){
	total_action += *_FIELD_AT(loc_action,j);
    }
    global_sum(&total_action,1);
    lprintf("DEBUG ACTION",0,"Action b4 momenta: %f\n",total_action);
#endif

  _MASTER_FOR(&glattice,i) {
    double a=0., tmp;
    /* Momenta */
    for (int j=0;j<4;++j) {
      suNg_algebra_vector *cmom=momenta->ptr+coord_to_index(i,j);
      _algebra_vector_sqnorm_g(tmp,*cmom); 
      a+=tmp; /* this must be positive */
    }
    double P2=0.0; //Scalar momentum squared
    if(u_scalar!=NULL){
	    suNg_vector P=*_FIELD_AT(momenta_s,i);
	    _vector_prod_re_g(P2,P,P);
    }
    a*=0.5*_FUND_NORM2;
    *_FIELD_AT(loc_action,i)+=a + P2;
  }

#ifdef DEBUGACTION
    total_action = 0.0;	
    _MASTER_FOR_SUM(&glattice,j,total_action){
	total_action += *_FIELD_AT(loc_action,j);
    }
    global_sum(&total_action,1);
    lprintf("DEBUG ACTION",0,"Action after momenta: %f\n",total_action);
#endif

  int nmon=num_mon();
  for (int i=0;i<nmon;++i) {
    const monomial *m=mon_n(i);
    
    m->add_local_action(m,loc_action);

#ifdef DEBUGACTION
    total_action = 0.0;	
    _MASTER_FOR_SUM(&glattice,j,total_action){
	total_action += *_FIELD_AT(loc_action,j);
    }
    global_sum(&total_action,1);
    lprintf("DEBUG ACTION",0,"Monomial type %d, monomial number %d, action after: %f\n",m->data.type, i, total_action);
#endif

  }
}

/* add the square of the pf field to the local action */
void pf_local_action(scalar_field *loc_action, spinor_field *pf) {
  if (pf!=NULL) {
    _MASTER_FOR(pf->type,i) {
      double a=0.;
      /* Fermions */
      _spinor_prod_re_f(a,*_FIELD_AT(pf,i),*_FIELD_AT(pf,i));
      *_FIELD_AT(loc_action,i)+=a;
    }
  }
}



