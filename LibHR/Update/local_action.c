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
void local_hmc_action(local_action_type_old type,
                      scalar_field *loc_action,
                      suNg_av_field *momenta
                      ) {

  
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

  _MASTER_FOR(&glattice,i) {
    double a=0., tmp;
    /* Momenta */
    for (int j=0;j<4;++j) {
      suNg_algebra_vector *cmom=momenta->ptr+coord_to_index(i,j);
      _algebra_vector_sqnorm_g(tmp,*cmom); 
      a+=tmp; /* this must be positive */
    }
    a*=0.5*_FUND_NORM2;
    *_FIELD_AT(loc_action,i)+=a;
  }

  int nmon=num_mon();
  for (int i=0;i<nmon;++i) {
    const monomial *m=mon_n(i);
    if ( m->type == PureGauge ) {
        double beta=((mon_pg_par*)m->par)->beta;
        _MASTER_FOR(&glattice,i) {
          /* Gauge action */
          *_FIELD_AT(loc_action,i) += -(beta/((double)NG))*local_plaq(i);
        }
    } else {
      spinor_field *phi=NULL;
      switch (m->type) {
        case HMC:
          phi=((mon_hmc_par*)m->par)->pf;
          break;
        case RHMC:
          phi=((mon_rhmc_par*)m->par)->pf;
          break;
        case Hasenbusch:
          phi=((mon_hasenbusch_par*)m->par)->pf;
          break;
        default:
          lprintf("MONOMIAL",0,"WARNING: unknown type!\n");
          break;
      }
      /* pseudo fermion action = phi^2 */
      _MASTER_FOR(phi->type,i) {
        double a=0.;
        /* Fermions */
        _spinor_prod_re_f(a,*_FIELD_AT(phi,i),*_FIELD_AT(phi,i));
        *_FIELD_AT(loc_action,i)+=a;
      }

    }
  }
}
