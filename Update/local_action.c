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
                      action_par *par,
                      scalar_field *loc_action,
                      suNg_av_field *momenta,
                      spinor_field *phi1,
                      spinor_field *phi2) {

  _DECLARE_INT_ITERATOR(i);
  int j;
  double a,tmp;

  /* check input types */
#ifndef CHECK_SPINOR_MATCHING
  _TWO_SPINORS_MATCHING(u_gauge,loc_action); /* check that action is defined on the global lattice */
  _TWO_SPINORS_MATCHING(loc_action,momenta);
  _TWO_SPINORS_MATCHING(loc_action,phi1);
  _TWO_SPINORS_MATCHING(loc_action,phi2);
#endif


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
    a=0.;
    /* Momenta */
    for (j=0;j<4;++j) {
      suNg_algebra_vector *cmom=momenta->ptr+coord_to_index(i,j);
      _algebra_vector_sqnorm_g(tmp,*cmom); 
      a+=tmp; /* this must be positive */
    }
    a*=0.5*_FUND_NORM2;
    *_FIELD_AT(loc_action,i)+=a;
  }

#ifndef ROTATED_SF
  _MASTER_FOR(&glattice,i) {

    /* Gauge action */
    *_FIELD_AT(loc_action,i) += -(par->beta/((double)NG))*local_plaq(i);
  }

#else /* ROTATED_SF */

  
  if(COORD[0]==0) {
    for(int x=0;x<X;x++) for(int y=0;y<Y;y++) for(int z=0;z<Z;z++) {
      i=ipt(1,x,y,z);
      *_FIELD_AT(loc_action,i) += -(.5*par->SF_ct*par->beta/((double)NG))*plaq(i,1,0);
      *_FIELD_AT(loc_action,i) += -(.5*par->SF_ct*par->beta/((double)NG))*plaq(i,2,0);
      *_FIELD_AT(loc_action,i) += -(.5*par->SF_ct*par->beta/((double)NG))*plaq(i,3,0);
    }
  } else{
    for(int t=0;t<2;t++) for(int x=0;x<X;x++) for(int y=0;y<Y;y++) for(int z=0;z<Z;z++) {
      i=ipt(t,x,y,z);
      *_FIELD_AT(loc_action,i) += -(par->beta/((double)NG))*local_plaq(i);
    }
  }
  
  for(int t=2;t<T-2;t++) for(int x=0;x<X;x++) for(int y=0;y<Y;y++) for(int z=0;z<Z;z++) {
    i=ipt(t,x,y,z);
    *_FIELD_AT(loc_action,i) += -(par->beta/((double)NG))*local_plaq(i);
  }
  
  if(COORD[0]==NP_T-1) {
    for(int x=0;x<X;x++) for(int y=0;y<Y;y++) for(int z=0;z<Z;z++) {
      i=ipt(T-2,x,y,z);
      *_FIELD_AT(loc_action,i) += -(.5*par->SF_ct*par->beta/((double)NG))*plaq(i,1,0);
      *_FIELD_AT(loc_action,i) += -(.5*par->SF_ct*par->beta/((double)NG))*plaq(i,2,0);
      *_FIELD_AT(loc_action,i) += -(.5*par->SF_ct*par->beta/((double)NG))*plaq(i,3,0);
      *_FIELD_AT(loc_action,i) += -(par->beta/((double)NG))*plaq(i,1,2);
      *_FIELD_AT(loc_action,i) += -(par->beta/((double)NG))*plaq(i,1,3);
      *_FIELD_AT(loc_action,i) += -(par->beta/((double)NG))*plaq(i,2,3);
    }
  } else {
    for(int t=T-2;t<T;t++) for(int x=0;x<X;x++) for(int y=0;y<Y;y++) for(int z=0;z<Z;z++) {
      i=ipt(t,x,y,z);
      *_FIELD_AT(loc_action,i) += -(par->beta/((double)NG))*local_plaq(i);
    }
  }

#endif /* ROTATED_SF */

  /* pseudofermion fields can be defined only on even sites is the preconditioning is used */
  _MASTER_FOR(phi1->type,i) {
    a=0.;
  /* Fermions */
    for (j=0;j<par->n_pf;++j) {
      _spinor_prod_re_f(tmp,*_FIELD_AT(&phi1[j],i),*_FIELD_AT(&phi2[j],i));
      a += tmp;
    }

    *_FIELD_AT(loc_action,i)+=a;
  }   
}
