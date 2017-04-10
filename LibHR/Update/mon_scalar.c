/***************************************************************************\
 * Copyright (c) 2017                                                     *
 * All rights reserved.                                                   *
 \***************************************************************************/

#include "global.h"
#include "update.h"
#include "logger.h"
#include "memory.h"
#include "dirac.h"
#include "linear_algebra.h"
#include "inverters.h"
#include "observables.h"
#include <stdlib.h>

static double total_action(scalar_field* la){
	double tot_action=0.0;
	_MASTER_FOR(&glattice,i){
		tot_action += *_FIELD_AT(la,i);
	}

	return tot_action;	
}
void blank(const struct _monomial *m){
	/* empty */
}

const spinor_field* scalar_pseudofermion(const struct _monomial *m)
{
	return NULL;
}

void scalar_add_local_action(const struct _monomial *m, scalar_field *loc_action)
{
	mon_scalar_par *par = (mon_scalar_par*)(m->data.par);
	double Msq = par->mass;
	Msq = Msq*Msq + 8.0;
	double lambda = par->lambda;
	double b4 = total_action(loc_action);
	_MASTER_FOR(&glattice,ix)
	{
//		printf("b4: %f\n", *_FIELD_AT(loc_action,ix));
		suNg_vector* Sx = _FIELD_AT(u_scalar,ix); //pu_scalar(ix);
		//printf("Sx = %f\n", *((double*) Sx));
		double SUSup=0.0;
		for(int mu=0; mu<4; mu++){
			suNg_vector *Sup = _FIELD_AT(u_scalar,iup(ix,mu)); //pu_scalar(iup(ix,mu));
			suNg *U = _4FIELD_AT(u_gauge, ix, mu); //pu_gauge(ix,mu);
			suNg_vector UtimesSup;
			SUSup=0.0;
			_suNg_multiply(UtimesSup,*U,*Sup); 	
			_vector_prod_re_g(SUSup,*Sx,UtimesSup);
			*_FIELD_AT(loc_action,ix) -= 2.0*SUSup;
		}
		double Ssq;
                _vector_prod_re_g(Ssq,*Sx,*Sx);
		*_FIELD_AT(loc_action,ix) += Msq*Ssq + lambda*Ssq*Ssq;
/*		printf("after: %f\n", *_FIELD_AT(loc_action,ix));
		printf("after: %f %f %f\n", Msq*Ssq, lambda*Ssq*Ssq, 2.0*SUSup);*/
	}
//	printf("Scalar action: %f\n", total_action(loc_action) - b4);
}

void scalar_free(struct _monomial *m)
{
	mon_scalar_par *par = (mon_scalar_par*)m->data.par;
	free(par);
	free(m);
}

struct _monomial* scalar_create(const monomial_data *data)
{
	if(u_scalar==NULL){
  		u_scalar=alloc_scalar_field(&glattice);
	}
	monomial *m = malloc(sizeof(*m));
	mon_scalar_par *par = (mon_scalar_par*)(data->par);

	// Copy data structure
	m->data = *data;

	//Setup force parameters
	par->force_par.mass = par->mass;
	par->force_par.lambda = par->lambda;
	par->force_par.momenta = &scalar_momenta;
	par->force_par.g_momenta = &suN_momenta;
	par->field_par.field = &u_scalar;
	par->field_par.momenta = &scalar_momenta;

	// Setup pointers to update functions
	m->free = &scalar_free;
	m->update_force = &force_scalar;
	m->force_par = &par->force_par;
	m->update_field = &update_scalar_field;
	m->field_par = &par->field_par;	

	m->pseudofermion = &scalar_pseudofermion;
	m->gaussian_pf = &blank;
	m->correct_pf = &blank;
	m->correct_la_pf = &blank;
	m->add_local_action = &scalar_add_local_action;

	return m;
}

