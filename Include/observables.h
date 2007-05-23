/*******************************************************************************
*
* File observables.h
* 
* Functions for measuring observables
*
*******************************************************************************/

#ifndef OBSERVABLES_H
#define OBSERVABLES_H

#include "suN.h"
#include "inverters.h"

float plaq(int ix,int mu,int nu);
double avr_plaquette();
double local_plaq(int ix);

void quark_propagator(unsigned int source, int nm, float *mass, suNf_spinor **out);
void quark_propagator_QMR(unsigned int source, int nm, float *mass, suNf_spinor_dble **out);
void pi_correlator(float *out, suNf_spinor *qp);
void pi_correlator_QMR(float *out, suNf_spinor_dble *qp);

#endif 
