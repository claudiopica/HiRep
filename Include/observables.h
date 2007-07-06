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
#include <stdio.h>

float plaq(int ix,int mu,int nu);
double avr_plaquette();
double local_plaq(int ix);

void quark_propagator(unsigned int source, int nm, float *mass, suNf_spinor **out);
/*void quark_propagator_QMR(unsigned int source, int nm, float *mass, suNf_spinor_dble **out);*/
void quark_propagator_QMR(FILE *propfile, unsigned int ssite, int nm, float *mass);

void id_correlator(float *out, suNf_spinor **qp);
void g0_correlator(float *out, suNf_spinor **qp);
void g5_correlator(float *out, suNf_spinor **qp);
void g0g5_correlator(float *out, suNf_spinor **qp);
void g1_correlator(float *out, suNf_spinor **qp);
void g2_correlator(float *out, suNf_spinor **qp);
void g3_correlator(float *out, suNf_spinor **qp);
void g0g1_correlator(float *out, suNf_spinor **qp);
void g0g2_correlator(float *out, suNf_spinor **qp);
void g0g3_correlator(float *out, suNf_spinor **qp);
void g5g1_correlator(float *out, suNf_spinor **qp);
void g5g2_correlator(float *out, suNf_spinor **qp);
void g5g3_correlator(float *out, suNf_spinor **qp);
void g0g5g1_correlator(float *out, suNf_spinor **qp);
void g0g5g2_correlator(float *out, suNf_spinor **qp);
void g0g5g3_correlator(float *out, suNf_spinor **qp);

#endif 
