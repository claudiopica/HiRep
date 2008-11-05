/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *   
* All rights reserved.                                                      * 
\***************************************************************************/

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

double plaq(int ix,int mu,int nu);
double avr_plaquette();
double local_plaq(int ix);

void pta_qprop_QMR_eo(spinor_field **pta_qprop, int nm, double *m, double acc);
void pta_qprop_QMR(spinor_field **pta_qprop, int nm, double *m, double acc);
void pta_qprop_MINRES(spinor_field **pta_qprop, int nm, double *m, double acc);

void id_correlator(double *out, spinor_field *qp);
void g0_correlator(double *out, spinor_field *qp);
void g5_correlator(double *out, spinor_field *qp);
void g0g5_correlator(double *out, spinor_field *qp);
void g1_correlator(double *out, spinor_field *qp);
void g2_correlator(double *out, spinor_field *qp);
void g3_correlator(double *out, spinor_field *qp);
void g0g1_correlator(double *out, spinor_field *qp);
void g0g2_correlator(double *out, spinor_field *qp);
void g0g3_correlator(double *out, spinor_field *qp);
void g5g1_correlator(double *out, spinor_field *qp);
void g5g2_correlator(double *out, spinor_field *qp);
void g5g3_correlator(double *out, spinor_field *qp);
void g0g5g1_correlator(double *out, spinor_field *qp);
void g0g5g2_correlator(double *out, spinor_field *qp);
void g0g5g3_correlator(double *out, spinor_field *qp);
void g5_g0g5_re_correlator(double *out, spinor_field *qp);
void g5_g0g5_im_correlator(double *out, spinor_field *qp);

/* void dublin_meson_correlators(double** correlator[], char corr_name[][256], int n_corr, int n_masses, double *mass); */

#endif 
