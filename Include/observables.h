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

#define SPIN_2D_INDEX(i,j) ( (i)*4 + (j) )

double plaq(int ix,int mu,int nu);
double avr_plaquette();
double local_plaq(int ix);
void full_plaquette();

void pta_qprop_QMR_eo(spinor_field **pta_qprop, int nm, double *m, double acc);
void pta_qprop_QMR(spinor_field **pta_qprop, int nm, double *m, double acc);
void pta_qprop_MINRES(spinor_field **pta_qprop, int nm, double *m, double acc);
void traced_ata_qprop_dublin_trunc(complex** ev_prop, complex*** prop, int n_points);
void ata_qprop_dublin_trunc_init(int nm, double *mptr, int pnev, int pnevt, int ord, int trnc, int pgns, int pgnsh, int pgnsc, int pdil, int pdilh, int pdilc);
void ata_qprop_dublin_trunc_modify_par(int ord, int trnc, int pgns, int pgnsh, int pgnsc, int pdil, int pdilh, int pdilc);
void ata_qprop_dublin_trunc_free();


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

void id_trace_H(complex* out, complex* smat);
void g0_trace_H(complex* out, complex* smat);
void g5_trace_H(complex* out, complex* smat);
void g0g5_trace_H(complex* out, complex* smat);
void g1_trace_H(complex* out, complex* smat);
void g2_trace_H(complex* out, complex* smat);
void g3_trace_H(complex* out, complex* smat);
void g0g1_trace_H(complex* out, complex* smat);
void g0g2_trace_H(complex* out, complex* smat);
void g0g3_trace_H(complex* out, complex* smat);
void g5g1_trace_H(complex* out, complex* smat);
void g5g2_trace_H(complex* out, complex* smat);
void g5g3_trace_H(complex* out, complex* smat);
void g0g5g1_trace_H(complex* out, complex* smat);
void g0g5g2_trace_H(complex* out, complex* smat);
void g0g5g3_trace_H(complex* out, complex* smat);

void id_debug(complex Gamma[4][4], int* sign);
void g0_debug(complex Gamma[4][4], int* sign);
void g5_debug(complex Gamma[4][4], int* sign);
void g0g5_debug(complex Gamma[4][4], int* sign);
void g1_debug(complex Gamma[4][4], int* sign);
void g2_debug(complex Gamma[4][4], int* sign);
void g3_debug(complex Gamma[4][4], int* sign);
void g0g1_debug(complex Gamma[4][4], int* sign);
void g0g2_debug(complex Gamma[4][4], int* sign);
void g0g3_debug(complex Gamma[4][4], int* sign);
void g5g1_debug(complex Gamma[4][4], int* sign);
void g5g2_debug(complex Gamma[4][4], int* sign);
void g5g3_debug(complex Gamma[4][4], int* sign);
void g0g5g1_debug(complex Gamma[4][4], int* sign);
void g0g5g2_debug(complex Gamma[4][4], int* sign);
void g0g5g3_debug(complex Gamma[4][4], int* sign);

#endif 
