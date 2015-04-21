
/*******************************************************************************
*
* File observables.h
* 
* Functions for measuring observables
*
*******************************************************************************/

#ifndef DISCONNECTED_H
#define DISCONNECTED_H

#include "suN.h"
#include <stdio.h>


void measure_bilinear_loops_spinorfield(spinor_field* prop,spinor_field* source,int k,int nm,int n_mom);
void measure_bilinear_loops_4spinorfield(spinor_field* prop,spinor_field* source,int k,int nm,int tau,int col,int eo);

void measure_loops(int nm, double* m, int nhits,int conf_num, double precision,int source_type,int n_mom);
void zero_even_or_odd_site_spinorfield(spinor_field *source,int nspinor,int eo);
#endif
