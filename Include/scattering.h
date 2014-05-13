
/*******************************************************************************
*
* File scattering.h
* 
* Functions for Finite Size Method 
*
*******************************************************************************/

#ifndef SCATTERING_H
#define SCATTERING_H

#include "suN.h"
#include <stdio.h>

void measure_pion_scattering(double* m, int nhits,int conf_num, double precision,int ts);
void contract_pion_scatt(spinor_field* phi_ts,spinor_field* phi_tsp1,int k,int ts);
void contract_pion_scatt_1spinorfield(spinor_field* phi_ts,spinor_field* phi_tsp1,int k,int ts);

#endif
