/***************************************************************************\
* Copyright (c) 2011, Agostino Patella                                      *   
* All rights reserved.                                                      * 
\***************************************************************************/

/*******************************************************************************
*
* File wilsonflow.h
* 
*******************************************************************************/

#ifndef WILSONFLOW_H
#define WILSONFLOW_H

#include "suN_types.h"
#include "complex.h"
#include "spinor_field.h"
#include "suN.h"

void WF_initialize();
void WF_free();

double max_distance(suNg_field* V, suNg_field* Vprime);

void WilsonFlow1(suNg_field* V, const double epsilon);
void WilsonFlow3(suNg_field* V, const double epsilon);
double WilsonFlow3_adaptative(suNg_field* V, double epsilon,double delta);

double WF_E(suNg_field* V);
double WF_Esym(suNg_field* V);
double WF_topo(suNg_field* V);

void WF_E_T(double* E,suNg_field* V);
void WF_Esym_T(double* Esym,suNg_field* V);

#endif /* WILSONFLOW_H */

