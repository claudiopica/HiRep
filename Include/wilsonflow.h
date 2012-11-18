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

void WilsonFlow1(suNg_field* V, const double epsilon);
void WilsonFlow3(suNg_field* V, const double epsilon);

double WF_E(suNg_field* V);
double WF_Esym(suNg_field* V);

#endif /* WILSONFLOW_H */

