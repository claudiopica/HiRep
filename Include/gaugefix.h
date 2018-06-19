/***************************************************************************\
* Copyright (c) 2011, Agostino Patella                                      *   
* All rights reserved.                                                      * 
\***************************************************************************/

/*******************************************************************************
*
* File gaugefix.h
* 
*******************************************************************************/

#ifndef GAUGEFIX_H
#define GAUGEFIX_H

#include "suN_types.h"
#include "complex.h"
#include "spinor_field.h"
#include "suN.h"



void unit_gauge(suNg_field *original);

void random_gauge_transform(suNg_field *gauge);

void reunit(suNg_field *fixed_gauge);

double calc_plaq(suNg_field* V);

void su2_hit(int fix_dir, int parity, double overrelax, suNg_field *fixed_gauge, int c );

double gaugefix_action(int fix_dir, suNg_field *gauge );

double gaugefixstep(int fix_dir,double overrelax, suNg_field *fixed_gauge );

double scalar_gaugefixstep(int fix_dir,double overrelax, suNg_field *fixed_gauge, suNg_scalar_field *fixed_scalar);

double gaugefix(int fix_dir,double overrelax,int max_it,
	      double fix_tol, suNg_field *fixed_gauge );

double scalar_gaugefix(int fix_dir,double overrelax,int max_it, double fix_tol, suNg_field *fixed_gauge, suNg_scalar_field *fixed_scalar);

void hmc_Landau_gaugefix();
void hmc_Coulomb_gaugefix();
void free_hmc_gaugefix();

#endif /* GAUGEFIX_H */

