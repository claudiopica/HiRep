/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *   
* All rights reserved.                                                      * 
\***************************************************************************/

/*******************************************************************************
*
* File utils.h
* 
* Some useful functions
*
*******************************************************************************/

#ifndef UTILS_H
#define UTILS_H

#include "suN_types.h"
#include "spinor_field.h"

void ExpX(double dt, suNg_algebra_vector *h, suNg *u);
void apply_bc();
void apply_bc_flt();

void cross_prod(suNg_vector *v1,suNg_vector *v2,suNg_vector *v3);
void cross_prod_flt(suNg_vector_flt *v1,suNg_vector_flt *v2,suNg_vector_flt *v3);
void project_to_suNg(suNg *u);
void project_to_suNg_flt(suNg_flt *u);
void project_cooling_to_suNg(suNg* g_out, suNg* g_in, int cooling);

void assign_u2ud(void);
void assign_ud2u(void);
void assign_ud2u_f(void);

/* void assign_s2sd(int len, suNf_spinor *out, suNf_spinor_flt *in); */
/* void assign_sd2s(int len, suNf_spinor_flt *out, suNf_spinor *in); */

void assign_s2sd(spinor_field *out, spinor_field_flt *in);
void assign_sd2s(spinor_field_flt *out, spinor_field *in);


/* HYP smearing */
void spatialHYP_smearing(suNg_field* out, suNg_field* in, double weight[3]);
void HYP_smearing(suNg_field* out, suNg_field* in, double weight[3]);
double min_tplaq(suNg_field* g);
void HYP_span_parameters(double mtp[6859]);
int HYP_best_parameters(double mtp[6859], double w[3]);


/* Timing */
#include <sys/time.h>
int timeval_subtract (struct timeval *result, struct timeval *x, struct timeval *y);

#endif 
