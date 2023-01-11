/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *   
* All rights reserved.                                                      * 
\***************************************************************************/

#ifndef RATIONAL_FUNCTIONS_H
#define RATIONAL_FUNCTIONS_H

#include "Inverters/linear_solvers.h"
#include "spinor_field.h"

#ifdef __cplusplus
	extern "C" {
#endif

typedef struct rational_app {
  /* approximate function = x^(n/d) */
  /* in the range [min,max] with relative precision rel_error */
	int n,d;
  unsigned int order;
  double rel_error;
  double *a;
  double *b;
  double min, max; /* eigenvalues */
	unsigned int _asize; /* allocated space for arrays */
} rational_app;


void r_app_alloc(rational_app *app);
void r_app_free(rational_app *app);
void r_app_realloc(rational_app *app);
void r_app_rescale(rational_app *app, double k);
void r_app_set(rational_app *app, double min, double max);

void rational_func(rational_app *coef, spinor_operator Q, spinor_field *out, spinor_field *in);

#ifdef __cplusplus
	}
#endif
#endif
