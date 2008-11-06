/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *   
* All rights reserved.                                                      * 
\***************************************************************************/

#ifndef RATIONAL_FUNCTIONS_H
#define RATIONAL_FUNCTIONS_H

#include "inverters.h"
#include "suN.h"

typedef struct _rational_app {
	int n,d; /* approximated function = x^(n/d) */
  unsigned int order;
  double rel_error;
  double *a;
  double *b;
	unsigned int _asize; /* allocated space for arrays */
} rational_app;


void r_app_alloc(rational_app *app);
void r_app_free(rational_app *app);
void r_app_realloc(rational_app *app);
void r_app_rescale(rational_app *app, double k);
void r_app_set(rational_app *app, double min, double max);

void rational_func(rational_app *coef, spinor_operator Q, spinor_field *out, spinor_field *in);


#endif
