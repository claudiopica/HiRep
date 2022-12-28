/***************************************************************************\
* Copyright (c) 2016, Martin Hansen                                         *
* All rights reserved.                                                      * 
\***************************************************************************/

#ifndef CLOVER_TOOLS_H
#define CLOVER_TOOLS_H

#include "spinor_field.h"

#ifdef __cplusplus
    extern "C" {
#endif

#if (defined(WITH_CLOVER) && defined(UPDATE_EO))
#define WITH_CLOVER_EO
#endif

double get_csw();
void compute_ldl_decomp(double);
void compute_clover_term();
void clover_la_logdet(double, double, scalar_field*);
void compute_force_logdet(double, double);
void clover_init(double);
void set_csw(double *);

#ifdef __cplusplus
    }
#endif
#endif
