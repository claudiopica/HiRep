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

// This is used in Cphi_inv_ which is duplicated
// so we need the CPU version of this to be visible,
// mostly for testing.
void compute_ldl_decomp_cpu(double);
void compute_clover_term_cpu();
void set_csw_cpu(double *);

extern double (*get_csw)();
extern void (*compute_ldl_decomp)(double);
extern void (*compute_clover_term)();
extern void (*clover_la_logdet)(double, double, scalar_field *);
extern void (*compute_force_logdet)(double, double);
extern void (*clover_init)(double);
extern void (*set_csw)(double *);

#ifdef __cplusplus
}
#endif
#endif
