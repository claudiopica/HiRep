/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *   
* All rights reserved.                                                      * 
\***************************************************************************/

#ifndef REPRESENTATION_H
#define REPRESENTATION_H

/*
#ifdef ADJOINT_FERMIONS
#elif defined SYMMETRIC_FERMIONS
#elif defined ANTISYMMETRIC_FERMIONS
#else
#endif
*/

#include "suN_repr_func.h"

#ifdef __cplusplus
extern "C" {
#endif


void _group_represent2(suNf* v, suNg *u);
void _group_represent2_flt(suNf_flt* v, suNg_flt *u);


void represent_gauge_field_cpu();
void represent_gauge_field_flt_cpu();

#ifdef WITH_GPU
void represent_gauge_field();
void represent_gauge_field_flt();
#else
extern void (*represent_gauge_field)();
extern void (*represent_gauge_field_flt)();
#endif

#ifdef __cplusplus
}
#endif //__cplusplus


#endif
