/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *   
* All rights reserved.                                                      * 
\***************************************************************************/

#ifndef REPRESENTATION_H
#define REPRESENTATION_H

#include "suN_types.h"
#include "suN_repr_func.h"

#ifdef __cplusplus
	extern "C" {
#endif

void _group_represent2(suNf* v, suNg *u);
void _group_represent2_flt(suNf_flt* v, suNg_flt *u);
void represent_gauge_field();
//void represent_gauge_field_measure(); //TODO: not defined in libhr

#ifdef __cplusplus
	}
#endif
#endif
