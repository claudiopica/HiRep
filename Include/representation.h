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

void _group_represent2(suNf* v, suNg *u);
void represent_gauge_field();
void represent_gauge_field_dble();


#endif
