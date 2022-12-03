/**
 * @file cpu_geometry.h
 * @brief This file contains macros to load elements of single sites of a field.
 */

#ifndef CPU_GEOMETRY_H
#define CPU_GEOMETRY_H

#include "field_ordering.h"

#define _FIELD_AT(s,i) (((s)->ptr) + i - (s)->type->master_shift)
#define _4FIELD_AT(s,i,mu) (((s)->ptr) + coord_to_index(i-(s)->type->master_shift,mu))
#define _6FIELD_AT(s,i,mu) (((s)->ptr) + (( i - (s)->type->master_shift)*6+mu))
#define _DFIELD_AT(s,i,mu,size) (size == 1) ? _FIELD_AT(s,i) : \
				((size == 4) ? _4FIELD_AT(s,i,mu) : \
				 _6FIELD_AT(s,i,mu))

#endif
