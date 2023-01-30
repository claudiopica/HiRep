/***************************************************************************\
* Copyright (c) 2023, Sofie Martins, Claudio Pica                           *
* All rights reserved.                                                      *
\***************************************************************************/

#ifndef GENERICS_H
#define GENERICS_H

#include <string.h>

// Construct suffixed function names
#define CONCAT(_name, _suffix) _name ## _suffix
#define _F_NAME(_name, _suffix) CONCAT(_name,_suffix)
#define _FUNC(_type, _name, _suffix, _args) _type _F_NAME(_name,_suffix) _args

// Construct suffixed function names with additional _gpu suffix
#define _CONCAT_GPU(_name,_suffix) _name ## _suffix ## _gpu
#define _GPU_F_NAME(_name,_suffix) _CONCAT_GPU(_name,_suffix)
#define _GPU_FUNC(_type, _name, _suffix, _args) _type _GPU_F_NAME(_name,_suffix) _args

// Construct number of sites from geom type
#define _N_SITES(_geom) CONCAT(f->type->gsize_,_geom)
#define _NUMBER_OF_SITES(_geom) _N_SITES(_geom)

// Construct number of buffers from geom type
#define _N_BUFFERS(_geom) CONCAT(type->nbuffers_,_geom)
#define _NUMBER_OF_BUFFERS(_geom) _N_BUFFERS(_geom)

#endif

