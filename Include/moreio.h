/***************************************************************************\
* Copyright (c) 2008, Agostino Patella, Claudio Pica                        *   
* All rights reserved.                                                      * 
\***************************************************************************/

#ifndef MOREIO_H
#define MOREIO_H
#include "input_par.h"

void read_gauge_field(char filename[]);
void write_gauge_field(char filename[]);

void read_gauge_field_eolexi(char filename[]);
void write_gauge_field_eolexi(char filename[]);

void read_gauge_field_henty(char filename[]);

#include "spinor_field.h"
/* void write_spinor_field_eo_lexi(char filename[],spinor_field *sp); */
void write_spinor_field(char filename[],spinor_field *sp);
void read_spinor_field(char filename[],spinor_field *sp);

#endif
