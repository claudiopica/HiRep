/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *   
* All rights reserved.                                                      * 
\***************************************************************************/

#ifndef IO_H
#define IO_H
#include "input_par.h"

void read_gauge_field(char filename[]);
void write_gauge_field(char filename[]);

void read_input(input_record_t irec[], char *filename);

void read_gauge_field_single_biagio(char filename[]);

#endif
