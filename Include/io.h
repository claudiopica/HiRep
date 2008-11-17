/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *   
* All rights reserved.                                                      * 
\***************************************************************************/

#ifndef IO_H
#define IO_H
#include "input_par.h"

void read_gauge_field(char filename[]);
void write_gauge_field(char filename[]);

/* void write_gauge_field_eo_lexi(char filename[]); */
void read_gauge_field_for_henty(char filename[]);

#include "spinor_field.h"
/* void write_spinor_field_eo_lexi(char filename[],spinor_field *sp); */
void write_spinor_field(char filename[],spinor_field *sp);
void read_spinor_field(char filename[],spinor_field *sp);


void write_ranlxd_state(char filename[]);
void read_ranlxd_state(char filename[]);

void read_input(input_record_t irec[], char *filename);

void read_gauge_field_single_biagio(char filename[]);

#endif
