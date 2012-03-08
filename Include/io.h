/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *   
* All rights reserved.                                                      * 
\***************************************************************************/

#ifndef IO_H
#define IO_H
#include "input_par.h"
#include "spinor_field.h"
#include <stdio.h>

int fwrite_BE_int(int* ptr, size_t n, FILE* fp);
int fwrite_LE_int(int* ptr, size_t n, FILE* fp);
int fwrite_BE_double(double* ptr, size_t n, FILE* fp);
int fwrite_LE_double(double* ptr, size_t n, FILE* fp);

int fread_BE_int(int* ptr, size_t n, FILE* fp);
int fread_LE_int(int* ptr, size_t n, FILE* fp);
int fread_BE_double(double* ptr, size_t n, FILE* fp);
int fread_LE_double(double* ptr, size_t n, FILE* fp);

void read_gauge_field(char filename[]);
void write_gauge_field(char filename[]);

void write_ranlxd_state(char filename[]);
void read_ranlxd_state(char filename[]);

void read_input(input_record_t irec[], char *filename);

void read_spinor_field_ascii(char filename[],spinor_field * sf);

#endif
