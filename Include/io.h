/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *   
* All rights reserved.                                                      * 
\***************************************************************************/

#ifndef IO_H
#define IO_H
#include "input_par.h"
#include "spinor_field.h"
#include <stdio.h>
#include "suN.h"

int fwrite_BE_int(int* ptr, size_t n, FILE* fp);
int fwrite_LE_int(int* ptr, size_t n, FILE* fp);
int fwrite_BE_double(double* ptr, size_t n, FILE* fp);
int fwrite_LE_double(double* ptr, size_t n, FILE* fp);

int fread_BE_int(int* ptr, size_t n, FILE* fp);
int fread_LE_int(int* ptr, size_t n, FILE* fp);
int fread_BE_double(double* ptr, size_t n, FILE* fp);
int fread_LE_double(double* ptr, size_t n, FILE* fp);
int fread_BE_float(float* ptr, size_t n, FILE* fp);

void read_gauge_field(char filename[]);
void read_gauge_field_nocheck(char filename[]);
void write_gauge_field(char filename[]);
void read_gauge_field_matrix(char filename[]);

#if NG==2
void read_gauge_field_su2(char filename[]);
void read_gauge_field_su2q(char filename[]);
void write_gauge_field_su2q(char filename[]);
#endif

void write_ranlxd_state(char filename[]);
void read_ranlxd_state(char filename[]);

void read_input(input_record_t irec[], char *filename);

void read_spinor_field_ascii(char filename[],spinor_field * sf);

void print_mat(suNg* mat);
#endif
