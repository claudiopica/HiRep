/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *   
* All rights reserved.                                                      * 
\***************************************************************************/

#ifndef ARCHIVE_H
#define ARCHIVE_H

#include <stdio.h>
#include "spinor_field.h"

#ifdef __cplusplus
  extern "C" {
#endif

// endian.h
int fwrite_BE_int(int* ptr, size_t n, FILE* fp);
int fwrite_LE_int(int* ptr, size_t n, FILE* fp);
int fwrite_BE_double(double* ptr, size_t n, FILE* fp);
int fwrite_LE_double(double* ptr, size_t n, FILE* fp);

int fread_BE_int(int* ptr, size_t n, FILE* fp);
int fread_LE_int(int* ptr, size_t n, FILE* fp);
int fread_BE_double(double* ptr, size_t n, FILE* fp);
int fread_LE_double(double* ptr, size_t n, FILE* fp);
int fread_BE_float(float* ptr, size_t n, FILE* fp);

//archive.h
void read_gauge_field(char filename[]);
void write_gauge_field(char filename[]);
void read_scalar_field(char filename[]);
void write_scalar_field(char filename[]);
void read_gauge_field_matrix(char filename[]);
void write_gauge_field_matrix(char filename[]);
void write_ranlxd_state(char filename[]);
void read_ranlxd_state(char filename[]);

//archive_su2quat.h
void read_gauge_field_su2(char filename[]);
void read_gauge_field_su2q(char filename[]);
void write_gauge_field_su2q(char filename[]);


//TODO: check this
void read_spinor_field_ascii(char filename[],spinor_field * sf); //this does not exist in library
void read_gauge_field_nocheck(char filename[]); //this does not exist in library
void print_mat(suNg* mat); //TODO: does not exist in library

#ifdef __cplusplus
  }
#endif
#endif
