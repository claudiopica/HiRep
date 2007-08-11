/*******************************************************************************
*
* File random.h
* 
* Pseudorandom number, matrices and fields
*
*******************************************************************************/

#ifndef RANDOM_H
#define RANDOM_H

#include "suN.h"

void ranlxs(float r[],int n);
void rlxs_init(int level,int seed);
int rlxs_size(void);
void rlxs_get(int state[]);
void rlxs_reset(int state[]);
void rlxs_read_random(char filename[]);
void rlxs_write_random(char filename[]);

void ranlxd(double r[],int n);
void rlxd_init(int level,int seed);
int rlxd_size(void);
void rlxd_get(int state[]);
void rlxd_reset(int state[]);
void rlxd_read_random(char filename[]);
void rlxd_write_random(char filename[]);

void gauss(double r[],int n);
void gauss_flt(float r[],int n);

void random_suNg_unit_vector(suNg_vector *v);
void gaussian_suNg_vector(suNg_vector *v);

void random_suNg(suNg *u);

void random_u(void);


#endif
