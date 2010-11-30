#ifndef EFFECTIVE_MASS_H
#define EFFECTIVE_MASS_H

#include "bs_type.h"

double fit(int left, int right,  Corr_t *data); 

/* Effective mass defined as in JHEP 0702:082,2007 [hep-lat/0701009] */
double effm_err();

/* this function solve the equation
 * cosh[ (a+1) M ] = k cosh[ a M ]
 */
double plain_eff_mass(int a, double k);

double hc(int t, double m, int Lt);


/* Effective masses using Prony's method
 * one and two effectie masses only
 */

/* compute the binomial */
double bin(int n, int k);

/* put together the modifies correlator */
void build_y(double *c, int n, double *dst);

/* only 1 effective mass */
int centered_prony_eff_mass_1(double *y, double *m1);
/* two effective masses */
int centered_eff_mass_2(double *y, double *m1, double *m2);



int shifted_prony_eff_mass_1(double* C, int t, int tmax, double* m1, int lt);
int shifted_prony_eff_mass_2(double* C, int t, int tmax, double* m1, double *m2, int lt);


#endif

