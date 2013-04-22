/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *   
* All rights reserved.                                                      * 
\***************************************************************************/

#ifndef INVERTERS_H
#define INVERTERS_H

#include "suN_types.h"
#include "complex.h"
#include "spinor_field.h"

typedef void (*spinor_operator_dbl)(spinor_field *out, spinor_field *in);
typedef void (*spinor_operator_flt)(spinor_field_flt *out, spinor_field_flt *in);
typedef void (*spinor_operator_m)(spinor_field *out, spinor_field *in, double m);

typedef struct _spinor_operator {
  spinor_operator_dbl dbl;
  spinor_operator_flt flt;
} spinor_operator;

typedef struct _inverter_par {
  double err2; /* maximum error on the solutions */
  int max_iter; /* maximum number of iterations: 0 => infinity */
  double err2_flt; /* maximum error on the solutions */
  int max_iter_flt; /* maximum number of iterations: 0 => infinity */
  int kry_dim; /* Size of krylov subspace */
  int kry_dim_flex; /* Size of krylov subspace */
  int n; /* number of shifts */
  double *shift; /* array of shifts */
} inverter_par;

typedef struct _mshift_par {
   int n; /* number of shifts */
   double *shift;
   double err2; /* relative error of the solutions */
   int max_iter; /* maximum number of iterations: 0 => infinity */
	 void *add_par; /* additional parameters for specific inverters */
} mshift_par;



// We might want to add a "trial" spinor to the argument list
typedef int (*inverter_ptr)(mshift_par* par, spinor_operator M, spinor_field *in, spinor_field *out);
void SAP_prec(int nu, inverter_ptr inv, mshift_par *par, spinor_operator M, spinor_field *in, spinor_field *out);


/*
 * performs the multi-shifted CG inversion:
 * out[i] = (M-(par->shift[i]))^-1 in
 * returns the number of cg iterations done.
 */
int cg_mshift(mshift_par *par, spinor_operator M, spinor_field *in, spinor_field *out);
int cg_mshift_def(mshift_par *par, spinor_operator M, spinor_operator P, spinor_operator_m Pinv, spinor_field *in, spinor_field *out);
int cg_mshift_flt(mshift_par *par, spinor_operator M, spinor_operator_flt F, spinor_field *in, spinor_field *out);


typedef struct {
  double err2; /* maximum error on the solutions */
  int max_iter; /* maximum number of iterations: 0 => infinity */
  double err2_flt; /* maximum error on the solutions */
  int max_iter_flt; /* maximum number of iterations: 0 => infinity */
} g5QMR_fltacc_par;

int g5QMR_mshift(mshift_par *par, spinor_operator M, spinor_field *in, spinor_field *out);
int g5QMR_mshift_trunc(mshift_par *par, int trunc_iter, spinor_operator M, spinor_field *in, spinor_field *out_trunc, spinor_field *out);
int g5QMR_fltacc(g5QMR_fltacc_par *par, spinor_operator M, spinor_operator_flt M_flt, spinor_field *in, spinor_field *out);


typedef struct _MINRES_par {
  double err2; /* maximum error on the solutions */
  int max_iter; /* maximum number of iterations: 0 => infinity */
} MINRES_par;

int HBiCGstab(MINRES_par *par, spinor_operator M, spinor_field *in, spinor_field *out);
int HBiCGstab_flt(MINRES_par *par, spinor_operator_flt M, spinor_field_flt *in, spinor_field_flt *out);

/*
int BiCGstab_mshift(mshift_par *par, spinor_operator M, spinor_field *in, spinor_field *out);
int HBiCGstab_mshift(mshift_par *par, spinor_operator M, spinor_field *in, spinor_field *out);
*/


int MINRES_mshift(mshift_par *par, spinor_operator M, spinor_field *in, spinor_field *out);

int MINRES(MINRES_par *par, spinor_operator M, spinor_field *in, spinor_field *out, spinor_field *trial);
int MINRES_flt(MINRES_par *par, spinor_operator_flt M, spinor_field_flt *in, spinor_field_flt *out, spinor_field_flt *trial);


int GMRES(inverter_par *par, spinor_operator M, spinor_field *in, spinor_field *out, spinor_field *trial);
int GMRES_flt(inverter_par *par, spinor_operator M, spinor_field_flt *in, spinor_field_flt *out, spinor_field_flt *trial);

int FGMRES(inverter_par *par, spinor_operator M, spinor_field *in, spinor_field *out, spinor_field *trial, spinor_operator precon);

void DDalphaAMG(mshift_par *par, spinor_operator M, spinor_field *in, spinor_field *out);
void DDalphaAMG_setup(mshift_par *par, spinor_operator M, int N, int nu, int n_inv);
void spinor_field_to_course_spinor_field(spinor_field *big);
void coarse_spinor_field_to_spinor_field(spinor_field *big);


int eva(int nev,int nevt,int init,int kmax,
               int imax,double ubnd,double omega1,double omega2,
               spinor_operator Op,spinor_field *ws,
               spinor_field *ev,double d[],int *status);

void jacobi1(int n,double a[],double d[],double v[]);
void jacobi2(int n,complex a[],double d[],complex v[]);



void empty_buffers(spinor_field *s);

#endif
