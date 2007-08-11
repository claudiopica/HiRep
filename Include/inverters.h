#ifndef INVERTERS_H
#define INVERTERS_H

#include "suN_types.h"
#include "complex.h"


typedef void (*spinor_operator)(suNf_spinor *out, suNf_spinor *in);
typedef void (*spinor_operator_flt)(suNf_spinor_flt *out, suNf_spinor_flt *in);

typedef struct _mshift_par {
   int n; /* number of shifts */
   double *shift;
   double err2; /* relative error of the solutions */
   int max_iter; /* maximum number of iterations: 0 => infinity */
	 void *add_par; /* additional parameters for specific inverters */
} mshift_par;

/*
 * performs the multi-shifted CG inversion:
 * out[i] = (M-(par->shift[i]))^-1 in
 * returns the number of cg iterations done.
 */
int cg_mshift(mshift_par *par, spinor_operator M, suNf_spinor *in, suNf_spinor **out);

int BiCGstab_mshift(mshift_par *par, spinor_operator M, suNf_spinor *in, suNf_spinor **out);
int HBiCGstab_mshift(mshift_par *par, spinor_operator M, suNf_spinor *in, suNf_spinor **out);

int g5QMR_mshift(mshift_par *par, spinor_operator M, suNf_spinor *in, suNf_spinor **out);
/*int g5QMR_mshift_flt(mshift_par *par, spinor_operator_flt M, suNf_spinor_flt *in, suNf_spinor_flt **out); */

int MINRES_mshift(mshift_par *par, spinor_operator M, suNf_spinor *in, suNf_spinor **out);

typedef struct _MINRES_par {
  double err2; /* maximum error on the solutions */
  int max_iter; /* maximum number of iterations: 0 => infinity */
} MINRES_par;
int MINRES(MINRES_par *par, spinor_operator M, suNf_spinor *in, suNf_spinor *out, suNf_spinor *trial);

int eva(int vol,int nev,int nevt,int init,int kmax,
               int imax,double ubnd,double omega1,double omega2,
               spinor_operator Op,
               suNf_spinor *ws[],suNf_spinor *ev[],double d[],int *status);

void jacobi1(int n,double a[],double d[],double v[]);
void jacobi2(int n,complex a[],double d[],complex v[]);

int eva_flt(int vol,int nev,int nevt,int init,int kmax,
               int imax,float ubnd,float omega1,float omega2,
               spinor_operator_flt Op,
               suNf_spinor_flt *ws[],suNf_spinor_flt *ev[],float d[],int *status);

void jacobi1_flt(int n,float a[],float d[],float v[]);
void jacobi2_flt(int n,complex_flt a[],float d[],complex_flt v[]);


#endif
