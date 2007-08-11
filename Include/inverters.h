#ifndef INVERTERS_H
#define INVERTERS_H

#include "suN_types.h"
#include "complex.h"


typedef void (*spinor_operator)(suNf_spinor *out, suNf_spinor *in);
typedef void (*spinor_operator_dble)(suNf_spinor_dble *out, suNf_spinor_dble *in);

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

int g5QMR_mshift(mshift_par *par, spinor_operator M, suNf_spinor *in, suNf_spinor_dble **out);
int g5QMR_mshiftd(mshift_par *par, spinor_operator_dble M, suNf_spinor_dble *in, suNf_spinor_dble **out);

int MINRES_mshift(mshift_par *par, spinor_operator M, suNf_spinor *in, suNf_spinor **out);

typedef struct _MINRES_par {
  double err2; /* maximum error on the solutions */
  int max_iter; /* maximum number of iterations: 0 => infinity */
} MINRES_par;
int MINRES(MINRES_par *par, spinor_operator M, suNf_spinor *in, suNf_spinor *out, suNf_spinor *trial);

int eva(int vol,int nev,int nevt,int init,int kmax,
               int imax,float ubnd,float omega1,float omega2,
               spinor_operator Op,
               suNf_spinor *ws[],suNf_spinor *ev[],float d[],int *status);

void jacobi1(int n,float a[],float d[],float v[]);
void jacobi2(int n,complex a[],float d[],complex v[]);


#endif
