#ifndef INVERTERS_H
#define INVERTERS_H

#include "suN_types.h"
#include "complex.h"


typedef void (*spinor_operator)(suNf_spinor *out, suNf_spinor *in);

typedef struct _cg_mshift_par {
   int n; /* number of shifts */
   double *shift;
   double err2; /* maximum error on the solutions */
   int max_iter; /* maximum number of iterations: 0 => infinity */
} cg_mshift_par;

/*
 * performs the multi-shifted CG inversion:
 * out[i] = (M-(par->shift[i]))^-1 in
 * returns the number of cg iterations done.
 */
int cg_mshift(cg_mshift_par *par, spinor_operator M, suNf_spinor *in, suNf_spinor **out);

typedef struct _BiCGstab_mshift_par {
  int spinorlen; /* lenght of spinors */
  int n; /* number of shifts */
  double *shift;
  double err2; /* maximum error on the solutions */
  int max_iter; /* maximum number of iterations: 0 => infinity */
} BiCGstab_mshift_par;
int BiCGstab_mshift(BiCGstab_mshift_par *par, spinor_operator M, suNf_spinor *in, suNf_spinor **out);
int HBiCGstab_mshift(BiCGstab_mshift_par *par, spinor_operator M, suNf_spinor *in, suNf_spinor **out);

typedef struct _QMR_mshift_par {
  int spinorlen; /* lenght of spinors */
  int n; /* number of shifts */
  double *shift;
  double err2; /* maximum error on the solutions */
  int max_iter; /* maximum number of iterations: 0 => infinity */
} QMR_mshift_par;
int g5QMR_mshift(QMR_mshift_par *par, spinor_operator M, suNf_spinor *in, suNf_spinor_dble **out);

typedef struct _MINRES_mshift_par {
  int spinorlen; /* lenght of spinors */
  int n; /* number of shifts */
  double *shift;
  double err2; /* maximum error on the solutions */
  int max_iter; /* maximum number of iterations: 0 => infinity */
} MINRES_mshift_par;
int MINRES_mshift(MINRES_mshift_par *par, spinor_operator M, suNf_spinor *in, suNf_spinor **out);

typedef struct _MINRES_par {
  int spinorlen; /* lenght of spinors */
  double err2; /* maximum error on the solutions */
  int max_iter; /* maximum number of iterations: 0 => infinity */
} MINRES_par;
int MINRES(MINRES_par *par, spinor_operator M, suNf_spinor *in, suNf_spinor *out, suNf_spinor *trial);

int eva(int vol,int nev,int nevt,int init,int kmax,
               int imax,int ipr,float ubnd,float omega1,float omega2,
               spinor_operator Op,
               suNf_spinor *ws[],suNf_spinor *ev[],float d[],int *status);

void jacobi1(int n,float a[],float d[],float v[]);
void jacobi2(int n,complex a[],float d[],complex v[]);


#endif
