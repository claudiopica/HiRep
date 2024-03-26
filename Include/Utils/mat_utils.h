/// Headerfile for:
/// - inv_Cmplx_Ng.c
/// - det_Cmplx_Ng.c
/// - diag_hmat.c

#ifndef MAT_UTILS_H
#define MAT_UTILS_H
#ifdef __cplusplus
extern "C" {
#endif

#include "Core/gpu.h"

//inv_Cmplx_Ng.c
#ifndef GAUGE_SON
visible void ludcmp(hr_complex *a, int *indx, double *d);
visible void lubksb(hr_complex *a, int *indx, hr_complex *b, int N);
visible void inv_Cmplx_Ng(suNg *a);
#endif

//det_Cmplx_Ng.c
#ifndef GAUGE_SON
visible void det_Cmplx_Ng(hr_complex *res, suNg *a);
#else
visible void det_Cmplx_Ng(double *res, suNg *a);
#endif

//diag_hmat.c
#ifdef GAUGE_SON
void tridiagonalize(suNg *hmat, double *diag, double *roffdiag);
void tridiagonalize(suNg *hmat, double *diag, double *roffdiag);
void diag_tridiag(suNg *hmat, double *diag, double *roffdiag);
void diag_tridiag(suNg *hmat, double *diag, double *offdiag);
double pythag(double a, double b);
double pythag(double a, double b);
double sign(double a, double b);
double sign(double a, double b);
void diag_hmat(suNg *hmat, double *diag);
void diag_hmat(suNg *hmat, double *dag);
#endif

#ifdef __cplusplus
}
#endif
#endif //MAT_UTILS_H
