#ifndef _MATRIX_ALG_AR_
#define _MATRIX_ALG_AR_

void lubksb(double *a, int n, int *indx, double *b);

void ludcmp(double *a, int n, int *indx, double *d);

double Det(double *a, int n);

double MinEigval(double *a, int n);

void inverse(double *a, double *ainv, int n);

int * orders(const int n,const double * vec);

void sorting(const int n,const double * vec, double * ordered_vec);

void quicksort ( int n, double *v);

void elmhes(double *a, int n);

void hqr(double *a, int n, double * wr, double * wi);

void getvec(double *a, double *EIG, double *EIGVEC, int dim, int max_ev);

void get_eigenvalues(double *mat, int mat_size, double *evalreal);

void get_eigenvectors(double *mat, int mat_size, int max_ev, double *evalreal, double *eigenvec);

#endif
