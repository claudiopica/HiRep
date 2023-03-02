#ifndef JACKNIFE_H
#define JACKNIFE_H
#ifdef __cplusplus
extern "C" {
#endif

double auto_corr_time(int n, int tmax, double g[], int *flag);
double sigma_bin(int n, int binsize, double a[]);
double sigma_replicas(int n, int r, double a[], double *tau, int *flag);
double sigma_jackknife(int nobs, int n, double a[], double *ave_j, double (*pobs)(double v[]));

#ifdef __cplusplus
}
#endif
#endif //JACKNIFE_H
