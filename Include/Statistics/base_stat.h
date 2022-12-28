#ifndef BASE_STAT_H
#define BASE_STAT_H
#ifdef __cplusplus
	extern "C" {
#endif

double average(int n, double a[]);
double sigma0(int n, double a[]);
void auto_corr(int n, double a[], int tmax, double gamma[]);
double sigma(int n, double a[], double *tau, int *flag);

#ifdef __cplusplus
	}
#endif
#endif //BASE_STAT_H
