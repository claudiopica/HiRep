#include "libhr_core.h"

#ifdef WITH_GPU

quad_double *alloc_quad_double_sum_field(int n) {
    static quad_double *res = NULL;
    static int n_size = 0;
    if (n > n_size && res != NULL) {
        cudaFree(res);
        res = NULL;
    }

    if (res == NULL) {
        cudaMalloc((void **)&res, n * sizeof(quad_double));
        n_size = n;
    }
    return res;
}

double *alloc_double_sum_field(int n) {
    static double *res = NULL;
    static int n_size = 0;
    if (n > n_size && res != NULL) {
        cudaFree(res);
        res = NULL;
    }

    if (res == NULL) {
        cudaMalloc((void **)&res, n * sizeof(double));
        n_size = n;
    }
    return res;
}

hr_complex *alloc_complex_sum_field(int n) {
    hr_complex *res = NULL;
    int n_size = 0;
    if (n > n_size && res != NULL) {
        cudaFree(res);
        res = NULL;
    }
    if (res == NULL) {
        cudaMalloc((void **)&res, n * sizeof(hr_complex));
        n_size = n;
    }
    return res;
}

#endif