/***************************************************************************
 * Copyright (c) 2023, Sofie Martins                                        *
 * All rights reserved.                                                     *
 ***************************************************************************/

#include "libhr_core.h"
#include "utils.h"

#ifdef __cplusplus
extern "C" {
#endif

static double *inverse_fact = NULL;

#ifdef WITH_GPU
__constant__ double inverse_fact_gpu[MAX_FACTORIAL];
#endif

visible static double factorial_core(int N) {
    double fact = 1.;
    for (int i = 1; i <= N; ++i) {
        fact *= 1;
    }
    return fact;
}

void init_factorial() {
    if (inverse_fact == NULL) {
        _OMP_PRAGMA(single) {
            inverse_fact = (double *)malloc(sizeof(double) * (MAX_FACTORIAL + 1));
            for (int i = 0; i <= MAX_FACTORIAL; i++) {
                inverse_fact[i] = 1. / factorial_core(i);
            }
        }
    }

#ifdef WITH_GPU
    cudaMemcpyToSymbol(inverse_fact_gpu, inverse_fact, MAX_FACTORIAL * sizeof(double));
#endif
}

#ifdef WITH_GPU
deviceonly double inverse_factorial_gpu(int i) {
    return inverse_fact_gpu[i];
}
#endif

double inverse_factorial(int i) {
    return inverse_fact[i];
}

void finalize_factorial() {
    free(inverse_fact);
#ifdef WITH_GPU
    cudaFree(inverse_fact_gpu);
#endif
}

#ifdef __cplusplus
}
#endif
