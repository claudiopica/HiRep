/******************************************************************************
*
* NOCOMPILE= !WITH_GPU
*
* Test of global_sum_gpu.c
*
******************************************************************************/

#include "libhr.h"

int reference_int(const int *vector, const int size) {
    int res = 0;
    for (int i = 0; i < size; i++) {
        res += vector[i];
    }
    return res;
}

double reference_double(const double *vector, const int size) {
    double res = 0;
    for (int i = 0; i < size; i++) {
        res += vector[i];
    }
    return res;
}

hr_complex reference_complex(const hr_complex *vector, const int size) {
    hr_complex res = 0;
    for (int i = 0; i < size; i++) {
        res += vector[i];
    }
    return res;
}

int main() {
    int size = 100;
    int vector_int[size];
    double vector_double[size];
    hr_complex vector_complex[size];
    int res_int, ref_int;
    double res_double, ref_double;
    hr_complex res_complex, ref_complex;
    for (int i = 0; i < size; i++) {
        vector_int[i] = i;
        vector_double[i] = i * 1.0;
        vector_complex[i] = i * 1.0 + (i * 1.0 + 1.0) * I;
    }

    res_int = global_sum_gpu_int(vector_int, size);
    ref_int = reference_int(vector_int, size);
    printf("%d %d\n", res_int, ref_int);

    res_double = global_sum_gpu_double(vector_double, size);
    ref_double = reference_double(vector_double, size);
    printf("%f %f\n", res_double, ref_double);

    res_complex = global_sum_gpu_complex(vector_complex, size);
    ref_complex = reference_complex(vector_complex, size);
    printf("%f+i*%f %f+i*%f\n", creal(res_complex), cimag(res_complex), creal(ref_complex), cimag(ref_complex));
    return 0;
}
