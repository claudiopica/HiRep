/***************************************************************************\
 * Copyright (c) 2024, Ari Hietanen, Ulrik Søndergaard, Sofie Martins       *
 * All rights reserved.                                                     *
\***************************************************************************/

#ifdef WITH_GPU

#include "inverters.h"
#include "libhr_core.h"
#ifndef HIP
#include <cub/cub.cuh>
#else
#include <hipcub/hipcub.hpp>
#endif

template <class T> T global_max_gpu(T *vector, int size) {
    T *maxval_d;
    T maxval;
    cudaMalloc((void **)&maxval_d, sizeof(T));
    void *d_tmp_storage = NULL;
    size_t temp_storage_bytes = 0;
    cub::DeviceReduce::Max(d_tmp_storage, temp_storage_bytes, vector, maxval_d, size);
    cudaMalloc((void **)&d_tmp_storage, temp_storage_bytes);
    cub::DeviceReduce::Max(d_tmp_storage, temp_storage_bytes, vector, maxval_d, size);
    cudaMemcpy(&maxval, maxval_d, sizeof(T), cudaMemcpyDeviceToHost);
    cudaFree(maxval_d);
    cudaFree(d_tmp_storage);
    return maxval;
}

template <class T> T global_sum_gpu(T *vector, int size) {
    T *sum_d;
    T sum;
    cudaMalloc((void **)&sum_d, sizeof(T));
    void *d_tmp_storage = NULL;
    size_t temp_storage_bytes = 0;
    cub::DeviceReduce::Sum(d_tmp_storage, temp_storage_bytes, vector, sum_d, size);
    cudaMalloc((void **)&d_tmp_storage, temp_storage_bytes);
    cub::DeviceReduce::Sum(d_tmp_storage, temp_storage_bytes, vector, sum_d, size);
    cudaMemcpy(&sum, sum_d, sizeof(T), cudaMemcpyDeviceToHost);
    cudaFree(sum_d);
    cudaFree(d_tmp_storage);
    return sum;
}

template int global_sum_gpu<int>(int *vector, int size);
template float global_sum_gpu<float>(float *vector, int size);
template double global_sum_gpu<double>(double *vector, int size);
template hr_complex_flt global_sum_gpu<hr_complex_flt>(hr_complex_flt *vector, int size);
template hr_complex global_sum_gpu<hr_complex>(hr_complex *vector, int size);

template int global_max_gpu<int>(int *vector, int size);
template float global_max_gpu<float>(float *vector, int size);
template double global_max_gpu<double>(double *vector, int size);

int global_sum_gpu_int(int *vector, int size) {
    int res;
    int *vector_d;
    cudaMalloc((void **)&vector_d, size * sizeof(int));
    cudaMemcpy(vector_d, vector, size * sizeof(int), cudaMemcpyHostToDevice);
    res = global_sum_gpu<int>(vector_d, size);
    cudaFree(vector_d);
    return res;
}

float global_sum_gpu_float(float *vector, int size) {
    float res;
    float *vector_d;
    cudaMalloc((void **)&vector_d, size * sizeof(float));
    cudaMemcpy(vector_d, vector, size * sizeof(float), cudaMemcpyHostToDevice);
    res = global_sum_gpu<float>(vector_d, size);
    cudaFree(vector_d);
    return res;
}

double global_sum_gpu_double(double *vector, int size) {
    double res;
    double *vector_d;
    cudaMalloc((void **)&vector_d, size * sizeof(double));
    cudaMemcpy(vector_d, vector, size * sizeof(double), cudaMemcpyHostToDevice);
    res = global_sum_gpu<double>(vector_d, size);
    cudaFree(vector_d);
    return res;
}

hr_complex_flt global_sum_gpu_complex_flt(hr_complex_flt *vector, int size) {
    hr_complex_flt res;
    hr_complex_flt *vector_d;
    cudaMalloc((void **)&vector_d, size * sizeof(hr_complex_flt));
    cudaMemcpy(vector_d, vector, size * sizeof(hr_complex_flt), cudaMemcpyHostToDevice);
    res = global_sum_gpu<hr_complex_flt>(vector_d, size);
    cudaFree(vector_d);
    return res;
}

hr_complex global_sum_gpu_complex(hr_complex *vector, int size) {
    hr_complex res;
    hr_complex *vector_d;
    cudaMalloc((void **)&vector_d, size * sizeof(hr_complex));
    cudaMemcpy(vector_d, vector, size * sizeof(hr_complex), cudaMemcpyHostToDevice);
    res = global_sum_gpu<hr_complex>(vector_d, size);
    cudaFree(vector_d);
    return res;
}

int global_max_gpu_int(int *vector, int size) {
    int res;
    int *vector_d;
    cudaMalloc((void **)&vector_d, size * sizeof(int));
    cudaMemcpy(vector_d, vector, size * sizeof(int), cudaMemcpyHostToDevice);
    res = global_max_gpu<int>(vector_d, size);
    cudaFree(vector_d);
    return res;
}

float global_max_gpu_float(float *vector, int size) {
    float res;
    float *vector_d;
    cudaMalloc((void **)&vector_d, size * sizeof(float));
    cudaMemcpy(vector_d, vector, size * sizeof(float), cudaMemcpyHostToDevice);
    res = global_max_gpu<float>(vector_d, size);
    cudaFree(vector_d);
    return res;
}

double global_max_gpu_double(double *vector, int size) {
    double res;
    double *vector_d;
    cudaMalloc((void **)&vector_d, size * sizeof(double));
    cudaMemcpy(vector_d, vector, size * sizeof(double), cudaMemcpyHostToDevice);
    res = global_max_gpu<double>(vector_d, size);
    cudaFree(vector_d);
    return res;
}

#endif
