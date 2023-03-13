#ifndef SINGLE_DOUBLE_UTILS_GPU_KERNELS_HPP
#define SINGLE_DOUBLE_UTILS_GPU_KERNELS_HPP

#ifdef WITH_GPU

__global__ void assign_ud2u_kernel(suNg_flt *gauge_flt, suNg *gauge, int N) {
    double d;
    float *f;
    const size_t n_doubles = 4 * N * sizeof(suNg) / sizeof(double);
    for (int ix = blockIdx.x * blockDim.x + threadIdx.x; ix < n_doubles; ix += blockDim.x * gridDim.x) {
        d = *(((double *)gauge) + ix);
        f = ((float *)gauge_flt) + ix;
        *f = (float)d;
    }
}

__global__ void assign_u2ud_kernel(suNg *gauge, suNg_flt *gauge_flt, int N) {
    double *d;
    float f;
    const size_t n_floats = 4 * N * sizeof(suNg_flt) / sizeof(float);
    for (int ix = blockIdx.x * blockDim.x + threadIdx.x; ix < n_floats; ix += blockDim.x * gridDim.x) {
        f = *(((float *)gauge_flt) + ix);
        d = ((double *)gauge) + ix;
        *d = (double)f;
    }
}

__global__ void assign_ud2u_f_kernel(suNf_flt *gauge_flt, suNf *gauge, int N) {
    double d;
    float *f;
    const size_t n_doubles = 4 * N * sizeof(suNf) / sizeof(double);
    for (int ix = blockIdx.x * blockDim.x + threadIdx.x; ix < n_doubles; ix += blockDim.x * gridDim.x) {
        d = *(((double *)gauge) + ix);
        f = ((float *)gauge_flt) + ix;
        *f = (float)d;
    }
}

__global__ void assign_u2ud_f_kernel(suNf *gauge, suNf_flt *gauge_flt, int N) {
    double *d;
    float f;
    const size_t n_floats = 4 * N * sizeof(suNf_flt) / sizeof(float);
    for (int ix = blockIdx.x * blockDim.x + threadIdx.x; ix < n_floats; ix += blockDim.x * gridDim.x) {
        f = *(((float *)gauge_flt) + ix);
        d = ((double *)gauge) + ix;
        *d = (double)f;
    }
}

__global__ void assign_s2sd_kernel(suNf_spinor *out, suNf_spinor_flt *in, int N) {
    float f;
    double *d;
    const size_t n_floats = N * sizeof(*in) / sizeof(float);
    for (int ix = blockIdx.x * blockDim.x + threadIdx.x; ix < n_floats; ix += blockDim.x * gridDim.x) {
        f = *(((float *)in) + ix);
        d = ((double *)out) + ix;
        *d = (double)f;
    }
}

__global__ void assign_sd2s_kernel(suNf_spinor_flt *out, suNf_spinor *in, int N) {
    double d;
    float *f;
    const size_t n_floats = N * sizeof(*in) / sizeof(double);
    for (int ix = blockIdx.x * blockDim.x + threadIdx.x; ix < n_floats; ix += blockDim.x * gridDim.x) {
        d = *((double *)in + ix);
        f = (float *)out + ix;
        *f = (float)d;
    }
}

__global__ void add_assign_s2sd_kernel(suNf_spinor *out, suNf_spinor_flt *in, int N) {
    float f;
    double *d;
    const size_t n_floats = N * sizeof(*in) / sizeof(float);
    for (int ix = blockIdx.x * blockDim.x + threadIdx.x; ix < n_floats; ix += blockDim.x * gridDim.x) {
        f = *(((float *)in) + ix);
        d = ((double *)out) + ix;
        *d += (double)f;
    }
}

__global__ void add_assign_sd2s_kernel(suNf_spinor_flt *out, suNf_spinor *in, int N) {
    double d;
    float *f;
    const size_t n_floats = N * sizeof(*in) / sizeof(double);
    for (int ix = blockIdx.x * blockDim.x + threadIdx.x; ix < n_floats; ix += blockDim.x * gridDim.x) {
        d = *((double *)in + ix);
        f = (float *)out + ix;
        *f += (float)d;
    }
}

#endif
#endif