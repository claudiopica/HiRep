/***************************************************************************\
* Copyright (c) 2023, Sofie Martins                                         *   
* All rights reserved.                                                      * 
\***************************************************************************/

/**
 * @file single_double_utils_gpu.cu
 * @brief Function for conversion from single to double precision and 
 *        vice versa, GPU version
*/

#include "libhr_core.h"
#include "utils.h"
#include "./single_double_utils_gpu_kernels.hpp"

#ifdef WITH_GPU

void assign_ud2u_gpu(void) {
    const int N = u_gauge->type->gsize_gauge;
    const int grid_size = (N - 1) / BLOCK_SIZE_LINEAR_ALGEBRA + 1;
    assign_ud2u_kernel<<<grid_size, BLOCK_SIZE_LINEAR_ALGEBRA, 0, 0>>>(u_gauge_flt->gpu_ptr, u_gauge->gpu_ptr, N);
    CudaCheckError();
}

void assign_u2ud_gpu(void) {
    const int N = u_gauge_flt->type->gsize_gauge;
    const int grid_size = (N - 1) / BLOCK_SIZE_LINEAR_ALGEBRA + 1;
    assign_u2ud_kernel<<<grid_size, BLOCK_SIZE_LINEAR_ALGEBRA, 0, 0>>>(u_gauge->gpu_ptr, u_gauge_flt->gpu_ptr, N);
    CudaCheckError();
}

void assign_ud2u_f_gpu(void) {
    const int N = u_gauge_f->type->gsize_gauge;
    const int grid_size = (N - 1) / BLOCK_SIZE_LINEAR_ALGEBRA + 1;
    assign_ud2u_f_kernel<<<grid_size, BLOCK_SIZE_LINEAR_ALGEBRA, 0, 0>>>(u_gauge_f_flt->gpu_ptr, u_gauge_f->gpu_ptr, N);
    CudaCheckError();
}

void assign_u2ud_f_gpu(void) {
    const int N = u_gauge_f_flt->type->gsize_gauge;
    const int grid_size = (N - 1) / BLOCK_SIZE_LINEAR_ALGEBRA + 1;
    assign_u2ud_f_kernel<<<grid_size, BLOCK_SIZE_LINEAR_ALGEBRA, 0, 0>>>(u_gauge_f->gpu_ptr, u_gauge_f_flt->gpu_ptr, N);
    CudaCheckError();
}

void assign_s2sd_gpu(spinor_field *out, spinor_field_flt *in) {
    const int N = in->type->gsize_spinor;
    const int grid_size = (N - 1) / BLOCK_SIZE_LINEAR_ALGEBRA + 1;
    assign_s2sd_kernel<<<grid_size, BLOCK_SIZE_LINEAR_ALGEBRA, 0, 0>>>(out->gpu_ptr, in->gpu_ptr, N);
    CudaCheckError();
}

void assign_sd2s_gpu(spinor_field_flt *out, spinor_field *in) {
    const int N = in->type->gsize_spinor;
    const int grid_size = (N - 1) / BLOCK_SIZE_LINEAR_ALGEBRA + 1;
    assign_sd2s_kernel<<<grid_size, BLOCK_SIZE_LINEAR_ALGEBRA, 0, 0>>>(out->gpu_ptr, in->gpu_ptr, N);
    CudaCheckError();
}

void add_assign_s2sd_gpu(spinor_field *out, spinor_field_flt *in) {
    const int N = in->type->gsize_spinor;
    const int grid_size = (N - 1) / BLOCK_SIZE_LINEAR_ALGEBRA + 1;
    add_assign_s2sd_kernel<<<grid_size, BLOCK_SIZE_LINEAR_ALGEBRA, 0, 0>>>(out->gpu_ptr, in->gpu_ptr, N);
    CudaCheckError();
}

void add_assign_sd2s_gpu(spinor_field_flt *out, spinor_field *in) {
    const int N = in->type->gsize_spinor;
    const int grid_size = (N - 1) / BLOCK_SIZE_LINEAR_ALGEBRA + 1;
    add_assign_sd2s_kernel<<<grid_size, BLOCK_SIZE_LINEAR_ALGEBRA, 0, 0>>>(out->gpu_ptr, in->gpu_ptr, N);
    CudaCheckError();
}

void (*assign_ud2u)(void) = assign_ud2u_gpu;
void (*assign_u2ud)(void) = assign_u2ud_gpu;
void (*assign_ud2u_f)(void) = assign_ud2u_f_gpu;
void (*assign_u2ud_f)(void) = assign_u2ud_f_gpu;
void (*assign_s2sd)(spinor_field *out, spinor_field_flt *in) = assign_s2sd_gpu;
void (*assign_sd2s)(spinor_field_flt *out, spinor_field *in) = assign_sd2s_gpu;
void (*add_assign_s2sd)(spinor_field *out, spinor_field_flt *in) = add_assign_s2sd_gpu;
void (*add_assign_sd2s)(spinor_field_flt *out, spinor_field *in) = add_assign_sd2s_gpu;
#endif