/***************************************************************************\
* Copyright (c) 2022 Sofie Martins                                          *   
* All rights reserved.                                                      * 
\***************************************************************************/

#include "geometry.h"
#include "libhr_core.h"
#include <string.h>

void init_neighbors_gpu() {
#ifdef WITH_GPU
    CudaCheckError();

    int N = glattice.gsize_gauge;

    CHECK_CUDA(cudaMalloc((void **)&iup_gpu, 4 * N * sizeof(int)));
    CHECK_CUDA(cudaMalloc((void **)&idn_gpu, 4 * N * sizeof(int)));
    CHECK_CUDA(cudaMalloc((void **)&imask_gpu, N * sizeof(char)));
    CHECK_CUDA(cudaMalloc((void **)&ipt_gpu,
                          (X + 2 * X_BORDER) * (Y + 2 * Y_BORDER) * (Z + 2 * Z_BORDER) * (T + 2 * T_BORDER) * sizeof(int)));
    CHECK_CUDA(cudaMemcpy(iup_gpu, iup, 4 * N * sizeof(int), cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpy(idn_gpu, idn, 4 * N * sizeof(int), cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpy(imask_gpu, imask, N * sizeof(*imask), cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpy(ipt_gpu, ipt,
                          (X + 2 * X_BORDER) * (Y + 2 * Y_BORDER) * (Z + 2 * Z_BORDER) * (T + 2 * T_BORDER) * sizeof(int),
                          cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpyToSymbol(&T_EXT_GPU, &T_EXT, sizeof(int), 0, cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpyToSymbol(&X_EXT_GPU, &X_EXT, sizeof(int), 0, cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpyToSymbol(&Y_EXT_GPU, &Y_EXT, sizeof(int), 0, cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpyToSymbol(&Z_EXT_GPU, &Z_EXT, sizeof(int), 0, cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpyToSymbol(&T_GPU, &T, sizeof(int), 0, cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpyToSymbol(&X_GPU, &X, sizeof(int), 0, cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpyToSymbol(&Y_GPU, &Y, sizeof(int), 0, cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpyToSymbol(&Z_GPU, &Z, sizeof(int), 0, cudaMemcpyHostToDevice));
#ifdef FERMION_THETA
    CHECK_CUDA(cudaMemcpyToSymbol(&eitheta_gpu[0], &eitheta[0], sizeof(hr_complex) * 4, 0, cudaMemcpyHostToDevice));
#endif

    box_t *L = geometryBoxes;
    int number_of_boxes = 0;
    do {
        number_of_boxes++;
    } while ((L = L->next));

    CHECK_CUDA(cudaMalloc((void **)&geometryBoxes_gpu, number_of_boxes * sizeof(box_t)));

    L = geometryBoxes;
    int i = 0;
    do {
        CHECK_CUDA(cudaMemcpy(&geometryBoxes_gpu[i], L, sizeof(box_t), cudaMemcpyHostToDevice));
        ++i;
    } while ((L = L->next));

    // adapted from new_geom.c
    size_t main_mem_volume, buf_mem_volume;
    geometryMemSize(geometryBoxes, &main_mem_volume, &buf_mem_volume);
    CHECK_CUDA(cudaMalloc((void **)&icoord_gpu, (main_mem_volume + buf_mem_volume) * sizeof(coord4)));
    CHECK_CUDA(cudaMemcpy(icoord_gpu, geometryBoxes->icoord, (main_mem_volume + buf_mem_volume) * sizeof(coord4),
                          cudaMemcpyHostToDevice));
    sb_icoord_gpu = icoord_gpu + main_mem_volume;

    CHECK_CUDA(cudaStreamCreate(&non_default_stream));

    input = (kernel_field_input **)malloc(glattice.nbuffers_spinor * sizeof(kernel_field_input) / 2);
    for (int i = 0; i < glattice.nbuffers_spinor / 2; i++) {
        CHECK_CUDA(cudaMalloc((void **)&input[i], sizeof(kernel_field_input)));
    }
#endif
}
