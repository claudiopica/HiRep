/***************************************************************************\
* Copyright (c) 2022 Sofie Martins                                          *   
* All rights reserved.                                                      * 
\***************************************************************************/

#include "geometry.h"
#include "libhr_core.h"
#include <string.h>

void init_neighbors_gpu() {
#ifdef WITH_GPU
    int N = glattice.gsize_gauge;

    cudaError_t error_id;
    error_id = cudaMalloc((void **)&iup_gpu, 4 * N * sizeof(int));
    error(error_id != cudaSuccess, 1, "init_neighbors_gpu", cudaGetErrorString(error_id));

    error_id = cudaMalloc((void **)&idn_gpu, 4 * N * sizeof(int));
    error(error_id != cudaSuccess, 1, "init_neighbors_gpu", "Error allocating idn_gpu neighbors array.\n");

    error_id = cudaMalloc((void **)&imask_gpu, N * sizeof(char));
    error(error_id != cudaSuccess, 1, "init_neighbors_gpu", "Error allocating imask_gpu lookup table.\n");

    error_id = cudaMalloc((void **)&ipt_gpu,
                          (X + 2 * X_BORDER) * (Y + 2 * Y_BORDER) * (Z + 2 * Z_BORDER) * (T + 2 * T_BORDER) * sizeof(int));
    error(error_id != cudaSuccess, 1, "init_neighbors_gpu", "Error allocating ipt_gpu lookup table.\n");

    error_id = cudaMemcpy(iup_gpu, iup, 4 * N * sizeof(int), cudaMemcpyHostToDevice);
    error(error_id != cudaSuccess, 1, "init_neighbors_gpu", "Error copying iup neighbors array to device memory.\n");

    error_id = cudaMemcpy(idn_gpu, idn, 4 * N * sizeof(int), cudaMemcpyHostToDevice);
    error(error_id != cudaSuccess, 1, "init_neighbors_gpu", "Error copying idn neighbors array to device memory.\n");

    error_id = cudaMemcpy(imask_gpu, imask, N * sizeof(*imask), cudaMemcpyHostToDevice);
    error(error_id != cudaSuccess, 1, "init_neighbors_gpu", "Error copying imask lookup table to device memory.\n");

    error_id = cudaMemcpy(ipt_gpu, ipt,
                          (X + 2 * X_BORDER) * (Y + 2 * Y_BORDER) * (Z + 2 * Z_BORDER) * (T + 2 * T_BORDER) * sizeof(int),
                          cudaMemcpyHostToDevice);
    error(error_id != cudaSuccess, 1, "init_neighbors_gpu", "Error copying ipt to device memory.\n");

    error_id = cudaMemcpyToSymbol(&T_EXT_GPU, &T_EXT, sizeof(int), 0, cudaMemcpyHostToDevice);
    error(error_id != cudaSuccess, 1, "init_neighbors_gpu", "Error adding T_EXT to global constant memory.\n");

    error_id = cudaMemcpyToSymbol(&X_EXT_GPU, &X_EXT, sizeof(int), 0, cudaMemcpyHostToDevice);
    error(error_id != cudaSuccess, 1, "init_neighbors_gpu", "Error adding X_EXT to global constant memory.\n");

    error_id = cudaMemcpyToSymbol(&Y_EXT_GPU, &Y_EXT, sizeof(int), 0, cudaMemcpyHostToDevice);
    error(error_id != cudaSuccess, 1, "init_neighbors_gpu", "Error adding Y_EXT to global constant memory.\n");

    error_id = cudaMemcpyToSymbol(&Z_EXT_GPU, &Z_EXT, sizeof(int), 0, cudaMemcpyHostToDevice);
    error(error_id != cudaSuccess, 1, "init_neighbors_gpu", "Error adding Z_EXT to global constant memory.\n");

#ifdef FERMION_THETA
    error_id = cudaMemcpyToSymbol(&eitheta_gpu[0], &eitheta[0], sizeof(hr_complex) * 4, 0, cudaMemcpyHostToDevice);
    error(error_id != cudaSuccess, 1, "init_neighbors_gpu", "Error adding Z_EXT to global constant memory.\n");
#endif

    box_t *L = geometryBoxes;
    int number_of_boxes = 0;
    do {
        number_of_boxes++;
    } while ((L = L->next));

    cudaMalloc((void **)&geometryBoxes_gpu, number_of_boxes * sizeof(box_t));

    L = geometryBoxes;
    int i = 0;
    do {
        cudaMemcpy(&geometryBoxes_gpu[i], L, sizeof(box_t), cudaMemcpyHostToDevice);
        ++i;
    } while ((L = L->next));

    // adapted from new_geom.c
    size_t main_mem_volume, buf_mem_volume;
    geometryMemSize(geometryBoxes, &main_mem_volume, &buf_mem_volume);
    cudaMalloc((void **)&icoord_gpu, (main_mem_volume + buf_mem_volume) * sizeof(coord4));
    cudaMemcpy(icoord_gpu, geometryBoxes->icoord, (main_mem_volume + buf_mem_volume) * sizeof(coord4), cudaMemcpyHostToDevice);
    sb_icoord_gpu = icoord_gpu + main_mem_volume;

    L = geometryBoxes;
    box_t *SB;
    do {
        L->icoord = icoord_gpu;
        if ((SB = L->sendBox)) { SB->icoord = sb_icoord_gpu; }
    } while ((L = L->next));

#endif
}
