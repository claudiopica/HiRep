/***************************************************************************\
* Copyright (c) 2023, Sofie Martins                                         *   
* All rights reserved.                                                      * 
\***************************************************************************/

// Ported functions from boundary_conditions.c by Agostino Patella
// and Claudio Pica

#include "libhr_core.h"
#include "utils.h"
#include "geometry.h"
#include "io.h"

#define ipt_ext_gpu_loc(t, x, y, z) ipt_gpu_d[_lexi(T_EXT_GPU, X_EXT_GPU, Y_EXT_GPU, Z_EXT_GPU, t, x, y, z)]

__global__ void apply_boundary_conditions_T(suNf *g, int border, int *ipt_gpu_d, int xmax, int ymax, int zmax) {
    for (int ix = blockDim.x * blockIdx.x + threadIdx.x; ix < xmax; ix += gridDim.x * blockDim.x) {
        for (int iy = blockDim.y * blockIdx.y + threadIdx.y; iy < ymax; iy += gridDim.y * blockDim.y) {
            for (int iz = blockDim.z * blockIdx.z + threadIdx.z; iz < zmax; iz += gridDim.z * blockDim.z) {
                int index = ipt_ext_gpu_loc(2 * border, ix, iy, iz);
                if (index != -1) {
                    suNf u;
                    read_gpu<double>(0, &u, g, index, 0, 4);
                    _suNf_minus(u, u);
                    write_gpu<double>(0, &u, g, index, 0, 4);
                }
            }
        }
    }
}

__global__ void apply_boundary_conditions_X(suNf *g, int border, int *ipt_gpu_d, int tmax, int ymax, int zmax) {
    for (int it = blockDim.x * blockIdx.x + threadIdx.x; it < tmax; it += gridDim.x * blockDim.x) {
        for (int iy = blockDim.y * blockIdx.y + threadIdx.y; iy < ymax; iy += gridDim.y * blockDim.y) {
            for (int iz = blockDim.z * blockIdx.z + threadIdx.z; iz < zmax; iz += gridDim.z * blockDim.z) {
                int index = ipt_ext_gpu_loc(it, 2 * border, iy, iz);
                if (index != -1) {
                    suNf u;
                    read_gpu<double>(0, &u, g, index, 1, 4);
                    _suNf_minus(u, u);
                    write_gpu<double>(0, &u, g, index, 1, 4);
                }
            }
        }
    }
}

__global__ void apply_boundary_conditions_Y(suNf *g, int border, int *ipt_gpu_d, int tmax, int xmax, int zmax) {
    for (int it = blockDim.x * blockIdx.x + threadIdx.x; it < tmax; it += gridDim.x * blockDim.x) {
        for (int ix = blockDim.y * blockIdx.y + threadIdx.y; ix < xmax; ix += gridDim.y * blockDim.y) {
            for (int iz = blockDim.z * blockIdx.z + threadIdx.z; iz < zmax; iz += gridDim.z * blockDim.z) {
                int index = ipt_ext_gpu_loc(it, ix, 2 * border, iz);
                if (index != -1) {
                    suNf u;
                    read_gpu<double>(0, &u, g, index, 2, 4);
                    _suNf_minus(u, u);
                    write_gpu<double>(0, &u, g, index, 2, 4);
                }
            }
        }
    }
}

__global__ void apply_boundary_conditions_Z(suNf *g, int border, int *ipt_gpu_d, int tmax, int xmax, int ymax) {
    for (int it = blockDim.x * blockIdx.x + threadIdx.x; it < tmax; it += gridDim.x * blockDim.x) {
        for (int ix = blockDim.y * blockIdx.y + threadIdx.y; ix < xmax; ix += gridDim.y * blockDim.y) {
            for (int iy = blockDim.z * blockIdx.z + threadIdx.z; iy < ymax; iy += gridDim.z * blockDim.z) {
                int index = ipt_ext_gpu_loc(it, ix, iy, 2 * border);
                if (index != -1) {
                    suNf u;
                    read_gpu<double>(0, &u, g, index, 3, 4);
                    _suNf_minus(u, u);
                    write_gpu<double>(0, &u, g, index, 3, 4);
                }
            }
        }
    }
}

void sp_T_antiperiodic_BCs_gpu(suNf_field *flt) {
    if (COORD[0] == 0) {
        int block_dim = 4;
        int grid_x = (X_EXT - 1) / block_dim + 1;
        int grid_y = (Y_EXT - 1) / block_dim + 1;
        int grid_z = (Z_EXT - 1) / block_dim + 1;
        dim3 block(block_dim, block_dim, block_dim);
        dim3 grid(grid_x, grid_y, grid_z);
        apply_boundary_conditions_T<<<grid, block>>>(flt->gpu_ptr, T_BORDER, ipt_gpu, X_EXT, Y_EXT, Z_EXT);
        cudaDeviceSynchronize();
    }
}

void sp_X_antiperiodic_BCs_gpu(suNf_field *flt) {
    if (COORD[1] == 0) {
        int block_dim = 4;
        int grid_t = (T_EXT - 1) / block_dim + 1;
        int grid_y = (Y_EXT - 1) / block_dim + 1;
        int grid_z = (Z_EXT - 1) / block_dim + 1;
        dim3 block(block_dim, block_dim, block_dim);
        dim3 grid(grid_t, grid_y, grid_z);
        apply_boundary_conditions_X<<<grid, block>>>(flt->gpu_ptr, X_BORDER, ipt_gpu, T_EXT, Y_EXT, Z_EXT);
        cudaDeviceSynchronize();
    }
}

void sp_Y_antiperiodic_BCs_gpu(suNf_field *flt) {
    if (COORD[2] == 0) {
        int block_dim = 4;
        int grid_t = (T_EXT - 1) / block_dim + 1;
        int grid_x = (X_EXT - 1) / block_dim + 1;
        int grid_z = (Z_EXT - 1) / block_dim + 1;
        dim3 block(block_dim, block_dim, block_dim);
        dim3 grid(grid_t, grid_x, grid_z);
        apply_boundary_conditions_Y<<<grid, block>>>(flt->gpu_ptr, Y_BORDER, ipt_gpu, T_EXT, X_EXT, Z_EXT);
        cudaDeviceSynchronize();
    }
}

void sp_Z_antiperiodic_BCs_gpu(suNf_field *flt) {
    if (COORD[3] == 0) {
        int block_dim = 4;
        int grid_t = (T_EXT - 1) / block_dim + 1;
        int grid_x = (X_EXT - 1) / block_dim + 1;
        int grid_y = (Y_EXT - 1) / block_dim + 1;
        dim3 block(block_dim, block_dim, block_dim);
        dim3 grid(grid_t, grid_x, grid_y);
        apply_boundary_conditions_Z<<<grid, block>>>(flt->gpu_ptr, Z_BORDER, ipt_gpu, T_EXT, X_EXT, Y_EXT);
        cudaDeviceSynchronize();
    }
}