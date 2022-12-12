#include "global.h"

extern int *iup_gpu, *idn_gpu;

void init_neighbors_gpu() 
{
#ifdef WITH_GPU
  #ifdef WITH_MPI
    int N = T_EXT*X_EXT*Y_EXT*Z_EXT;
  #else
    int N = T*X*Y*Z;
  #endif

  cudaError_t error_id;
  error_id = cudaMalloc((void **)&iup_gpu, 4 * N * sizeof(int));
  error(error_id != cudaSuccess, 1, "init_neighbors_gpu", "Error allocating iup_gpu neighbors array.\n");

  error_id = cudaMalloc((void **)&idn_gpu, 4 * N * sizeof(int));
  error(error_id != cudaSuccess, 1, "init_neighbors_gpu", "Error allocating idn_gpu neighbors array.\n");

  error_id = cudaMalloc((void **)&imask_gpu, N * sizeof(*imask));
  error(error_id != cudaSuccess, 1, "init_neighbors_gpu", "Error allocating imask_gpu lookup table.\n");

  error_id = cudaMalloc((void **)&ipt_gpu, (X+2*X_BORDER)*(Y+2*Y_BORDER)*(Z+2*Z_BORDER)*(T+2*T_BORDER));
  error(error_id != cudaSuccess, 1, "init_neighbors_gpu", "Error allocating ipt_gpu lookup table.\n");

  error_id = cudaMemcpy(iup_gpu, iup, 4 * N * sizeof(int), cudaMemcpyHostToDevice);
  error(error_id != cudaSuccess, 1, "init_neighbors_gpu", "Error copying iup neighbors array to device memory.\n");

  error_id = cudaMemcpy(idn_gpu, idn, 4 * N * sizeof(int), cudaMemcpyHostToDevice);
  error(error_id != cudaSuccess, 1, "init_neighbors_gpu", "Error copying idn neighbors array to device memory.\n");

  error_id = cudaMemcpy(imask_gpu, imask, N * sizeof(*imask), cudaMemcpyHostToDevice);
  error(error_id != cudaSuccess, 1, "init_neighbors_gpu", "Error copying imask lookup table to device memory.\n");

  error_id = cudaMemcpy(ipt_gpu, ipt, (X+2*X_BORDER)*(Y+2*Y_BORDER)*(Z+2*Z_BORDER)*(T+2*T_BORDER), cudaMemcpyHostToDevice);
  error(error_id != cudaSuccess, 1, "init_neighbors_gpu", "Error copying ipt to device memory.\n");
#endif
}