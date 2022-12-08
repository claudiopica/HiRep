#include "global.h"

extern int *iup_gpu, *idn_gpu;

void init_neighbors_gpu() 
{
#ifdef WITH_GPU
  int N = T*X*Y*Z;
  cudaError_t error_id;
  error_id = cudaMalloc((void **)&iup_gpu, 4 * N * sizeof(int));
  error(error_id != cudaSuccess, 1, "init_neighbors_gpu", "Error allocating iup_gpu neighbors array.\n");

  error_id = cudaMalloc((void **)&idn_gpu, 4 * N * sizeof(int));
  error(error_id != cudaSuccess, 1, "init_neighbors_gpu", "Error allocating idn_gpu neighbors array.\n");

  error_id = cudaMemcpy(iup_gpu, iup, 4 * N * sizeof(int), cudaMemcpyHostToDevice);
  error(error_id != cudaSuccess, 1, "init_neighbors_gpu", "Error copying iup neighbors array to device memory.\n");

  error_id = cudaMemcpy(idn_gpu, idn, 4 * N * sizeof(int), cudaMemcpyHostToDevice);
  error(error_id != cudaSuccess, 1, "init_neighbors_gpu", "Error copying idn neighbors array to device memory.\n");
#endif
}