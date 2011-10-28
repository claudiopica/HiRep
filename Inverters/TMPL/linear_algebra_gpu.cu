#ifndef LINEAR_ALGEBRA_GPU_CU
#define LINEAR_ALGEBRA_GPU_CU

#include <cuda.h>
#include "global.h"

template <typename SPINOR_TYPE, REAL>
__global__ void spinor_field_mul_add_assign_gpu(SPINOR_TYPE *s1, REAL r, SPINOR_TYPE *s2,int N){
  int i = blockIdx.x*BLOCK_SIZE + threadId.x;
  _spinor_mul_add_assign_f(s1[i],r,s2[i]);
}

#endif
