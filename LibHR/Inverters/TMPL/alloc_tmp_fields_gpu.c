/***************************************************************************\
 * Copyright (c) 2012, Ari Hietanen                                         *   
 * All rights reserved.                                                     * 
\***************************************************************************/

#ifndef ALLOC_TMP_FIELDS_GPU_C
#define ALLOC_TMP_FIELDS_GPU_C

#include <gpu.h>

double* alloc_double_sum_field(int n){
  static double* res = NULL;
  static int n_size=0;
  if (n>n_size && res!=NULL){
    cudaFree(res);
    res = NULL;
  }

  if (res == NULL){
    cudaMallocManaged((void **) & res,n*sizeof(double),cudaMemAttachGlobal);
    n=n_size;
  }
  return res;
}

complex* alloc_complex_sum_field(int n){
  static complex* res = NULL;
  static int n_size = 0;
  if (n>n_size && res!=NULL){
    cudaFree(res);
    res = NULL;
  }
  if (res == NULL){
    cudaMallocManaged((void **) & res,n*sizeof(complex),cudaMemAttachGlobal);
  }
  return res;
}

#endif
