/***************************************************************************\
 * Copyright (c) 2012, Ari Hietanen                                         *   
 * All rights reserved.                                                     * 
\***************************************************************************/

#ifndef ALLOC_TMP_FIELDS_GPU_C
#define ALLOC_TMP_FIELDS_GPU_C
#include "hr_complex.h"
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

hr_complex* alloc_complex_sum_field(int n){
  static hr_complex* res = NULL;
  static int n_size = 0;
  if (n>n_size && res!=NULL){
    cudaFree(res);
    res = NULL;
  }
  if (res == NULL){
    cudaMallocManaged((void **) & res,n*sizeof(hr_complex),cudaMemAttachGlobal);
  }
  return res;
}

#endif
