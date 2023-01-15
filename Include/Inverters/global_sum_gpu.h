/**
 * @file 
 * @brief Functions for gloabl redution operations on MPI and GPU
 */

#ifndef GLOBAL_SUM_GPU_H
#define GLOBAL_SUM_GPU_H

#ifdef WITH_GPU

#include "hr_complex.h"

#ifdef __cplusplus
   extern "C" {
#endif

/**
* @brief Sums across GPU nodes after finding the local sum (integer)
*
* @param vector		Vector with local results
* @param size		   Size of vector
*
* @return int 		   Result of sum
*/
int global_sum_gpu_int(int *vector, int size);

/**
* @brief Sums across GPU nodes after finding the local sum (single precision reals)
*
* @param vector		Vector with local results
* @param size		   Size of vector
*
* @return float		Result of sum
*/
float global_sum_gpu_float(float *vector, int size);

/**
* @brief Sums across GPU nodes after finding the local sum (double precision reals)
*
* @param vector 		Vector with local results
* @param size		   Size of vector
*
* @return double		Result of sum
*/
double global_sum_gpu_double(double *vector, int size);

/**
* @brief Sums across GPU nodes after finding the local sum (single precision complex)
*
* @param vector 		Vector with local results
* @param size		   Size of vector
*
* @return hr_complex_flt	Result of sum
*/
hr_complex_flt global_sum_gpu_complex_flt(hr_complex_flt *vector, int size);

/**
* @brief Sums across GPU nodes after finding the local sum (double precision complex)
*
* @param vector		Vector with local results
* @param size		   Size of vector
*
* @return hr_complex	Result of sum
*/
hr_complex global_sum_gpu_complex(hr_complex *vector, int size);

    
#ifdef __cplusplus
   }
#endif
#endif
#endif 
