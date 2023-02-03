/***************************************************************************\
* Copyright (c) 2008, 2022, Claudio Pica, Sofie Martins                     *   
* All rights reserved.                                                      * 
\***************************************************************************/

/**
 * @file global_sum.h
 * @brief Functions for global reduction operations across nodes
 */

#ifndef GLOBAL_SUM_H
#define GLOBAL_SUM_H

#include "hr_complex.h"
#include "spinor_field.h"

/**
 * @brief Collects sum results from the local lattices and sums over all nodes (double).
 *
 * @param d			Pointer that contains the local result
 * @param n			Size of the array
 */
void global_sum(double *d, int n);

/**
 * @brief Collects sum results from the local lattices and sums over all nodes (integer).
 *
 * @param d			Pointer that contains the local result
 * @param n			Size of the array
 */
void global_sum_int(int *d, int n);

/**
 * @brief Finds maximum across nodes after finding the local maximum.
 *
 * @param d			Pointer that contains the local result
 * @param n			Size of the array
 */
void global_max(double *d, int n);
void global_max_flt(float *d, int n);

/**
 * @brief Finds minimum across nodes after finding the local minimum.
 *
 * @param d			Pointer that contains the local result
 * @param n			Size of the array
 */
void global_min(double *d, int n);

/**
 * @brief FIXME: add docs
 *
 * @param d
 * @param n
 */
void bcast(double *d, int n);

/**
 * @brief FIXME: add docs
 *
 * @param i
 * @param n
 */
void bcast_int(int *i, int n);

#ifdef WITH_GPU

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
