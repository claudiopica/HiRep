/**
 * @file 
 * @brief Functions for MPI reductions
 */

#ifndef REDUCTION_H
#define REDUCTION_H

#include "hr_complex.h"

// This needs to be there due to the cross-compilation that 
// we need to do if we are compiling WITH_GPU.
#ifdef __cplusplus
   extern "C" {
#endif



#if defined(__clang__)
    #define _ACCURATE_MATH optnone
#elif defined(__GNUC__)
    #define _ACCURATE_MATH optimize("-fno-fast-math")
#elif defined(__FAST_MATH__)
#warning Compensated summation is unsafe with -ffast-math
#endif

/**
 * @brief Perform Kahan summation
 *
 * @param sum		Pointer that contains sum
 * @param c			Pointer that contains the sum compensation
 * @param y			Number to add
 */
static inline void __attribute__((always_inline,_ACCURATE_MATH)) kadd(double *sum, double *c, double y) {
  y -= (*c);
  double t = (*sum) + y;
  (*c) = (t - (*sum)) - y;
  (*sum) = t;
}

/**
 * @brief Collects sum results from the local lattices and sums over all nodes (double).
 *
 * @param d			Pointer that contains the local result
 * @param n			Number of nodes
 */
void global_sum(double *d, int n);

/**
 * @brief Collects sum results from the local lattices and sums over all nodes (integer).
 *
 * @param d			Pointer that contains the local result
 * @param n			Number of nodes
 */
void global_sum_int(int *d, int n);

/**
 * @brief Finds maximum across nodes after finding the local maximum.
 *
 * @param d			Pointer that contains the local result
 * @param n			Number of nodes
 */
void global_max(double *d, int n);

/**
 * @brief Finds minimum across nodes after finding the local minimum.
 *
 * @param d			Pointer that contains the local result
 * @param n			Number of nodes
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

#ifdef __cplusplus
   }
#endif

#ifdef WITH_GPU
   // TODO: Here the cross compilation does not work (SAM)
   #ifdef __cplusplus
      /**
      * @brief Sums across GPU nodes after finding the local sum (generics)
      *
      * @param vector		Vector with local results
      * @param size		Size of vector
      *
      * @return T		Result of sum of generic type T.
      */
      template <class T>
      T global_sum_gpu(T *vector, int size);
      extern "C" {
   #endif

   /**
   * @brief Sums across GPU nodes after finding the local sum (integer)
   *
   * @param vector		Vector with local results
   * @param size		Size of vector
   *
   * @return int 		Result of sum
   */
   int global_sum_gpu_int(int *vector, int size);

   /**
   * @brief Sums across GPU nodes after finding the local sum (single precision reals)
   *
   * @param vector		Vector with local results
   * @param size		Size of vector
   *
   * @return float		Result of sum
   */
   float global_sum_gpu_float(float *vector, int size);

   /**
   * @brief Sums across GPU nodes after finding the local sum (double precision reals)
   *
   * @param vector 		Vector with local results
   * @param size		Size of vector
   *
   * @return double		Result of sum
   */
   double global_sum_gpu_double(double *vector, int size);

   /**
   * @brief Sums across GPU nodes after finding the local sum (single precision complex)
   *
   * @param vector 		Vector with local results
   * @param size		Size of vector
   *
   * @return hr_complex_flt	Result of sum
   */
   hr_complex_flt global_sum_gpu_complex_flt(hr_complex_flt *vector, int size);

   /**
   * @brief Sums across GPU nodes after finding the local sum (double precision complex)
   *
   * @param vector		Vector with local results
   * @param size		Size of vector
   *
   * @return hr_complex	Result of sum
   */
   hr_complex global_sum_gpu_complex(hr_complex *vector, int size);
    
   #ifdef __cplusplus
      }
   #endif

#endif
#endif 
