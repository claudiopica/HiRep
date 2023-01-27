/***************************************************************************\
* Copyright (c) 2008, 2022, Claudio Pica, Sofie Martins                     *   
* All rights reserved.                                                      * 
\***************************************************************************/

/**
 * @file communications.h
 * @brief Communications to send and receive sites from other nodes using MPI
 */

#ifndef COMMUNICATIONS_H
#define COMMUNICATIONS_H

#include "spinor_field.h"

#ifdef __cplusplus
    extern "C" {
#endif

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

#ifdef __cplusplus
    }
#endif
#endif 
