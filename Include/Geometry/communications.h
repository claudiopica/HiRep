/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *   
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

/**
 * @brief Wait for communications to finish before continuing.
 *
 * @param suNg_field		Gauge field that needs to be synchronized across nodes.
 */
void complete_gf_sendrecv(suNg_field*);

/**
 * @brief Load buffers and start MPI requests to send and receive.
 *
 * @param suNg_field		Gauge field that needs to be synchronized across nodes.
 */
void start_gf_sendrecv(suNg_field*);

/**
 * @brief Wait for communications to finish before continuing.
 *
 * @param spinor_field		Spinor field that needs to be synchronized across nodes.
 */
void complete_sf_sendrecv(spinor_field*);

/**
 * @brief Load buffers and start MPI requests to send and receive.
 *
 * @param spinor_field 		Spinor field that needs to be synchronized across nodes.
 */
void start_sf_sendrecv(spinor_field*);

/**
 * @brief Wait for communications to finish before continuing.
 *
 * @param suNg_scalar_field	Scalar field that needs to be synchronized across nodes.
 */
void complete_sc_sendrecv(suNg_scalar_field*);

/**
 * @brief Load buffers and start MPI requests to send and receive.
 *
 * @param suNg_scalar_field 	Scalar field that needs to be synchronized across nodes.
 */
void start_sc_sendrecv(suNg_scalar_field*);

/**
 * @brief Wait for communications to finish before continuing.
 *
 * @param suNf_field 		Represented gauge field that needs to be synchronized 
 * 				across nodes.
 */
void complete_clover_force_sendrecv(suNf_field*);

/**
 * @brief Load buffers and start MPI requests to send and receive.
 *
 * @param suNf_field		Represented gauge field that needs to be synchronized 
 * 				across nodes.
 */
void start_clover_force_sendrecv(suNf_field*);

/**
 * @brief Wait for communications to finish before continuing.
 *
 * @param suNg_field		Gauge transformation that needs to be synchronized
 * 				across nodes.
 */
void complete_gt_sendrecv(suNg_field*);

/**
 * @brief Fill buffers and start MPI requests to send and receive.
 *
 * @param 			Gauge transformation that needs to be synchronized across
 * 				nodes
 */
void start_gt_sendrecv(suNg_field*);

/**
 * @brief Wait for communications to finish before continuing.
 *
 * @param suNg_field_flt	Single precision gauge field that needs to be synchronized
 * 				across nodes.
 */
void complete_gf_sendrecv_flt(suNg_field_flt*);

/**
 * @brief Fill buffers and start MPI requests to send and receive.
 *
 * @param suNg_field_flt	Single precision gauge field that needs to be synchronized
 * 				across nodes.
 */
void start_gf_sendrecv_flt(suNg_field_flt*);

/**
 * @brief Wait for communications to finish before continuing.
 *
 * @param spinor_field_flt	Single precision spinor field that needs to be synchronized
 * 				across nodes.
 */
void complete_sf_sendrecv_flt(spinor_field_flt*);

/**
 * @brief Fill buffers and start MPI requests to send and receive.
 *
 * @param spinor_field_flt	Single precision spinor field that needs to be synchronized
 * 				across nodes.
 */
void start_sf_sendrecv_flt(spinor_field_flt*);


void complete_staple_field_sendrecv(suNg_field *);

/**
 * @brief Wait for communications to finish before continuing.
 *
 * @param staple field that needs to be synchronized across
 * 				nodes
 */
void start_staple_field_sendrecv(suNg_field*);

/**
 * @brief Fill buffers and start MPI requests to send and receive.
 *
 * @param staple field that needs to be synchronized across
 * 				nodes
 */

#ifdef __cplusplus
    }
#endif
#endif 
