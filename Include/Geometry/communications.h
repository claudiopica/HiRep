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

/**
 * @brief Wait for communications to finish before continuing.
 *
 * @param suNg_field		Gauge field that needs to be synchronized across nodes.
 */
//extern void (*complete_sendrecv_gfield) (suNg_field*);

/**
 * @brief Load buffers and start MPI requests to send and receive.
 *
 * @param suNg_field		Gauge field that needs to be synchronized across nodes.
 */
//extern void (*start_sendrecv_gfield) (suNg_field*);

/**
 * @brief Wait for communications to finish before continuing.
 *
 * @param spinor_field		Spinor field that needs to be synchronized across nodes.
 */
//extern void (*complete_sendrecv_spinor_field_f) (spinor_field*);

/**
 * @brief Load buffers and start MPI requests to send and receive.
 *
 * @param spinor_field 		Spinor field that needs to be synchronized across nodes.
 */
//extern void (*start_sendrecv_spinor_field_f) (spinor_field*);

/**
 * @brief Wait for communications to finish before continuing.
 *
 * @param suNg_scalar_field	Scalar field that needs to be synchronized across nodes.
 */
//extern void (*complete_sendrecv_suNg_scalar_field) (suNg_scalar_field*);

/**
 * @brief Load buffers and start MPI requests to send and receive.
 *
 * @param suNg_scalar_field 	Scalar field that needs to be synchronized across nodes.
 */
//extern void (*start_sendrecv_suNg_scalar_field) (suNg_scalar_field*);

/**
 * @brief Wait for communications to finish before continuing.
 *
 * @param suNf_field 		Represented gauge field that needs to be synchronized 
 * 				across nodes.
 */
//extern void (*complete_clover_force_sendrecv) (suNf_field*);

/**
 * @brief Load buffers and start MPI requests to send and receive.
 *
 * @param suNf_field		Represented gauge field that needs to be synchronized 
 * 				across nodes.
 */
//extern void (*start_clover_force_sendrecv) (suNf_field*);

/**
 * @brief Wait for communications to finish before continuing.
 *
 * @param suNg_field		Gauge transformation that needs to be synchronized
 * 				across nodes.
 */
//extern void (*complete_sendrecv_gtransf) (suNg_field*);

/**
 * @brief Fill buffers and start MPI requests to send and receive.
 *
 * @param 			Gauge transformation that needs to be synchronized across
 * 				nodes
 */
//extern void (*start_sendrecv_gtransf) (suNg_field*);

/**
 * @brief Wait for communications to finish before continuing.
 *
 * @param suNg_field_flt	Single precision gauge field that needs to be synchronized
 * 				across nodes.
 */
//extern void (*complete_sendrecv_gfield_flt) (suNg_field_flt*);

/**
 * @brief Fill buffers and start MPI requests to send and receive.
 *
 * @param suNg_field_flt	Single precision gauge field that needs to be synchronized
 * 				across nodes.
 */
//extern void (*start_sendrecv_gfield_flt) (suNg_field_flt*);

/**
 * @brief Wait for communications to finish before continuing.
 *
 * @param spinor_field_flt	Single precision spinor field that needs to be synchronized
 * 				across nodes.
 */
//extern void (*complete_sendrecv_spinor_field_f_flt) (spinor_field_flt*);

/**
 * @brief Fill buffers and start MPI requests to send and receive.
 *
 * @param spinor_field_flt	Single precision spinor field that needs to be synchronized
 * 				across nodes.
 */
//extern void (*start_sendrecv_spinor_field_f_flt) (spinor_field_flt*);

/**
 * @brief Wait for communications to finish before continuing.
 *
 * @param staple field that needs to be synchronized across
 * 				nodes
 */
//extern void (*complete_sendrecv_staple_field) (suNg_field *);

/**
 * @brief Fill buffers and start MPI requests to send and receive.
 *
 * @param staple field that needs to be synchronized across
 * 				nodes
 */
//extern void (*start_sendrecv_staple_field) (suNg_field*);

void complete_sendrecv_gfield_cpu(suNg_field*);
void start_sendrecv_gfield_cpu(suNg_field*);
void complete_sendrecv_spinor_field_f_cpu(spinor_field*);
void start_sendrecv_spinor_field_f_cpu(spinor_field*);
void complete_sendrecv_suNg_scalar_field_cpu(suNg_scalar_field*);
void start_sendrecv_suNg_scalar_field_cpu(suNg_scalar_field*);
void complete_clover_force_sendrecv_cpu(suNf_field*);
void start_clover_force_sendrecv_cpu(suNf_field*);
void complete_sendrecv_gtransf_cpu(suNg_field*);
void start_sendrecv_gtransf_cpu(suNg_field*);
void complete_sendrecv_gfield_flt_cpu(suNg_field_flt*);
void start_sendrecv_gfield_flt_cpu(suNg_field_flt*);
void complete_sendrecv_spinor_field_f_flt_cpu(spinor_field_flt*);
void start_sendrecv_spinor_field_f_flt_cpu(spinor_field_flt*);
void complete_sendrecv_staple_field_cpu(suNg_field*);
void start_sendrecv_staple_field_cpu(suNg_field*);

#ifdef __cplusplus
    }
#endif
#endif 
