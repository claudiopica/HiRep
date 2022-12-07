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
#include "reduction.h"

#ifdef __cplusplus
    extern "C" {
#endif

#include "spinor_field.h"

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

// TODO: What is this? (SAM)
//void test_spinor_field(spinor_field *p);

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

#ifdef WITH_GPU

    #define _DECLARE_COMMS(_name, _field_type, _human_readable) \
	/** \
	  @brief Wait for communications between GPUs to finish before continuing. \
	  \
	  @param ##_field_type		##_human_readable that needs to be synchronized \
	  				            across nodes \
	 */ \
	void complete_sendrecv_gpu_##_name(_field_type*); \
	/** \
	  @brief Fill buffers and start MPI requests to send and receive.\
	 \
	  @param ##_field_type		##_human_readable that needs to be synchonized \
	  				            across nodes \
	 */ \
	void start_sendrecv_gpu_##_name(_field_type*);

    _DECLARE_COMMS(spinor_field_f, spinor_field, "Spinor field");
    _DECLARE_COMMS(spinor_field_f_flt, spinor_field_flt, "Single precision spinor field");
    _DECLARE_COMMS(sfield, scalar_field, "Scalar field");
    _DECLARE_COMMS(gfield, suNg_field, "Gauge field");
    _DECLARE_COMMS(gfield_flt, suNg_field_flt, "Single precision gauge field");
    _DECLARE_COMMS(gfield_f, suNf_field, "Represented gauge field");
    _DECLARE_COMMS(gfield_f_flt, suNf_field_flt, "Represented single precision gauge field");
    _DECLARE_COMMS(scalar_field, suNg_scalar_field, "SU(N_g) scalar field");
    _DECLARE_COMMS(avfield, suNg_av_field, "SU(N_g) algebra vector field");
    _DECLARE_COMMS(gtransf, suNg_field, "Gauge transformation");
    _DECLARE_COMMS(clover_term, suNfc_field, "Clover term");
    _DECLARE_COMMS(clover_force, suNf_field, "Clover force");

    #undef _DECLARE_COMMS

#endif 


#ifdef __cplusplus
}
#endif
#endif 
