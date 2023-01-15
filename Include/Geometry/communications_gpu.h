/***************************************************************************\
* Copyright (c) 2022, Claudio Pica, Sofie Martins                           *   
* All rights reserved.                                                      * 
\***************************************************************************/

/**
 * @file communications_gpu.h
 * @brief Communications to send and receive sites from other nodes using MPI with GPU
 */

#ifndef COMMUNICATIONS_GPU_H
#define COMMUNICATIONS_GPU_H

#ifdef WITH_GPU

#include "spinor_field.h"

#ifdef __cplusplus
    extern "C" {
#endif

#define _DECLARE_COMMS(_name, _field_type, _human_readable) \
/** \
    /* @brief Wait for communications between GPUs to finish before continuing. */ \
    /**/ \
    /* @param ##_field_type		##_human_readable that needs to be synchronized */ \
    /*                          across nodes */ \
    /**/ \
void complete_sendrecv_gpu_##_name(_field_type*); \
/** \
    /* @brief Fill buffers and start MPI requests to send and receive. */\
    /**/ \
    /* @param ##_field_type		##_human_readable that needs to be synchronized */ \
    /*                          across nodes */ \
    /**/ \
void start_sendrecv_gpu_##_name(_field_type*); \
/** \
    /* @brief Sync field before communications. This can mean different things */ \
    /*    depending on geometry implementation. */ \
    /**/ \
    /* @param ##_field_type      ##_human_readable that needs to be synchronized */ \
    /*                           on the local lattice. */ \
*/ \
void sync_gpu_##_name(_field_type*); \
void fill_buffers_##_name(_field_type*); \
void fill_buffers_with_zeroes_##_name(_field_type*); 

_DECLARE_COMMS(spinor_field_f, spinor_field, "Spinor field");
_DECLARE_COMMS(spinor_field_f_flt, spinor_field_flt, "Single precision spinor field");
_DECLARE_COMMS(sfield, scalar_field, "Scalar field");
_DECLARE_COMMS(gfield, suNg_field, "Gauge field");
_DECLARE_COMMS(gfield_flt, suNg_field_flt, "Single precision gauge field");
_DECLARE_COMMS(gfield_f, suNf_field, "Represented gauge field");
_DECLARE_COMMS(gfield_f_flt, suNf_field_flt, "Represented single precision gauge field");
_DECLARE_COMMS(suNg_scalar_field, suNg_scalar_field, "SU(N_g) scalar field");
_DECLARE_COMMS(avfield, suNg_av_field, "SU(N_g) algebra vector field");
_DECLARE_COMMS(gtransf, suNg_field, "Gauge transformation");
_DECLARE_COMMS(clover_ldl, ldl_field, "Clover ldl field");
_DECLARE_COMMS(clover_term, suNfc_field, "Clover term");
_DECLARE_COMMS(clover_force, suNf_field, "Clover force");

#undef _DECLARE_COMMS

#ifdef __cplusplus
}
#endif
#endif 
#endif 
