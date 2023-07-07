/***************************************************************************\
* Copyright (c) 2008, 2022, Claudio Pica, Sofie Martins                     *   
* All rights reserved.                                                      * 
\***************************************************************************/

/**
 * @file
 * @brief The elementary site structures defined in suN.h are used in this 
 *        file to define field structures that contain arrays to save field
 *        data and information on the geometry
 */

#ifndef SPINOR_FIELD_H
#define SPINOR_FIELD_H

#include "Geometry/geometry_descriptor.h"
#include "suN_types.h"
#ifdef WITH_MPI
#include "hr_mpi.h"
#endif

/* Communications */
typedef enum comm_t { CPU_COMM = 1 << 0, GPU_COMM = 1 << 1 } comm_t;
#define ALL_COMMS (CPU_COMM | GPU_COMM)

/* Memory */
typedef enum mem_t { CPU_MEM = 1 << 0, GPU_MEM = 1 << 1 } mem_t;

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief This macro corresponds to a field declaration line to be used 
 *        inside a struct declaration. It only adds a GPU field data
 *        array, if the code is compiled with the flag WITH_GPU.
 *
 * @param _type Elementary site/link type of the field data
 */
#define _GPU_FIELD_DATA(_type)
#ifdef WITH_GPU
#undef _GPU_FIELD_DATA
#define _GPU_FIELD_DATA(_type) _type *gpu_ptr;
#endif

/**
 * @brief This macro corresponds to a field declaration line to be used
 * 	  inside a struct declaration. It only add MPI communication 
 * 	  handles if the code is compiled with the flag WITH_MPI.
 */
#define _MPI_FIELD_DATA(_type)
#ifdef WITH_MPI
#undef _MPI_FIELD_DATA
#define _MPI_FIELD_DATA(_type) \
    MPI_Request *comm_req;     \
    _type *sendbuf_ptr;        \
    _type *sendbuf_gpu_ptr;    \
    _type *recvbuf_gpu_ptr;

//
#endif

/**
 * @brief This macro declares a field struct that contains all necessary 
 *        field data arrays, geometry information and communication handles.
 *
 * @param _name The name of the field type
 * @param _type The elementary type struct of data stored at each site/link.
 */
#define _DECLARE_FIELD_STRUCT(_name, _type) \
    typedef struct _name {                  \
        _type *ptr;                         \
        geometry_descriptor *type;          \
        _MPI_FIELD_DATA(_type)              \
        _GPU_FIELD_DATA(_type)              \
        comm_t comm_type;                   \
        mem_t alloc_t;                      \
    } _name

/**
 * @struct suNg_field
 * @brief Gauge field of SU(N_g) matrices
 *
 * @var suNg_field::ptr	
 * 				Array of SU(N_g) matrices, that contains the field data 
 * 				for the local lattice (in host memory)
 * @var suNg_field::type	
 * 				geometry_descriptor contains geometry information for 
 * 				example pertaining the block decomposition
 * @var suNg_field::gpu_ptr	
 * 				Only active for WITH_GPU. Contains field data on GPU 
 * 				memory.
 * @var suNg_field::comm_req	
 * 				Only active for WITH_MPI. Communication request 
 * 				handles for inter-node communication.
 */
_DECLARE_FIELD_STRUCT(suNg_field, suNg);

/**
 * @struct suNg_scalar_field
 * @brief SU(N_g) scalar field of SU(N_g) vectors 
 *
 * @var suNg_scalar_field::ptr			
 * 				Array of SU(N_g) vectors, that contains the field data for 
 * 				the local lattice (in host memory)
 * @var suNg_scalar_field::type			
 * 				geometry_descriptor contains geometry information, for 
 * 				example pertaining block decomposition
 * @var suNg_scalar_field::gpu_ptr			
 * 				Only active for WITH_GPU. Contains field data in GPU 
 * 				memory.
 * @var comm_req		
 * 				Only active for WITH_MPI. Communication request handles for
 * 				inter-node communication.
 */
_DECLARE_FIELD_STRUCT(suNg_scalar_field, suNg_vector);

/**
 * @struct suNg_field_flt
 * @brief Gauge field of single precision SU(N_g) matrices
 *
 * @var suNg_field_flt::ptr			
 * 				Array of single precision SU(N_g) vectors, that contains
 * 				the field data for the local lattice (in host memory)
 * @var suNg_field_flt::type			
 * 				geometry_descriptor that contains geometry information, for
 * 				example pertaining the block decomposition
 * @var suNg_field_flt::gpu_ptr			
 * 				Only active for WITH_GPU. Contains field data in GPU memory.
 * @var suNg_field_flt::comm_req		
 * 				Only active for WITH_MPI. Communication request handles for
 * 				inter-node communication. 
 */
_DECLARE_FIELD_STRUCT(suNg_field_flt, suNg_flt);

/**
 * @struct suNf_field
 * @brief Gauge field in chosen fermion representation
 *
 * @var suNf_field::ptr			
 * 				Array of SU(N_f) matrices, that contains the field data
 * 				for the local lattice (in host memory)
 * @var suNf_field::type			
 * 				geometry_descriptor that contains geometry information, for
 * 				example pertaining the block decomposition
 * @var suNf_field::gpu_ptr			
 * 				Only active for WITH_GPU. Contains field data in GPU memory.
 * @var suNf_field::comm_req		
 * 				Only active for WITH_MPI. Communication request handles for
 * 				inter-node communication.
 */
_DECLARE_FIELD_STRUCT(suNf_field, suNf);

/**
 * @struct suNf_field_flt
 * @brief Single precision gauge field in the chosen fermion representation
 *
 * @var suNf_field_flt::ptr			
 * 				Array of single precision SU(N_f) matrices, that contains
 * 				the field data for the local lattice (in host memory)
 * @var suNf_field_flt::type			
 * 				geometry_descriptor that contains geometry information, 
 * 				for example pertaining the block decomposition
 * @var suNf_field_flt::gpu_ptr			
 * 				Only active for WITH_GPU. Contains field data in GPU memory
 * @var suNf_field_flt::comm_req		
 * 				Only active for WITH_MPI. Communication request handles for 
 * 				inter-node communication.
 */
_DECLARE_FIELD_STRUCT(suNf_field_flt, suNf_flt);

/**
 * @struct spinor_field
 * @brief Spinor field array containing SU(N_f) spinors in chosen fermion representation
 *
 * @var spinor_field::ptr			
 * 				Array of SU(N_f) spinors, that contains the field data for
 * 				the local lattice (in host memory)
 * @var spinor_field::type			
 * 				geometry_descriptor that contains geometry information,
 * 				for example pertaining the block decomposition
 * @var spinor_field::gpu_ptr			
 * 				Only active for WITH_GPU. Contains field data in GPU memory
 * @var spinor_field::comm_req		
 * 				Only active for WITH_MPI. Communication request handles for
 * 				inter-node communication.
 */
_DECLARE_FIELD_STRUCT(spinor_field, suNf_spinor);

/**
 * @struct spinor_field_flt
 * @brief Spinor field array containing single precision SU(N_f) spinors in chosen
 * 	  fermion representation
 *
 * @var spinor_field_flt::ptr			
 * 				Array of single precision SU(N_f) spinors that contains the
 * 				field data for the local lattice (in host memory)
 * @var spinor_field_flt::type			
 * 				geometry_descriptor that contains geometry information, 
 * 				for example pertaining the block decomposition
 * @var spinor_field_flt::gpu_ptr			
 * 				Only active for WITH_GPU. Contains field data in GPU memory.
 * @var spinor_field_flt::comm_req		
 * 				Only active for WITH_MPI
 */
_DECLARE_FIELD_STRUCT(spinor_field_flt, suNf_spinor_flt);

/**
 * @struct suNg_av_field
 * @brief Field of SU(N_g) algebra vectors
 * @var suNg_av_field::ptr			
 * 				Array of single precision SU(N_g) algebra vectors, that
 * 				contains the field data for the local lattice (in host 
 * 				memory).
 * @var suNg_av_field::type			
 * 				geometry_descriptor that contains geometry information, 
 * 				for example pertaining the block decomposiiton
 * @var suNg_av_field::gpu_ptr			
 * 				Only active for WITH_GPU. Contains field data in GPU memory
 * @var suNg_av_field::comm_req		
 * 				Only active for WITH_MPI
 */
_DECLARE_FIELD_STRUCT(suNg_av_field, suNg_algebra_vector);

/**
 * @struct scalar_field
 * @brief Scalar field of double precision real values
 *
 * @var scalar_field::ptr	Array of double precision reals, that contains the field
 * 				data for the local lattice (in host memory)
 * @var scalar_field::type	
 * 				geometry_descriptor that contains geometry information, 
 * 				for example pertaining the block decomposition
 * @var scalar_field::gpu_ptr			
 * 				Only active for WITH_GPU. Contains field daa in GPU memory
 * @var scalar_field::comm_req		
 * 				Only active for WITH_MPI
 */
_DECLARE_FIELD_STRUCT(scalar_field, double);

/**
 * @struct ldl_field
 * @brief FIXME: Add docs
 */
_DECLARE_FIELD_STRUCT(ldl_field, ldl_t);

/**
 * @struct suNfc_field
 * @brief FIXME: Add docs
 */
_DECLARE_FIELD_STRUCT(suNfc_field, suNfc);

#define _FIELD_NAME spinor_field_f
#define _FIELD_TYPE spinor_field
#define _FIELD_DIM 1
#include "TMPL/spinor_field.h.tmpl"

#define _FIELD_NAME spinor_field_f_flt
#define _FIELD_TYPE spinor_field_flt
#define _FIELD_DIM 1
#include "TMPL/spinor_field.h.tmpl"

#define _FIELD_NAME sfield
#define _FIELD_TYPE scalar_field
#define _FIELD_DIM 1
#include "TMPL/spinor_field.h.tmpl"

#define _FIELD_NAME gfield
#define _FIELD_TYPE suNg_field
#define _FIELD_DIM 4
#include "TMPL/spinor_field.h.tmpl"

#define _FIELD_NAME gfield_flt
#define _FIELD_TYPE suNg_field_flt
#define _FIELD_DIM 4
#include "TMPL/spinor_field.h.tmpl"

#define _FIELD_NAME gfield_f
#define _FIELD_TYPE suNf_field
#define _FIELD_DIM 4
#include "TMPL/spinor_field.h.tmpl"

#define _FIELD_NAME gfield_f_flt
#define _FIELD_TYPE suNf_field_flt
#define _FIELD_DIM 4
#include "TMPL/spinor_field.h.tmpl"

#define _FIELD_NAME suNg_scalar_field
#define _FIELD_TYPE suNg_scalar_field
#define _FIELD_DIM 1
#include "TMPL/spinor_field.h.tmpl"

#define _FIELD_NAME avfield
#define _FIELD_TYPE suNg_av_field
#define _FIELD_DIM 4
#include "TMPL/spinor_field.h.tmpl"

#define _FIELD_NAME gtransf
#define _FIELD_TYPE suNg_field
#define _FIELD_DIM 1
#include "TMPL/spinor_field.h.tmpl"

#define _FIELD_NAME clover_ldl
#define _FIELD_TYPE ldl_field
#define _FIELD_DIM 1
#include "TMPL/spinor_field.h.tmpl"

#define _FIELD_NAME clover_term
#define _FIELD_TYPE suNfc_field
#define _FIELD_DIM 4
#include "TMPL/spinor_field.h.tmpl"

#define _FIELD_NAME clover_force
#define _FIELD_TYPE suNf_field
#define _FIELD_DIM 6
#include "TMPL/spinor_field.h.tmpl"

#define _FIELD_NAME staple_field
#define _FIELD_TYPE suNg_field
#define _FIELD_DIM 3
#include "TMPL/spinor_field.h.tmpl"

#ifdef __cplusplus
}
#endif
#endif
