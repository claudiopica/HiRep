/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *   
* All rights reserved.                                                      * 
\***************************************************************************/

/**************************************************************************\
 *
 * File spinor_field.h
 *
 * Definitions of field structures
 *
 **************************************************************************/

/**
 * @file
 * @brief The elementary site structures defined in suN.h are used in this 
 *        file to define field structures that contain arrays to save field
 *        data and information on the geometry
 */

#ifndef SPINOR_FIELD_H
#define SPINOR_FIELD_H

#include "geometry.h"
#include "suN_types.h"
#include "error.h"
#include "hr_complex.h"
#ifdef WITH_MPI
#include <mpi.h>
#endif

/* GPU data */
#define _GPU_FIELD_DATA(_type)
#ifdef WITH_GPU
#undef _GPU_FIELD_DATA
#define _GPU_FIELD_DATA(_type) _type *gpu_ptr;
#endif //WITH_GPU

/* MPI data */
#define _MPI_FIELD_DATA(_type)
#ifdef WITH_MPI
#undef _MPI_FIELD_DATA
#define _MPI_FIELD_DATA(_type) MPI_Request *comm_req;
#endif //WITH_MPI

typedef struct {// TODO: this is probably not the right complex type
	double _Complex up[NF*(2*NF+1)];
	double _Complex dn[NF*(2*NF+1)];
} ldl_t;

#define _DECLARE_FIELD_STRUCT(_name,_type) \
	typedef struct { \
		_type *ptr; \
		geometry_descriptor *type; \
		_MPI_FIELD_DATA(_type) \
		_GPU_FIELD_DATA(_type) \
	} _name

//typedef struct _##_name { \
//_type *ptr; \
//geometry_descriptor *type;\
//_MPI_FIELD_DATA(_type) \
//_GPU_FIELD_DATA(_type) \
//} _name


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
 * @struct _ldl_field
 * @brief FIXME: Add docs
 */
_DECLARE_FIELD_STRUCT(ldl_field, ldl_t);

/**
 * @struct _suNfc_field
 * @brief FIXME: Add docs
 */
_DECLARE_FIELD_STRUCT(suNfc_field, suNfc);


/* LOOPING MACRO */

#ifdef CHECK_SPINOR_MATCHING

#define _TWO_SPINORS_MATCHING(s1,s2) \
	error((s1)->type!=(s2)->type,1,__FILE__ ": ", "Spinors don't match!");

#define _ARRAY_SPINOR_MATCHING(s,n) \
	for(int _i=0; _i<n; _i++) \
		error((s)->type!=((s)+_i)->type,1,__FILE__ ": ", "Spinors don't match!");

#else /* CHECK_SPINOR_MATCHING */

#define _TWO_SPINORS_MATCHING(s1,s2)

#define _ARRAY_SPINOR_MATCHING(s,n)

#endif /* CHECK_SPINOR_MATCHING */


#define _ONE_SPINOR_FOR_RED(s,redop1,redop2) _MASTER_FOR_RED((s)->type,_spinor_for_is,redop1,redop2)
#define _ONE_SPINOR_FOR(s) _ONE_SPINOR_FOR_RED(s,,)
#define _ONE_SPINOR_FOR_SUM(s,...) _ONE_SPINOR_FOR_RED(s,_omp_sum(__VA_ARGS__),)

#define _TWO_SPINORS_FOR_RED(s1,s2,redop1,redop2) \
  _TWO_SPINORS_MATCHING(s1,s2); \
  _ONE_SPINOR_FOR_RED(s1,redop1,redop2)
#define _TWO_SPINORS_FOR(s1,s2) _TWO_SPINORS_FOR_RED(s1,s2,,)
#define _TWO_SPINORS_FOR_SUM(s1,s2,...) _TWO_SPINORS_FOR_RED(s1,s2,_omp_sum(__VA_ARGS__),)

#define _THREE_SPINORS_FOR_RED(s1,s2,s3,redop1,redop2) \
  _TWO_SPINORS_MATCHING(s1,s2); \
  _TWO_SPINORS_MATCHING(s1,s3); \
  _ONE_SPINOR_FOR_RED(s1,redop1,redop2)
#define _THREE_SPINORS_FOR(s1,s2,s3) _THREE_SPINORS_FOR_RED(s1,s2,s3,,)

#include "field_ordering.h"

#define _FIELD_AT(s,i) (((s)->ptr)+i-(s)->type->master_shift)
#define _4FIELD_AT(s,i,mu) (((s)->ptr)+coord_to_index(i-(s)->type->master_shift,mu))
#define _6FIELD_AT(s,i,mu) (((s)->ptr)+((i-(s)->type->master_shift)*6+mu))
#define _DFIELD_AT(s,i,mu,size) (size == 1) ? _FIELD_AT(s,i) : ((size == 4) ? _4FIELD_AT(s,i,mu) : _6FIELD_AT(s,i,mu))

#define _SPINOR_PTR(s) _FIELD_AT(s,_spinor_for_is)

#ifdef  WITH_GPU
  #define _GPU_FIELD_AT(s,i) (((s)->gpu_ptr)+i-(s)->type->master_shift)
  #define _GPU_4FIELD_AT(s,i,mu) (((s)->gpu_ptr)+coord_to_index(i-(s)->type->master_shift,mu))
  #define _GPU_6FIELD_AT(s,i,mu) (((s)->gpu_ptr)+((i-(s)->type->master_shift)*6+mu))
#endif

#endif
