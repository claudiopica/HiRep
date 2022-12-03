/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *
* All rights reserved.                                                      *
\***************************************************************************/

/**
 * @file memory.h
 * @brief Memory handling functions
 */

#ifndef MEMORY_H
#define MEMORY_H

#include <stdlib.h>
#include "suN.h"
#include "spinor_field.h"

#ifdef P4
#  define ALIGN 7
#else
#define ALIGN 8
#endif

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Allocated memory aligned, because this improves bandwidth
 *
 * @param size			size to be allocated
 * @param p			alignment
 */
void *amalloc(size_t size,int p);

/**
 * @brief Free memory that was allocated aligned using amalloc
 *
 * @param addr			Free this pointer
 */
void afree(void *addr);

#define _DECLARE_MEMORY_FUNC_GAUGE(_name, _field_type, _site_type, _size, _human_readable) \
	/**
	  @brief Copy _human_readable from the host to the device by synchronizing the GPU
	         field data with the CPU field data

          @param _field_type		Field to be copied.
	 */ \
        void copy_to_gpu_##_name(_field_type*); \
	/**
	  @brief Copy _human_readable from the device to the host by synchronizing the CPU
	         field data with the GPU field data

          @param _field_type		Field to be copied.
	 */ \
        void copy_from_gpu_##_name(_field_type*); \
	/**
	  @brief Convert _human_readable field data host geometry layout to the device 
	         geometry layout. Read more on GPU geometry in the corresponding section
		 in the development manual.

	  @param _field_type		_human_readable that will contain the converted
	  				field data.
	  @param _field_type		_human_readable that contains the initial field data.
	 */ \
        void to_gpu_format_##_name(_field_type*, _field_type*); \
	/**
	  @brief Convert _human_readable field data device geometry layout to host geometry
	         layout. Read more on GPU geometry in the corresponding section in the 
		 development manual. 
	
	  @param _field_type		_human_readable that will contain the converted
	   				field_data
	  @param _field_type		_human_readable that contains the initial field 
	   				data
	 */ \
        void to_cpu_format_##_name(_field_type*, _field_type*); \
	/**
	  @brief Free field data, other struct fields and struct pointer of field struct.

	  @param _field_type		Field to be freed.
	 */ \
        void free_##_name(_field_type*); \
	/**
	  @brief Allocate field struct pointer

	  @param geometry_descriptor	Underlying lattice geometry to allocate on.
	 */
        _field_type *alloc_##_name(geometry_descriptor*);

#define _DECLARE_MEMORY_FUNC_SPINOR(_name, _field_type, _site_type, _size, _human_readable) \
	/**
	  @brief Copy _human_readable from the host to the device by synchronizing the
	         GPU field data with the CPU field data
	 
	  @param _field_type		Field to be copied.
	 */ \
        void copy_to_gpu_##_name(_field_type*); \
	/**
	  @brief Copy _human_readable from the device to the host by synchronizing the CPU
	         field data with the GPU field data
	 
	  @param _field_type		Field to be copied.
	 */ \
        void copy_from_gpu_##_name(_field_type*); \
	/**
	  @brief Convert _human_readable field data host geometry layout to device geometry
	          layout. Read more on GPU geometry in the corresponding section in the 
	          development manual.
	 
          @param _field_type		_human_readable that will contain the converted 
	  				field data.
	  @param _field_type		_human_readable that contains the initial field data.
	 */ \
        void to_gpu_format_##_name(_field_type*, _field_type*); \
	/**
	  @brief Convert _human_readable field data device geometry layout to host 
	  	  geometry layout. Read more on GPU geometry in the corresponding section
	  	  in the development manual.

	  @param _field_type		_human_readable that will contain the converted
	  				field data
	  @param _field_type		_human_readable that contains the initial field
	  				data
	 */ \
        void to_cpu_format_##_name(_field_type*, _field_type*); \
	/**
	  @brief Free field data, other struct fields and struct pointer of field struct.
	 
	  @param _field_type		Field to be freed.
	 */ \
     	void free_##_name(_field_type*); \
	/**
	  @brief Allocate field struct pointer and struct fields, in particular the field
	          data arrays.
	 
	  @param unsigned int		Number of spinors
	  @param geometry_descriptor	Underlying lattice geometry to allocate on.
	 */ \
         _field_type *alloc_##_name(unsigned int, geometry_descriptor*);

// Gauge Fields
_DECLARE_MEMORY_FUNC_GAUGE(gfield, suNg_field, suNg, 4);
_DECLARE_MEMORY_FUNC_GAUGE(gfield_flt, suNg_field_flt, suNg_flt, 4);
_DECLARE_MEMORY_FUNC_GAUGE(gfield_f, suNf_field, suNf, 4);
_DECLARE_MEMORY_FUNC_GAUGE(gfield_f_flt, suNf_field_flt, suNf_flt, 4);
_DECLARE_MEMORY_FUNC_GAUGE(scalar_field, suNg_scalar_field, suNg_vector, 1);
_DECLARE_MEMORY_FUNC_GAUGE(avfield, suNg_av_field, suNg_algebra_vector, 4);
_DECLARE_MEMORY_FUNC_GAUGE(gtransf, suNg_field, suNg, 1);
//_DECLARE_MEMORY_FUNC_GAUGE(clover_ldl, ldl_field, ldl_t, 1);
_DECLARE_MEMORY_FUNC_GAUGE(clover_term, suNfc_field, suNfc, 4);
_DECLARE_MEMORY_FUNC_GAUGE(clover_force, suNf_field, suNf, 6);

// Matter Fields
_DECLARE_MEMORY_FUNC_SPINOR(spinor_field_f, spinor_field, suNf_spinor, 1, Spinor field);
_DECLARE_MEMORY_FUNC_SPINOR(spinor_field_f_flt, spinor_field_flt, suNf_spinor_flt, 1, Spinor field single precision);
_DECLARE_MEMORY_FUNC_SPINOR(sfield, scalar_field, double, 1, Double precision real field);

#undef _DECLARE_MEMORY_FUNC

#ifdef __cplusplus
}
#endif

#endif
