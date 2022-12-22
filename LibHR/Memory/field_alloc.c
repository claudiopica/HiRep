/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *   
* All rights reserved.                                                      * 
\***************************************************************************/

/**
 * @file field_alloc.c
 * @brief Allocation and free functions for all fields defined in spinor_field.h
 */

#include <stdlib.h>
#include "suN.h"
#include "error.h"
#include "memory.h"
#include "global.h"
#include "spinor_field.h"
#include "geometry.h"
#include "geometry_check.h"
#include "gpu_geometry.h"
#include "gpu.h"
#ifdef WITH_MPI
#include <mpi.h>
#endif

#include "alloc_cpu_field_data.c"
#include "alloc_gpu_field_data.c"
#include "alloc_mpi_data.c"

/**
 * @brief Generate free function for all fields defined in spinor_field.h
 *        considering whether HiRep is compiled WITH_GPU and/or WITH_MPI.
 *
 * @param _name                     Name of the field, this should describe
 *                                  the function this field has in the code
 *                                  given the dimension (_size)
 * @param _field_type               Type of field, use definitions in
 *                                  spinor_field.h
 * @param _site_type                Elementary site type from suN.h
 */
#define _DECLARE_FREE_FUNC(_name, _field_type, _site_type)                                                  \
    /** 
	  @brief Free field data, other struct fields and struct pointer of field struct. 
	  @param _field_type		Field to be freed. \
	 */ \
    void free_##_name(_field_type *f)                                                                       \
    {                                                                                                       \
        if (f != NULL)                                                                                      \
        {                                                                                                   \
            if (f->ptr != NULL)                                                                             \
                afree(f->ptr);                                                                              \
            _FREE_GPU_FIELD_DATA(_name, _site_type);                                                        \
            _FREE_MPI_FIELD_DATA;                                                                           \
            afree(f);                                                                                       \
            f = NULL;                                                                                       \
        }                                                                                                   \
    }

/**
 * @brief Generate allocation function for all fields defined in spinor_field.h
 *        considering whether HiRep is compiled WITH_GPU and/or WITH_MPI.
 *
 * @param _name                     Name of the field, this should describe
 *                                  the function this field has in the code
 *                                  given the dimension (_size)
 * @param _field_type               Type of field, use definitions in 
 *                                  spinor_field.h
 * @param _site_type                Elementary site type from suN.h
 */
#define _DECLARE_ALLOC_FUNC(_name, _field_type, _site_type, _size, _geom)                                   \
    /** 
	  @brief Allocate field struct pointer 
	  @param geometry_descriptor	Underlying lattice geometry to allocate on.
	 */ \
    _field_type *_ALLOCATE(_name)                                                                           \
    {                                                                                                       \
        _field_type *f;                                                                                     \
                                                                                                            \
        _ALLOC_FIELD_STRUCT(_name);                                                                         \
        _ALLOC_CPU_FIELD_DATA(_name, _size, _geom);                                                         \
        _ALLOC_GPU_FIELD_DATA(_name, _site_type, _size, _geom);                                             \
        _ALLOC_MPI_FIELD_DATA(_name, _size, _geom);                                                         \
                                                                                                            \
        return f;                                                                                           \
    }

/**
 * @brief Use this macro to declare both allocation and free functions
 *        for the given field type.
 *
 * @param _name                     Name of the field, this should describe the
 *                                  function this field has in the code
 *                                  given the dimension (_size)
 * @param _field_type               Type of field, use definitions in 
 *                                  spinor_field.h
 * @param _site_type                Elementary site type from suN.h
 * @param _size                     Number of elementary site types saved
 *                                  per site
 * @param _geom                     The geometries are slightly different depending 
 *                                  on whether the field is located on the sites
 *                                  or on the links. Fields that are located
 *                                  on the sites are "spinor-like", so put 
 *                                  'spinor' (without quotation marks), fields
 *                                  that are located on the links are "gauge-like"
 *                                  so put 'gauge'.
 */
#define _DECLARE_MEMORY_FUNC(_name, _field_type, _site_type, _size, _geom)                                  \
    _DECLARE_FREE_FUNC(_name, _field_type, _site_type)                                                      \
    _DECLARE_ALLOC_FUNC(_name, _field_type, _site_type, _size, _geom)                                              

/* Spinor allocation and free */
// We need to use the macros _n and _ALLOCATE here, because spinor-like 
// fields have the additional parameter n, that is also and argument to
// allocation.
#define _n n
#define _ALLOCATE(_name) alloc_##_name(unsigned int n, geometry_descriptor *type)
_DECLARE_MEMORY_FUNC(spinor_field_f, spinor_field, suNf_spinor, 1, spinor);
_DECLARE_MEMORY_FUNC(spinor_field_f_flt, spinor_field_flt, suNf_spinor_flt, 1, spinor);
_DECLARE_MEMORY_FUNC(sfield, scalar_field, double, 1, spinor);
#undef _n
#undef _ALLOCATE

/* Gauge allocation and free */
#define _n 1
#define _ALLOCATE(_name) alloc_##_name(geometry_descriptor *type) 
_DECLARE_MEMORY_FUNC(gfield, suNg_field, suNg, 4, gauge);
_DECLARE_MEMORY_FUNC(gfield_flt, suNg_field_flt, suNg_flt, 4, gauge);
_DECLARE_MEMORY_FUNC(gfield_f, suNf_field, suNf, 4, gauge);
_DECLARE_MEMORY_FUNC(gfield_f_flt, suNf_field_flt, suNf_flt, 4, gauge);
_DECLARE_MEMORY_FUNC(suNg_scalar_field, suNg_scalar_field, suNg_vector, 1, gauge);
_DECLARE_MEMORY_FUNC(avfield, suNg_av_field, suNg_algebra_vector, 4, gauge);
_DECLARE_MEMORY_FUNC(gtransf, suNg_field, suNg, 1, gauge);
_DECLARE_MEMORY_FUNC(clover_ldl, ldl_field, ldl_t, 1, gauge);
_DECLARE_MEMORY_FUNC(clover_term, suNfc_field, suNfc, 4, gauge);
_DECLARE_MEMORY_FUNC(clover_force, suNf_field, suNf, 6, gauge);
#undef _n
#undef _ALLOCATE

#undef _DECLARE_MEMORY_FUNC
#undef _DECLARE_FREE_FUNC
#undef _DECLARE_ALLOC_FUNC
#undef _ALLOC_FIELD_STRUCT
#undef _ALLOC_CPU_FIELD_DATA
#undef _ALLOC_GPU_FIELD_DATA
#undef _ALLOC_MPI_FIELD_DATA
#undef _FREE_GPU_FIELD_DATA
#undef _FREE_MPI_FIELD_DATA