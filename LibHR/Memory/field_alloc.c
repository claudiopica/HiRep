/***************************************************************************\
* Copyright (c) 2008, 2022, Claudio Pica, Sofie Martins                     *   
* All rights reserved.                                                      * 
\***************************************************************************/

/**
 * @file field_alloc.c
 * @brief Allocation and free functions for all fields defined in spinor_field.h
 */


#include "memory.h"
#include "libhr_core.h"

/**
 * @file alloc_cpu_field_data.c
 * @brief Allocation code snippets for the CPU/host that can be added to allocation 
 *        functions and are general to the field type.
 */

/**
 * @brief Allocate the poiner to the spinor field structure.
 *
 * @param _name                     Name of the field, this should describe
 *                                  the function this field has in the code
 *                                  given the dimension (_size)
 */
#define _ALLOC_FIELD_STRUCT(_name) \
    f = amalloc(_n*sizeof(*f), ALIGN);                                      \
    error(f == NULL, 1, "alloc_" #_name " [" __FILE__ "]",                  \
                    "Could not allocate memory space for field structure"); \
    for (int i = 0; i < _n; ++i) { f[i].type=type; } 


/**
 * @brief Allocate space for the field data of the local lattice on the CPU/host.
 *
 * @param _name                     Name of the field, this should describe
 *                                  the function this field has in the code
 *                                  given the dimension (_size)
 * @param _size                     Number of elementary site types saved
 *                                  per site
 * @param _geom                     The geometries are slightly different depending
 *                                  on whether the field is located on the sites or
 *                                  on the links. Fields that are located
 *                                  on the sites are "spinor-like", so put 
 *                                  'spinor' (without quotation marks), fields
 *                                  that are located on the links are "gauge-like"
 *                                  so put 'gauge'.
 */
#define _ALLOC_CPU_FIELD_DATA(_name, _size, _geom)                              \
    if (alloc_mem_t & CPU_MEM)                                                  \
    {                                                                           \
        /* For spinors: Allocate for all spinor array elements */               \
        int bytes_per_site = sizeof(*(f->ptr));                                 \
        int number_of_sites = _n * (_size) * type->gsize_##_geom;               \
        int field_size = bytes_per_site * number_of_sites;                      \
        f->ptr = amalloc(field_size, ALIGN);                                    \
        error((f->ptr) == NULL, 1, "alloc_" #_name " [" __FILE__ "]",           \
                    "Could not allocate memory space for field (data)");        \
        /* For spinors: Map the elements of the spinor arrays to the */         \
        /* starting points in the previously allocated space. */                \
        for (int i = 1; i < _n; ++i)                                            \
            f[i].ptr = f[i-1].ptr + type->gsize_##_geom * (_size);              \
    } else {                                                                    \
        for (int i = 0; i < _n; ++i) { f[i].ptr = NULL; }                       \
    }

/**
 * @file alloc_gpu_field_data.c
 * @brief Allocation code snippets for the GPU/device that can be added to allocation 
 *        functions and are general to the field type.
 */
#ifdef WITH_GPU

    /**
     * @brief Code snippet to free GPU field data.
     *
     * @param _name                 Name of the field, this should
     *                              describe the function this field has in
     *                              the code given the dimension (_size)
     * @param _site_type            Elementary site type from suN.h
     */
    #define _FREE_GPU_FIELD_DATA(_name, _site_type)  \
        if (f->gpu_ptr != NULL) { cudaFree(f->gpu_ptr); }

    /**
     * @brief Code snipped to allocate GPU field data.
     *
     * @param _name                 Name of the field, this should
     *                              describe the function this field has in
     *                              the code given the dimension (_size)
     * @param _site_type            Elementary site type from suN.h
     * @param _size                 Number of elementary site types saved
     *                              per site
     * @param _geom                 The geometries are slightly different depending 
     *                              on whether the field is located on the sites
     *                              or on the links. Fields that are located
     *                              on the sites are "spinor-like", so put 
     *                              'spinor' (without quotation marks), fields
     *                              that are located on the links are "gauge-like"
     *                              so put 'gauge'.
     */
    #define _ALLOC_GPU_FIELD_DATA(_name, _site_type, _size, _geom)             \
        if (alloc_mem_t & GPU_MEM)                                             \
        {                                                                      \
            cudaError_t err;                                                   \
            /* For spinors: Allocate for all spinor array elements */          \
            int bytes_per_site = sizeof(*(f->gpu_ptr));                        \
            int number_of_sites = _n * _size * type->gsize_##_geom;            \
            int field_size = number_of_sites * bytes_per_site;                 \
            err = cudaMalloc((void **)&(f->gpu_ptr), field_size);              \
            error(err != cudaSuccess, 1, "alloc_" #_name " [" __FILE__ "]",    \
                            "Could not allocate GPU memory space for field");  \
            /* For spinors: Map the elements of the spinor arrays to the */    \
            /* starting points in the previously allocated space. */           \
            for (int i = 1; i < _n; ++i)                                       \
                f[i].gpu_ptr = f[i-1].gpu_ptr + type->gsize_##_geom * _size;   \
        } else {                                                               \
            for (int i = 0; i < _n; ++i) { f[i].gpu_ptr = NULL; }              \
        }

#else

    /**
     * @brief Empty macro if code compiled without flag WITH_GPU
     */
    #define _FREE_GPU_FIELD_DATA(_name, _site_type) do {} while (0)
    /**
     * @brief Empty macro if code compiled without flag WITH_GPU
     */
    #define _ALLOC_GPU_FIELD_DATA(_name, _site_type, _size, _geom) do {} while (0)

#endif

/**
 * @file alloc_mpi_data.c
 * @brief Allocation of sendbuffers and communication handles necessary
 *        for MPI communications for both CPU and GPU.
 */

 //TODO: sendbuf_alloc analogy for GPU
 //TODO: Deallocation of senbuffers

#ifdef WITH_MPI

    #ifdef WITH_NEW_GEOMETRY
        #ifdef WITH_GPU
            #define _SENDBUF_ALLOC(_size, _i) \
                f[_i].sendbuf_ptr = sendbuf_alloc((_size)*sizeof(*(f[_i].ptr))); \
                f[_i].sendbuf_gpu_ptr = sendbuf_alloc_gpu((_size)*sizeof(*(f[_i].gpu_ptr)));
        #else
            #define _SENDBUF_ALLOC(_size, _i) \
                f[_i].sendbuf_ptr = sendbuf_alloc((_size)*sizeof(*(f[_i].ptr)));
        #endif
    #else
        #ifdef WITH_GPU
            #define _SENDBUF_ALLOC(_size, _i) \
                f[_i].sendbuf_ptr = f[_i].ptr; \
                f[_i].sendbuf_gpu_ptr = f[_i].gpu_ptr;
        #else
            #define _SENDBUF_ALLOC(_size, _i) \
                f[_i].sendbuf_ptr = f[_i].ptr; 
        #endif
    #endif

    /**
     * @brief Free memory allocated for MPI communications
     */
    #define _FREE_MPI_FIELD_DATA \
        if (f->comm_req != NULL) { afree(f->comm_req); }

    #define _ALLOC_MPI_FIELD_DATA(_name, _size, _geom)                                           \
        if (type->nbuffers_##_geom > 0)                                                          \
        {                                                                                        \
            f->comm_req = amalloc(_n * 2 * type->nbuffers_##_geom * sizeof(MPI_Request), ALIGN); \
            error((f->comm_req) == NULL, 1, "alloc_" #_name " [" __FILE__ "]",                   \
                "Could not allocate memory space for field (MPI)");                              \
            for (int ix = 0; ix < _n * 2 * type->nbuffers_##_geom; ++ix)                         \
                f->comm_req[ix] = MPI_REQUEST_NULL;                                              \
            for (int i = 1; i < _n; ++i)                                                         \
                f[i].comm_req = f[i-1].comm_req + 2 * type->nbuffers_##_geom;                    \
            for (int i = 0; i < _n; ++i) { _SENDBUF_ALLOC(_size, i); }                           \
        } else {                                                                                 \
            f->comm_req = NULL;                                                                  \
        }
  
#else
    /**
     * @brief Empty macro if compiled without WITH_MPI compilation flag
     */
    #define _FREE_MPI_FIELD_DATA do {} while (0)
    /**
     * @brief Empty macro if compiled without WITH_MPI compilation flag
     */ 
    #define _ALLOC_MPI_FIELD_DATA(_name, _size, _geom) do {} while (0)

#endif

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
#define _DECLARE_FREE_FUNC(_name, _field_type, _site_type) \
    /** 
	  @brief Free field data, other struct fields and struct pointer of field struct. 
	  @param _field_type		Field to be freed.
	 */ \
    void free_##_name(_field_type *f) {              \
        if (f != NULL) {                             \
            if (f->ptr != NULL)                      \
                afree(f->ptr);                       \
            _FREE_GPU_FIELD_DATA(_name, _site_type); \
            _FREE_MPI_FIELD_DATA;                    \
            afree(f);                                \
            f = NULL;                                \
        }                                            \
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
#define _DECLARE_ALLOC_FUNC(_name, _field_type, _site_type, _size, _geom) \
    /** 
	  @brief Allocate field struct pointer 
	  @param geometry_descriptor	Underlying lattice geometry to allocate on.
	 */ \
    _field_type *_ALLOCATE(_name) {                             \
        _field_type *f;                                         \
        _ALLOC_FIELD_STRUCT(_name);                             \
        _ALLOC_CPU_FIELD_DATA(_name, _size, _geom);             \
        _ALLOC_GPU_FIELD_DATA(_name, _site_type, _size, _geom); \
        _ALLOC_MPI_FIELD_DATA(_name, _size, _geom);             \
        return f;                                               \
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
#define _DECLARE_MEMORY_FUNC(_name, _field_type, _site_type, _size, _geom) \
    _DECLARE_FREE_FUNC(_name, _field_type, _site_type)                     \
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
_DECLARE_MEMORY_FUNC(staple_field, suNg_field, suNg, 3, gauge);
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