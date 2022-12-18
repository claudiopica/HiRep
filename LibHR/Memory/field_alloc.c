/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *   
* All rights reserved.                                                      * 
\***************************************************************************/

/*******************************************************************************
*
* File gauge_field_alloc.c
*
* Functions for gauge field allocation and memory operations
*
*******************************************************************************/

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

// TODO: doxygen docstrings
// TODO: inline comments
// TODO: More encapsulation?

/*
 * Here we use for all macros the parameters:
 *
 *  _name = suffix of the allocation and deallocation functions
 *  _field_type = field type to allocate/deallocate
 *  _site_type = type of elementary objects on the lattice sites
 *  _size = the number of elementary objects per lattice site
 *
 */

/* ================================================= MPI ================================================= */

#if defined(WITH_MPI)

        #define _FREE_MPI_CODE  if (u->comm_req != NULL) afree(u->comm_req)

        #ifdef WITH_NEW_GEOMETRY
            #ifdef WITH_GPU
                #define _SENDBUF_ALLOC(_size, _i) \
                /*TODO: GPU sendbuf not allocated correctly. (SAM)*/ \
                    f[_i].sendbuf_ptr = sendbuf_alloc((_size)*sizeof(*(f[_i].ptr))); \
                    int alloc_length = (_size)*sizeof(*(f[_i].ptr))*(glattice.gsize_gauge - boxVolume(geometryBoxes)); \
                    cudaMalloc((void **)&(f[_i].sendbuf_gpu_ptr), alloc_length);
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


        #define _FREE_MPI_FIELD_DATA                                                                        \
            if (f->comm_req != NULL)                                                                        \
                afree(f->comm_req)                                                                          \

        #define _ALLOC_MPI_FIELD_DATA(_name, _size, _geom)                                                  \
            if (type->nbuffers_##_geom > 0)                                                                 \
            {                                                                                               \
                f->comm_req = amalloc(_n * 2 * type->nbuffers_##_geom * sizeof(MPI_Request), ALIGN);        \
                error((f->comm_req) == NULL, 1, "alloc_" #_name " [" __FILE__ "]",                          \
                    "Could not allocate memory space for field (MPI)");                                     \
                for (int ix = 0; ix < _n * 2 * type->nbuffers_##_geom; ++ix)                                \
                    f->comm_req[ix] = MPI_REQUEST_NULL;                                                     \
                for (int i = 1; i < _n; ++i)                                                                 \
                {                                                                                           \
                    f[i].comm_req = f[i-1].comm_req + 2 * type->nbuffers_##_geom;                           \
                }                                                                                           \
                for (int i = 1; i < _n; ++i)                                                                 \
                {                                                                                           \
                    _SENDBUF_ALLOC(_size, i);                                                               \
                }                                                                                           \
                                                                                                            \
            }                                                                                               \
            else                                                                                            \
            {                                                                                               \
                f->comm_req = NULL;                                                                         \
            }                                                                                               \
  
#endif

/* ================================================= GPU ================================================= */

#if defined(WITH_GPU)

        /* Free device memory */
        /* Note: to be used inside function declaration */
        #define _FREE_GPU_FIELD_DATA(_name, _site_type)                                                     \
            if (f->gpu_ptr != NULL)                                                                         \
                cudaFree(f->gpu_ptr);                                                                       \

        /* Allocate device memory */
        /* Note: to be used inside function declaration */
        #define _ALLOC_GPU_FIELD_DATA(_name, _site_type, _size, _geom)                                             \
            if (alloc_mem_t & GPU_MEM)                                                                      \
            {                                                                                               \
                cudaError_t err;                                                                            \
                int field_size = _n * _size * type->gsize_##_geom * sizeof(*(f->gpu_ptr));                         \
                err = cudaMalloc((void **)&(f->gpu_ptr), field_size);                                       \
                error(err != cudaSuccess, 1, "alloc_" #_name " [" __FILE__ "]",                             \
                                "Could not allocate GPU memory space for field");                           \
                for (int i = 1; i < _n; ++i)                                                                \
                    f[i].gpu_ptr = f[i-1].gpu_ptr + type->gsize_spinor * _size;                             \
            }                                                                                               \
            else                                                                                            \
                for (int i = 0; i < _n; ++i)                                                                \
                    f[i].gpu_ptr = NULL;

#endif

/* ================================================= Empty defs ========================================== */

#ifndef WITH_GPU

    #define _FREE_GPU_FIELD_DATA(_name, _site_type) do {} while (0)
    #define _ALLOC_GPU_FIELD_DATA(_name, _site_type, _size, _geom) do {} while (0)

#endif

#ifndef WITH_MPI

    #define _FREE_MPI_FIELD_DATA do {} while (0)
    #define _ALLOC_MPI_FIELD_DATA(_name, _size, _geom) do {} while (0)

#endif

/* ============================================== CPU ==================================================== */

#define _ALLOC_FIELD_STRUCT(_name)                                                                   \
    f = amalloc(sizeof(*f), ALIGN);                                                                         \
        error(f == NULL, 1, "alloc_" #_name " [" __FILE__ "]",                                              \
                    "Could not allocate memory space for field (structure)");                               \
        f->type = type;

#define _ALLOC_CPU_FIELD_DATA(_name, _size, _geom)                                                                 \
    if (alloc_mem_t & CPU_MEM)                                                                              \
    {                                                                                                       \
        int field_size = _n * (_size) * type->gsize_##_geom * sizeof(*(f->ptr));                            \
        f->ptr = amalloc(field_size, ALIGN);                                                                \
        error((f->ptr) == NULL, 1, "alloc_" #_name " [" __FILE__ "]",                                       \
                    "Could not allocate memory space for field (data)");                                    \
        for (int i = 1; i < _n; ++i)                                                                        \
            f[i].ptr = f[i-1].ptr + type->gsize_##_geom * (_size);                                           \
    }                                                                                                       \
    else \
    { \
        for (int i = 0; i < _n; ++i)  \
        { \
            f[i].ptr = NULL; \
        } \
    }

/* ============================================== All cases ============================================== */

#define _DECLARE_FREE_FUNC(_name, _field_type, _site_type)                                                  \
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

#define _DECLARE_ALLOC_FUNC(_name, _field_type, _site_type, _size, _geom)                                   \
    _field_type *_ALLOCATE(_name)                                                                                  \
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

#define _DECLARE_MEMORY_FUNC(_name, _field_type, _site_type, _size, _geom)                                  \
    _DECLARE_FREE_FUNC(_name, _field_type, _site_type)                                                      \
    _DECLARE_ALLOC_FUNC(_name, _field_type, _site_type, _size, _geom)                                              

#define _n n
#define _ALLOCATE(_name) alloc_##_name(unsigned int n, geometry_descriptor *type)
_DECLARE_MEMORY_FUNC(spinor_field_f, spinor_field, suNf_spinor, 1, spinor);
_DECLARE_MEMORY_FUNC(spinor_field_f_flt, spinor_field_flt, suNf_spinor_flt, 1, spinor);
_DECLARE_MEMORY_FUNC(sfield, scalar_field, double, 1, spinor);
#undef _n
#undef _ALLOCATE

#define _n 1
#define _ALLOCATE(_name) alloc_##_name(geometry_descriptor *type) 
_DECLARE_MEMORY_FUNC(gfield, suNg_field, suNg, 4, gauge);
_DECLARE_MEMORY_FUNC(gfield_flt, suNg_field_flt, suNg_flt, 4, gauge);
_DECLARE_MEMORY_FUNC(gfield_f, suNf_field, suNf, 4, gauge);
_DECLARE_MEMORY_FUNC(gfield_f_flt, suNf_field_flt, suNf_flt, 4, gauge);
_DECLARE_MEMORY_FUNC(scalar_field, suNg_scalar_field, suNg_vector, 1, gauge);
_DECLARE_MEMORY_FUNC(avfield, suNg_av_field, suNg_algebra_vector, 4, gauge);
_DECLARE_MEMORY_FUNC(gtransf, suNg_field, suNg, 1, gauge);
//_DECLARE_MEMORY_FUNC(clover_ldl, ldl_field, ldl_t, 1);
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