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

// TODO: Doxygen docs
// TODO: Remove brackets from macro params
// TODO: spinor n parameter
// TODO: block handle allocation right? It seems that we have a handle for every block/process combination, 
//       which we do not need, because every process only deals with a subset of blocks

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
            #define _SENDBUF_ALLOC(_size) \
                    f->sendbuf_ptr = sendbuf_alloc((_size)*sizeof(*(f->ptr))); \
                    int alloc_length = (_size)*sizeof(*(f->ptr))*(glattice.gsize_gauge - boxVolume(geometryBoxes)); \
                    cudaMalloc((void **)&(f->sendbuf_gpu_ptr), alloc_length);

        #else
            #define _SENDBUF_ALLOC(_size) \
                    f->sendbuf_ptr = f->ptr; \
                    f->sendbuf_gpu_ptr = f->gpu_ptr;

        #endif


        #define _FREE_MPI_FIELD_DATA                                                                        \
            if (f->comm_req != NULL)                                                                        \
                afree(f->comm_req)                                                                          \

        #define _ALLOC_MPI_FIELD_DATA(_name, _size)                                                         \
            if (type->nbuffers_gauge > 0)                                                                   \
            {                                                                                               \
                f->comm_req = amalloc(2 * type->nbuffers_gauge * sizeof(MPI_Request), ALIGN);               \
                error((f->comm_req) == NULL, 1, "alloc_" #_name " [" __FILE__ "]",                          \
                    "Could not allocate memory space for field (MPI)");                                     \
                for (int ix = 0; ix < 2 * type->nbuffers_gauge; ++ix)                                       \
                    f->comm_req[ix] = MPI_REQUEST_NULL;                                                     \
                _SENDBUF_ALLOC(_size);                                                                      \
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
        #define _ALLOC_GPU_FIELD_DATA(_name, _site_type, _size)                                             \
            if (alloc_mem_t & GPU_MEM)                                                                      \
            {                                                                                               \
                cudaError_t err;                                                                            \
                int field_size = _size * type->gsize_gauge * sizeof(*(f->gpu_ptr));                         \
                err = cudaMalloc((void **)&(f->gpu_ptr), field_size);                                       \
                error(err != cudaSuccess, 1, "alloc_" #_name " [" __FILE__ "]",                             \
                                "Could not allocate GPU memory space for field");                           \
            }                                                                                               \
            else                                                                                            \
                f->gpu_ptr = NULL;                                                                          \

        /* Declare function to copy field from host to device */
        #define _DECLARE_COPY_TO(_name, _field_type, _site_type, _size)                                     \
            void copy_to_gpu_##_name(_field_type *f)                                                        \
            {                                                                                               \
                _field_type *tmp = alloc_##_name(f->type);                                                  \
                to_gpu_format_##_name(tmp, f);                                                              \
                int field_size = _size * f->type->gsize_gauge * sizeof(*(f->gpu_ptr));                      \
                cudaMemcpy(f->gpu_ptr, tmp->ptr, field_size, cudaMemcpyHostToDevice);                       \
                free_##_name(tmp);                                                                          \
            }

        /* Declare function to copy field from device to host */
        #define _DECLARE_COPY_FROM(_name, _field_type, _site_type, _size)                                   \
            void copy_from_gpu_##_name(_field_type *f)                                                      \
            {                                                                                               \
                _field_type *tmp = alloc_##_name(f->type);                                                  \
                int field_size = _size * f->type->gsize_gauge * sizeof(*(f->gpu_ptr));                      \
                cudaMemcpy(tmp->ptr, f->gpu_ptr, field_size, cudaMemcpyDeviceToHost);                       \
                to_cpu_format_##_name(f, tmp);                                                              \
                free_##_name(tmp);                                                                          \
            }

    /* Declare function to convert GPU to CPU format */
    /* This is necessary, because GPU performance is very sensitive to */
    /* memory access patterns, see documentation on GPU Geometry */
    #define _DECLARE_CONVERT_TO_GPU_FORMAT(_name, _field_type, _site_type, _size)                           \
        void to_gpu_format_##_name(_field_type *out, _field_type *in)                                       \
        {                                                                                                   \
            _site_type *source, *target;                                                                    \
            _CHECK_GEOMETRY_MATCHING(out, in);                                                              \
                                                                                                            \
            int stride = 0;                                                                                 \
            int ix_loc = 0;                                                                                 \
            _PIECE_FOR(in->type, ixp)                                                                       \
            {                                                                                               \
                stride = in->type->master_end[ixp] - in->type->master_start[ixp] + 1;                       \
                target = _DFIELD_BLK(out, ixp, (_size));                                                    \
                _SITE_FOR(in->type, ixp, ix)                                                                \
                {                                                                                           \
                    ix_loc = _GPU_IDX_TO_LOCAL(in, ix, ixp);                                                \
                    for (int comp = 0; comp < _size; ++comp)                                                \
                    {                                                                                       \
                        source = _DFIELD_AT(in, ix, comp, (_size));                                         \
                        write_gpu_##_site_type(stride, (*source), target, ix_loc, comp);                    \
                    }                                                                                       \
                }                                                                                           \
            }                                                                                               \
        }

    #define _DECLARE_CONVERT_TO_CPU_FORMAT(_name, _field_type, _site_type, _size)                           \
        void to_cpu_format_##_name(_field_type *out, _field_type *in)                                       \
        {                                                                                                   \
            _site_type *target, *source;                                                                    \
            _CHECK_GEOMETRY_MATCHING(out, in);                                                              \
                                                                                                            \
            int ix_loc = 0;                                                                                 \
            int stride = 0;                                                                                 \
            _PIECE_FOR(in->type, ixp)                                                                       \
            {                                                                                               \
                stride = in->type->master_end[ixp] - in->type->master_start[ixp] + 1;                       \
                source = _DFIELD_BLK(in, ixp, (_size));                                                     \
                _SITE_FOR(in->type, ixp, ix)                                                                \
                {                                                                                           \
                    ix_loc = _GPU_IDX_TO_LOCAL(in, ix, ixp);                                                \
                    for (int comp = 0; comp < _size; ++comp)                                                \
                    {                                                                                       \
                        target = _DFIELD_AT(out, ix, comp, (_size));                                        \
                        read_gpu_##_site_type(stride, (*target), source, ix_loc, comp);                     \
                    }                                                                                       \
                }                                                                                           \
            }                                                                                               \
        }

#endif

/* ================================================= Empty defs ========================================== */

#ifndef WITH_GPU

    #define _FREE_GPU_FIELD_DATA(_name, _site_type) do {} while (0)
    #define _ALLOC_GPU_FIELD_DATA(_name, _site_type, _size) do {} while (0)
    #define _DECLARE_COPY_TO(_name, _field_type, _site_type, _size)
    #define _DECLARE_COPY_FROM(_name, _field_type, _site_type, _size)
    #define _DECLARE_CONVERT_TO_GPU_FORMAT(_name, _field_type, _site_type, _size) 
    #define _DECLARE_CONVERT_TO_CPU_FORMAT(_name, _field_type, _site_type, _size)

#endif

#ifndef WITH_MPI

    #define _FREE_MPI_FIELD_DATA do {} while (0)
    #define _ALLOC_MPI_FIELD_DATA(_name, _size) do {} while (0)

#endif

/* ============================================== CPU ==================================================== */

#define _ALLOC_FIELD_STRUCT(_name)                                                                          \
    f = amalloc(sizeof(*f), ALIGN);                                                                         \
        error(f == NULL, 1, "alloc_" #_name " [" __FILE__ "]",                                              \
                    "Could not allocate memory space for field (structure)");                               \
        f->type = type;

#define _ALLOC_CPU_FIELD_DATA(_name, _size)                                                                 \
    if (alloc_mem_t & CPU_MEM)                                                                              \
    {                                                                                                       \
        int field_size = (_size) * type->gsize_gauge * sizeof(*(f->ptr));                                   \
        f->ptr = amalloc(field_size, ALIGN);                                                                \
        error((f->ptr) == NULL, 1, "alloc_" #_name " [" __FILE__ "]",                                       \
                    "Could not allocate memory space for field (data)");                                    \
    }                                                                                                       \
    else f->ptr = NULL;

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

#define _DECLARE_ALLOC_FUNC(_name, _field_type, _site_type, _size)                                          \
    _field_type *alloc_##_name(geometry_descriptor *type)                                                   \
    {                                                                                                       \
        _field_type *f;                                                                                     \
                                                                                                            \
        _ALLOC_FIELD_STRUCT(_name);                                                                         \
        _ALLOC_CPU_FIELD_DATA(_name, _size);                                                                \
        _ALLOC_GPU_FIELD_DATA(_name, _site_type, _size);                                                    \
        _ALLOC_MPI_FIELD_DATA(_name, _size);                                                                       \
                                                                                                            \
        return f;                                                                                           \
    }

#define _DECLARE_MEMORY_FUNC(_name, _field_type, _site_type, _size)                                         \
    _DECLARE_FREE_FUNC(_name, _field_type, _site_type)                                                      \
    _DECLARE_ALLOC_FUNC(_name, _field_type, _site_type, _size)                                              \
    _DECLARE_COPY_TO(_name, _field_type, _site_type, _size)                                                 \
    _DECLARE_COPY_FROM(_name, _field_type, _site_type, _size)                                               \
    _DECLARE_CONVERT_TO_GPU_FORMAT(_name, _field_type, _site_type, _size)                                   \
    _DECLARE_CONVERT_TO_CPU_FORMAT(_name, _field_type, _site_type, _size)

_DECLARE_MEMORY_FUNC(gfield, suNg_field, suNg, 4);
_DECLARE_MEMORY_FUNC(gfield_flt, suNg_field_flt, suNg_flt, 4);
_DECLARE_MEMORY_FUNC(gfield_f, suNf_field, suNf, 4);
_DECLARE_MEMORY_FUNC(gfield_f_flt, suNf_field_flt, suNf_flt, 4);
_DECLARE_MEMORY_FUNC(scalar_field, suNg_scalar_field, suNg_vector, 1);
_DECLARE_MEMORY_FUNC(avfield, suNg_av_field, suNg_algebra_vector, 4);
_DECLARE_MEMORY_FUNC(gtransf, suNg_field, suNg, 1);
//_DECLARE_MEMORY_FUNC(clover_ldl, ldl_field, ldl_t, 1);
_DECLARE_MEMORY_FUNC(clover_term, suNfc_field, suNfc, 4);
_DECLARE_MEMORY_FUNC(clover_force, suNf_field, suNf, 6);

#undef _DECLARE_MEMORY_FUNC
#undef _DECLARE_FREE_FUNC
#undef _DECLARE_ALLOC_FUNC
#undef _DECLARE_COPY_TO
#undef _DECLARE_COPY_FROM
#undef _DECLARE_CONVERT_TO_GPU_FORMAT
#undef _DECLARE_CONVERT_TO_CPU_FORMAT
#undef _ALLOC_FIELD_STRUCT
#undef _ALLOC_CPU_FIELD_DATA
#undef _ALLOC_GPU_FIELD_DATA
#undef _ALLOC_MPI_FIELD_DATA
#undef _FREE_GPU_FIELD_DATA
#undef _FREE_MPI_FIELD_DATA