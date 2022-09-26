/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *   
* All rights reserved.                                                      * 
\***************************************************************************/

/*******************************************************************************
*
* File field_alloc.c
*
* Functions for fields allocation and memory operations
*
*******************************************************************************/

#include <stdlib.h>
#include "suN.h"
#include "error.h"
#include "memory.h"
#include "global.h"
#include "spinor_field.h"
#include "geometry.h"
#ifdef WITH_MPI
#include <mpi.h>
#endif

/*
 * Here we use for all macros the parameters:
 *
 *  _name = suffix of the allocation and deallocation functions
 *  _type = field type to allocate/deallocate
 *  _size = the number of elementary objects per lattice site
 *
 */

/* ================================================= MPI ================================================= */
/* MPI allocation and deallocation code */
#ifdef WITH_MPI

    #ifdef WITH_GPU
        // TODO:  Implement allocation for Cudaware MPI
    #else

        #define _FREE_MPI_CODE                                                                              \
            if (u->comm_req != NULL)                                                                        \
            afree(u->comm_req)

        #define _ALLOC_MPI_CODE(_name)                                                                      \
            if (type->nbuffers_gauge > 0)                                                                   \
            {                                                                                               \
                f->comm_req = amalloc(2 * type->nbuffers_gauge * sizeof(MPI_Request), ALIGN);               \
                error((f->comm_req) == NULL, 1, "alloc_" #_name " [" __FILE__ "]",                          \
                    "Could not allocate memory space for field (MPI)");                                     \
                for (int ix = 0; ix < 2 * type->nbuffers_gauge; ++ix)                                       \
                    f->comm_req[ix] = MPI_REQUEST_NULL;                                                     \
            }                                                                                               \
            else                                                                                            \
            {                                                                                               \
                f->comm_req = NULL;                                                                         \
            }                                                                                               \
            do {} while (0)                                                                                 \
  
    #endif

#else /* WITH_MPI */

    #define _FREE_MPI_CODE do {} while (0)
    #define _ALLOC_MPI_CODE(_name) do {} while (0) 

#endif /* WITH_MPI */

/* ================================================= GPU ================================================= */
#ifdef WITH_GPU

    /* Free device memory */
    /* Note: to be used inside function declaration */
    #define _FREE_GPU_CODE                                                                                  \
        if (u->gpu_ptr != NULL)                                                                             \
            cudaFree(u->gpu_ptr)

    /* Allocate device memory */
    /* Note: to be used inside function declaration */
    #define _ALLOC_GPU_CODE(_name, _size)                                                                   \
        if (alloc_mem_t & GPU_MEM)                                                                          \
        {                                                                                                   \
            cudaError_t err;                                                                                \
            int field_size = _size * type->gsize_gauge * sizeof(*(f->gpu_ptr));                             \
            err = cudaMalloc((void **)&f->gpu_ptr, field_size);                                             \
            error(err != cudaSuccess, 1, "alloc_" #_name " [" __FILE__ "]",                                 \
                            "Could not allocate GPU memory space for field");                               \
        }                                                                                                   \
        else                                                                                                \
            f->gpu_ptr = NULL

    /* Declare function to copy field from host to device */
    #define _DECLARE_COPY_TO(_name, _type, _size)                                                           \
        void copy_to_gpu_##_name(_type *u)                                                                  \
        { \
            _type *tmp = alloc_##_name(u->type);                                                            \
            to_gpu_format_##_name(tmp, u);                                                                  \
            int field_size = _size * u->type->gsize_gauge * sizeof(*(u->gpu_ptr));                          \
            cudaMemcpy(u->gpu_ptr, tmp->ptr, field_size, cudaMemcpyHostToDevice);                           \
            free_##_name(tmp);                                                                              \
        }

    /* Declare function to copy field from device to host */
    #define _DECLARE_COPY_FROM(_name, _type, _size)                                                         \
        void copy_from_gpu_##_name(_type *u)                                                                \
        {                                                                                                   \
            _type *tmp = alloc_##_name(u->type);                                                            \
            int field_size = _size * u->type->gsize_gauge * sizeof(*(u->gpu_ptr));                          \
            cudaMemcpy(tmp->ptr, u->gpu_ptr, field_size, cudaMemcpyDeviceToHost);                           \
            to_cpu_format_##_name(u, tmp);                                                                  \
            free_##_name(tmp);                                                                              \
        }

    /* Declare function to convert GPU to CPU format */
    /* This is necessary, because GPU performance is very sensitive to */
    /* memory access patterns, see documentation Doc/gpu_geometry.tex */
    #define _DECLARE_CONVERT_TO_GPU_FORMAT(_name, _field_type, _site_type, _size)                           \
        void to_gpu_format_##_name(_type *out, _type *in)                                                   \
        {                                                                                                   \
            _site_type *r = 0;                                                                              \
            int number_of_elements;                                                                         \
            error(out->type!=in->type, 1, "to_gpu_format_" ##_name " " __FILE__,                            \
                    "Gauge field geometries do not match!\n");                                              \
            _PIECE_FOR(in->type, ixp)                                                                       \
            {                                                                                               \
                const int start = in->type->master_start[ixp];                                              \
                const int N = in->type->master_end[ixp] - in->type->master_start[ixp]+1;                    \
                double *cout = (double*)(_4FIELD_AT(out, start, 0));                                        \
                _SITE_FOR(in->type, ixp, ix)                                                                \
                {                                                                                           \
                    r = _4FIELD_AT(in, ix, 0);                                                              \
                                                                                                            \
                    number_of_elements = _size * sizeof(*r)/sizeof(double);                                 \
                    for (int j = 0; j < number_of_elements, ++j)                                            \
                    {                                                                                       \
                        cout[j*N] = ((double*)(r))[j];                                                      \
                    }                                                                                       \
                    ++cout;                                                                                 \
                }                                                                                           \
            }                                                                                               \
        }

    #define _DECLARE_CONVERT_TO_CPU_FORMAT(_name, _field_type, _site_type, _size)                           \
        void to_cpu_format_##_name(_type *out, _type *in)                                                   \
        {                                                                                                   \
            _site_type *r = 0;                                                                              \
            int number_of_elements;                                                                         \
            error(out->type!=in->type, 1, "to_cpu_format_" ##_name " " __FILE__,                            \
                        "Gauge field geometries do not match!\n");                                          \
            _PIECE_FOR(in->type, ixp)                                                                       \
            {                                                                                               \
                const int start = in->type->master_start[ixp];                                              \
                const int N = in->type->master_end[ixp] - in->type->master_start[ixp] + 1;                  \
                double *cin = (double*)(_4FIELD_AT(in, start, 0));                                          \
                _SITE_FOR(in->type, ixp, ix)                                                                \
                {                                                                                           \
                    r = _4FIELD_AT(out, ix, 0);                                                             \
                                                                                                            \
                    number_of_elements = _size * sizeof(*r)/sizeof(double);                                 \
                    for (int j = 0; j < number_of_elements; ++j)                                            \
                    {                                                                                       \
                        ((double*)(r))[j] = cin[j*N];                                                       \
                    }                                                                                       \
                    ++cin;                                                                                  \
                }                                                                                           \
            }                                                                                               \
        }

#else /* WITH_GPU */

    #define _FREE_GPU_CODE do {} while (0)
    #define _ALLOC_GPU_CODE(_name, _size) do {} while (0)
    #define _DECLARE_COPY_TO(_name, _size)
    #define _DECLARE_COPY_FROM(_name, _size) 
    #define _DECLARE_CONVERT_TO_GPU_FORMAT(_name, _field_type, _site_type, _size) 
    #define _DECLARE_CONVERT_TO_CPU_FORMAT(_name, _field_type, _site_type, _size)

#endif /* WITH_GPU */

/* ============================================== All cases ============================================== */

/* deallocation function declaration */
#define _DECLARE_FREE_FUNC(_name, _type)                                                                    \
    void free_##_name(_type *u)                                                                             \
    {                                                                                                       \
        if (u != NULL)                                                                                      \
        {                                                                                                   \
            if (u->ptr != NULL)                                                                             \
                afree(u->ptr);                                                                              \
            _FREE_GPU_CODE;                                                                                 \
            _FREE_MPI_CODE;                                                                                 \
            afree(u);                                                                                       \
            u = NULL;                                                                                       \
        }                                                                                                   \
    }

/* allocation function declaration */
#define _DECLARE_ALLOC_FUNC(_name, _type, _size)                                                            \
    _type *alloc_##_name(geometry_descriptor *type)                                                         \
    {                                                                                                       \
        /* Allocate field struct pointer */                                                                 \
        _type *f;                                                                                           \
                                                                                                            \
        f = amalloc(sizeof(*f), ALIGN);                                                                     \
        error(f == NULL, 1, "alloc_" #_name " [" __FILE__ "]",                                              \
                "Could not allocate memory space for field (structure)");                                   \
        f->type = type;                                                                                     \
                                                                                                            \
        /* Allocate CPU field data */                                                                       \
        if (alloc_mem_t & CPU_MEM)                                                                          \
        {                                                                                                   \
            int field_size = _size * type->gsize_gauge * sizeof(*(f->ptr), ALIGN);                          \
            f->ptr = amalloc(field_size, ALIGN);                                                            \
            error((f->ptr) == NULL, 1, "alloc_" #_name " [" __FILE__ "]",                                   \
                            "Could not allocate memory space for field (data)");                            \
        }                                                                                                   \
        else f->ptr = NULL;                                                                                 \
                                                                                                            \
        /* Allocate GPU field, if compiling with GPU */                                                     \
        _ALLOC_GPU_CODE(_name, _size);                                                                      \
                                                                                                            \
        /* Allocate buffers for MPI comms, if compiling with MPI */                                         \
        _ALLOC_MPI_CODE(_name);                                                                             \
                                                                                                            \
        return f;                                                                                           \
    }



#define _DECLARE_MEMORY_FUNC(_name, _type, _site_type, _size)                                               \
    _DECLARE_FREE_FUNC(_name, _type)                                                                        \
    _DECLARE_ALLOC_FUNC(_name, _type, _size)                                                                \
    _DECLARE_COPY_TO(_name, _type, _size)                                                                   \
    _DECLARE_COPY_FROM(_name, _type, _size)                                                                 \
    _DECLARE_CONVERT_TO_GPU_FORMAT(_name, _type, _site_type, _size)                                         \
    _DECLARE_CONVERT_TO_CPU_FORMAT(_name, _type, _site_type, size)

_DECLARE_MEMORY_FUNC(gfield, suNg_field, suNg, 4);
_DECLARE_MEMORY_FUNC(gfield_flt, suNg_field_flt, suNg_flt, 4);
_DECLARE_MEMORY_FUNC(gfield_f, suNf_field, suNf, 4);
_DECLARE_MEMORY_FUNC(gfield_f_flt, suNf_field_flt, suNf_flt, 4);
_DECLARE_MEMORY_FUNC(scalar_field, suNg_scalar_field, suNf_spinor, 1);
_DECLARE_MEMORY_FUNC(avfield, suNg_av_field, suNg_algebra_vector, 4);
_DECLARE_MEMORY_FUNC(gtransf, suNg_field, suNg, 1);
_DECLARE_MEMORY_FUNC(clover_ldl, ldl_field, ldl_t, 1);
_DECLARE_MEMORY_FUNC(clover_term, suNfc_field, suNfc, 4);
_DECLARE_MEMORY_FUNC(clover_force, suNf_field, suNf, 6);

#undef _DECLARE_MEMORY_FUNC
