/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *
* All rights reserved.                                                      *
\***************************************************************************/

/*******************************************************************************
*
* File field_alloc.c
*
* Functions for matter field allocation and memory operations
*
*******************************************************************************/

#include <stdlib.h>
#include "suN.h"
#include "error.h"
#include "memory.h"
#include "geometry.h"
#include "linear_algebra.h"
#include "global.h"
#ifdef WITH_MPI
#include <mpi.h>
#endif

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
/* MPI allocation and deallocation code */
#ifdef WITH_MPI

    #ifdef WITH_GPU
        //TODO: Implement allocation for CUDA-aware MPI
    #else

        #define _FREE_MPI_CODE                                                                              \
            if(u->comm_req != NULL)                                                                         \
            afree(u->comm_req)

        #define _ALLOC_MPI_CODE(_name)                                                                      \
            if (type->nbuffers_spinor > 0)                                                                  \
            {                                                                                               \
                f->comm_req = amalloc( n * 2 * type->nbuffers_spinor * sizeof(MPI_Request), ALIGN);         \
                error((f->comm_req) == NULL, 1, "alloc_" #_name " [" __FILE__ "]",                          \
                    "Could not allocate memory space for field (MPI)");                                     \
                for (int ix = 0; ix < n * 2 * type->nbuffers_spinor; ++ix)                                  \
                    f->comm_req[ix]=MPI_REQUEST_NULL;                                                       \
                for (int i = 1; i < n; ++i) f[i].comm_req=f[i-1].comm_req + 2 * type->nbuffers_spinor;      \
            }                                                                                               \
            else                                                                                            \
            {                                                                                               \
                for (int i = 0; i < n; ++i) f[i].comm_req = NULL;                                           \
            }                                                                                               \
            do {} while (0)
    #endif

#else /* WITH_MPI */

    #define _FREE_MPI_CODE do {} while(0)
    #define _ALLOC_MPI_CODE(_name) do {} while(0)

#endif /* WITH_MPI */

/* ================================================= GPU ================================================= */
/* GPU allocation and deallocation code */
#ifdef WITH_GPU

    /* Free device memory */
    /* Note: to be used inside function declaration */
    #define _FREE_GPU_CODE                                                                                  \
        if(u->gpu_ptr != NULL)                                                                              \
            cudaFree(u->gpu_ptr);                                                                           \
        do {} while (0)

    #define _ALLOC_GPU_CODE(_name,_size)                                                                    \
        if(alloc_mem_t & GPU_MEM)                                                                           \
        {                                                                                                   \
            cudaError_t err;                                                                                \
            int field_size = n * _size * type->gsize_spinor * sizeof(*(f->gpu_ptr));                        \
            err = cudaMalloc((void **) &(f->gpu_ptr), field_size);                                          \
            error(err != cudaSuccess, 1, "alloc_" #_name " [" __FILE__ "]",                                 \
                                "Could not allocate GPU memory space for field");                           \
            for (int i = 1; i < n; ++i)                                                                     \
                f[i].gpu_ptr = f[i-1].gpu_ptr + type->gsize_spinor * _size;                                 \
        } else                                                                                              \
            for (int i = 0; i < n; ++i)                                                                     \
                f[i].gpu_ptr = NULL;                                                                        \
        do {} while (0);

    /* Declare function to copy field from host to device */
    #define _DECLARE_COPY_TO(_name, _field_type, _size)                                                     \
        void copy_to_gpu_##_name(_field_type *f)                                                            \
        {                                                                                                   \
            /* FIXME: Do not only copy one layer. */                                                        \
            _field_type *tmp = alloc_##_name(1, f->type);                                                   \
            to_gpu_format_##_name(tmp, f);                                                                  \
            int field_size = _size * f->type->gsize_gauge * sizeof(*(f->gpu_ptr));                          \
            cudaMemcpy(f->gpu_ptr, tmp->ptr, field_size, cudaMemcpyHostToDevice);                           \
            free_##_name(tmp);                                                                              \
        }

    /* Declare function to copy field from device to host */
    #define _DECLARE_COPY_FROM(_name, _field_type, _size)                                                   \
        void copy_from_gpu_##_name(_field_type *f)                                                          \
        {                                                                                                   \
            _field_type *tmp = alloc_##_name(1, f->type);                                                   \
            int field_size = _size * f->type->gsize_gauge * sizeof(*(f->gpu_ptr));                          \
            cudaMemcpy(tmp->ptr, f->gpu_ptr, field_size, cudaMemcpyDeviceToHost);                           \
            to_cpu_format_##_name(f, tmp);                                                                  \
            free_##_name(tmp);                                                                              \
        }

    /* Declare function to convert GPU to CPU format */
    /* This is necessary, because GPU performance is very sensitive to */
    /* memory access patterns, see documentation Doc/gpu_geometry.tex */
    #define _DECLARE_CONVERT_TO_GPU_FORMAT(_name, _field_type, _site_type, _size)                           \
        void to_gpu_format_##_name(_field_type *out, _field_type *in)                                       \
        {                                                                                                   \
            _site_type *r = 0;                                                                              \
            int number_of_elements;                                                                         \
            error(out->type != in->type, 1, "to_gpu_format_" #_name " " __FILE__,                           \
                    "Gauge field geometries do not match!\n");                                              \
            _PIECE_FOR(in->type, ixp)                                                                       \
            {                                                                                               \
                const int start = in->type->master_start[ixp];                                              \
                const int N = in->type->master_end[ixp] - in->type->master_start[ixp];                      \
                /* Does not work for the single precision types */                                          \
                hr_complex *cout = (hr_complex*)(_FIELD_AT(out, start));                                    \
                _SITE_FOR(in->type, ixp, ix)                                                                 \
                {                                                                                           \
                    r = _FIELD_AT(in, ix);                                                                   \
                                                                                                            \
                    number_of_elements = sizeof(*r)/sizeof(hr_complex);                                     \
                    for (int j = 0; j < number_of_elements; ++j)                                            \
                    {                                                                                       \
                        cout[j*N] = ((hr_complex*)(r))[j];                                                  \
                    }                                                                                       \
                    ++cout;                                                                                 \
                }                                                                                           \
            }                                                                                               \
        }

    #define _DECLARE_CONVERT_TO_CPU_FORMAT(_name, _field_type, _site_type, _size)                           \
        void to_cpu_format_##_name(_field_type *out, _field_type *in)                                       \
        {                                                                                                   \
            _site_type *r = 0;                                                                              \
            int number_of_elements;                                                                         \
            error(out->type != in->type, 1, "to_cpu_format_" #_name " " __FILE__,                           \
                        "Gauge field geometries do not match!\n");                                          \
            _PIECE_FOR(in->type, ixp)                                                                       \
            {                                                                                               \
                const int start = in->type->master_start[ixp];                                              \
                const int N = in->type->master_end[ixp] - in->type->master_start[ixp];                      \
                double *cin = (double*)(_4FIELD_AT(in, start, 0));                                          \
                                                                                                            \
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

    #define _FREE_GPU_CODE do {} while(0)
    #define _ALLOC_GPU_CODE(_name,_size) do {} while(0)
    #define _DECLARE_CONVERT_TO_GPU_FORMAT(_name, _field_type, _site_type, _size) 
    #define _DECLARE_CONVERT_TO_CPU_FORMAT(_name, _field_type, _site_type, _size)
    #define _DECLARE_COPY_TO(_name, _size)
    #define _DECLARE_COPY_FROM(_name, _size) 

#endif /* WITH_GPU */

/* ============================================== All cases ============================================== */

/* deallocation function declaration */
#define _DECLARE_FREE_FUNC(_name,_field_type)                                                                     \
    void free_##_name(_field_type *u)                                                                             \
    {                                                                                                       \
        if (u!=NULL) {                                                                                      \
            if (u->ptr!=NULL)                                                                               \
                afree(u->ptr);                                                                              \
            _FREE_GPU_CODE;                                                                                 \
            _FREE_MPI_CODE;                                                                                 \
            afree(u);                                                                                       \
            u = NULL;                                                                                       \
        }                                                                                                   \
    }

/* allocation function declaration */
#define _DECLARE_ALLOC_FUNC(_name, _field_type, _size)                                                              \
    _field_type *alloc_##_name(unsigned int n, geometry_descriptor *type)                                         \
    { \
        /* Allocate field struct pointer */                                                                 \
        _field_type *f;                                                                                           \
                                                                                                            \
        if (n == 0)                                                                                         \
            return NULL;                                                                                    \
        f = amalloc(n * sizeof(*f), ALIGN);                                                                 \
        error(f == NULL, 1, "alloc_" #_name " [" __FILE__ "]",                                              \
                "Could not allocate memory space for field (structure)");                                   \
                                                                                                            \
        /* Allocate CPU field data */\
        for (int i = 0; i < n; ++i)                                                                         \
            f[i].type=type;                                                                                 \
                                                                                                            \
        if(alloc_mem_t & CPU_MEM)                                                                           \
        {                                                                                                   \
            int field_size = _size * type->gsize_gauge * sizeof(*(f->ptr), ALIGN);                          \
            f->ptr = amalloc(field_size, ALIGN);                                                            \
            for(int i = 1; i < n; ++i)                                                                      \
                f[i].ptr=f[i-1].ptr + type->gsize_spinor * _size;                                           \
        }                                                                                                   \
        else                                                                                                \
        {                                                                                                   \
            for (int i = 0; i < n; ++i)                                                                     \
                f[i].ptr=NULL;                                                                              \
        }	                                                                                                \
                                                                                                            \
        /* Allocate GPU field, if compiling with GPU */                                                     \
        _ALLOC_GPU_CODE(_name,_size);                                                                       \
                                                                                                            \
        /* Allocate buffers for MPI comms, if compiling with MPI */                                         \
        _ALLOC_MPI_CODE(_name);                                                                             \
                                                                                                            \
        return f;                                                                                           \
    }

#define _DECLARE_MEMORY_FUNC(_name, _field_type, _site_type, _size)                                                             \
    _DECLARE_FREE_FUNC(_name,_field_type)                                                                         \
    _DECLARE_ALLOC_FUNC(_name,_field_type,_size)                                                                  \
    _DECLARE_COPY_TO(_name, _field_type, _size)                                                             \
    _DECLARE_COPY_FROM(_name, _field_type, _size)                                                           \
    _DECLARE_CONVERT_TO_GPU_FORMAT(_name, _field_type, _site_type, _size)                                   \
    _DECLARE_CONVERT_TO_CPU_FORMAT(_name, _field_type, _site_type, _size)

_DECLARE_MEMORY_FUNC(spinor_field_f, spinor_field, suNf_spinor, 1);
_DECLARE_MEMORY_FUNC(spinor_field_f_flt, spinor_field_flt, suNf_spinor_flt, 1);
_DECLARE_MEMORY_FUNC(sfield, scalar_field, double, 1);

#undef _DECLARE_MEMORY_FUNC
#undef _DECLARE_FREE_FUNC
#undef _DECLARE_ALLOC_FUNC
#undef _DECLARE_COPY_TO
#undef _DECLARE_COPY_FROM
#undef _DECLARE_CONVERT_TO_GPU_FORMAT
#undef _DECLARE_CONVERT_TO_CPU_FORMAT