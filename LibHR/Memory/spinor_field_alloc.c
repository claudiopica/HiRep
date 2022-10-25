/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *
* All rights reserved.                                                      *
\***************************************************************************/

/*******************************************************************************
*
* File spinor_field_alloc.c
*
* Functions for spinor field allocation and memory operations
*
*******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
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

/*
 * Here we use for all macros the parameters:
 *
 *  _name = suffix of the allocation and deallocation functions
 *  _field_type = field type to allocate/deallocate
 *  _site_type = type of elementary objects on the lattice sites
 *  _size = the number of elementary objects per lattice site
 *
 */


/* ================================================= MPI and GPU ================================================= */

#if defined(WITH_MPI) && defined(WITH_GPU)

    #define _FREE_GPU_FIELD_DATA(_name, _site_type) \
        if (f->gpu_ptr != NULL) \
        {\
            cudaError_t err;\
            _PIECE_FOR_MPI(f->type, ixp)\
            {\
                CHECK(cudaFree(f->block_handles[ixp][PID]));\
            }\
            /* TODO: Deallocate handles. */\
        }
    
    /* This only works for n=1 */
    #define _ALLOC_GPU_FIELD_DATA(_name, _site_type, _size) \
        cudaError_t err;\
        int block_size = 0;\
        _PIECE_FOR_MPI(f->type, ixp) \
        {\
            f->block_handles[ixp] = (_site_type**)malloc(MPI_WORLD_SIZE*sizeof(_site_type*));\
            block_size = f->type->master_end[ixp] - f->type->master_start[ixp] + 1;\
            f->block_handles[ixp][PID] = _GPU_FIELD_BLK(f, ixp);\
            int mem_size = (_size)*block_size*sizeof(*(f->gpu_ptr));\
            CHECK(cudaMalloc((void **)&(f->block_handles[ixp][PID]), mem_size));\
        }


    /*#define _FREE_MPI_FIELD_DATA */
        /*if (f->comm_req != NULL) \
            cudaFree(f->comm_req);*/

    /*#define _ALLOC_MPI_FIELD_DATA(_name) */
        /*if (type->nbuffers_gauge > 0) \
        {\
            cudaError_t err;\
            err = cudaMalloc((void **)f->comm_req, 2 * type->nbuffers_gauge * sizeof(MPI_Request));\
            error(err != cudaSuccess, 1, "alloc_" #_name " [" __FILE__ "]",\
                            "Could not allocate buffers for multi-GPU communication.\n");\
            for(int ix = 0; ix < 2 * type->nbuffers_gauge; ++ix) \
                f->comm_req[ix] = MPI_REQUEST_NULL;\
        }\
        else\
        {\
            f->comm_req = NULL;\
        }*/

    #define _DECLARE_COPY_TO(_name, _field_type, _site_type, _size) \
        void copy_to_gpu_##_name(_field_type *f)\
        {\
            /*_field_type *tmp = alloc_##_name(1, f->type);\
            to_gpu_format_##_name(tmp, f);\
            cudaError_t err;\
            int block_size = 0;\
            int active_device = 0;\
            _site_type *block_start_in, *block_start_tmp;\
            _PIECE_FOR(f->type, ixp) \
            {\
                block_size = f->type->master_end[ixp] - f->type->master_start[ixp] + 1;\
                block_start_tmp = _FIELD_BLK(tmp, tmp->type->master_start[ixp]);\
                block_start_in = _GPU_FIELD_BLK(f, f->type->master_start[ixp]);\
\
                active_device = ixp / ngpus;\
                _CHANGE_DEVICE(_name);\
                cudaMemcpy(block_start_in, block_start_tmp, block_size, cudaMemcpyHostToDevice);\
            }\
            free_##_name(tmp);*/\
        }
    
    #define _DECLARE_COPY_FROM(_name, _field_type, _site_type, _size) \
        void copy_from_gpu_##_name(_field_type *f) \
        {\
            /*_field_type *tmp = alloc_##_name(1, f->type);\
            cudaError_t err;\
            _QUERY_NGPUS(_name);\
            int block_size = 0;\
            int active_device = 0;\
            _site_type *block_start_in, *block_start_tmp;\
            _PIECE_FOR(f->type, ixp) \
            {\
                block_size = f->type->master_end[ixp] - f->type->master_start[ixp] + 1;\
                block_start_tmp = _FIELD_BLK(tmp, tmp->type->master_start[ixp]);\
                block_start_in = _FIELD_BLK(f, f->type->master_start[ixp]);\
\
                active_device = ixp / ngpus; \
                _CHANGE_DEVICE(_name); \
                cudaMemcpy(block_start_tmp, block_start_in, block_size, cudaMemcpyDeviceToHost);\
            }\
            to_cpu_format_##_name(f, tmp);\
            free_##_name(tmp);*/\
        }

#endif

/* ================================================= MPI and CPU ========================================= */

#if defined(WITH_MPI) //&& !defined(WITH_GPU)

        #define _FREE_MPI_FIELD_DATA                                                                        \
            if(f->comm_req != NULL)                                                                         \
                afree(f->comm_req)

        #define _ALLOC_MPI_FIELD_DATA(_name)                                                                \
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
            }

#endif 

/* ================================================= Single GPU ========================================== */

#if defined(WITH_GPU) && !defined(WITH_MPI)

    /* Free device memory */
    /* Note: to be used inside function declaration */
    #define _FREE_GPU_FIELD_DATA(_name, _site_type)                                                         \
        if(f->gpu_ptr != NULL)                                                                              \
            cudaFree(f->gpu_ptr);                                                                           \

    #define _ALLOC_GPU_FIELD_DATA(_name, _site_type, _size)                                                 \
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

    /* Declare function to copy field from host to device */
    #define _DECLARE_COPY_TO(_name, _field_type, _site_type, _size)                                         \
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
    #define _DECLARE_COPY_FROM(_name, _field_type, _site_type, _size)                                       \
        void copy_from_gpu_##_name(_field_type *f)                                                          \
        {                                                                                                   \
            _field_type *tmp = alloc_##_name(1, f->type);                                                   \
            int field_size = _size * f->type->gsize_gauge * sizeof(*(f->gpu_ptr));                          \
            cudaMemcpy(tmp->ptr, f->gpu_ptr, field_size, cudaMemcpyDeviceToHost);                           \
            to_cpu_format_##_name(f, tmp);                                                                  \
            free_##_name(tmp);                                                                              \
        }
#endif

/* ================================================= Single and Multi-GPU ================================ */

/* These macros work with or without MPI! */
#ifdef WITH_GPU

    /* Declare function to convert GPU to CPU format */
    /* This is necessary, because GPU performance is very sensitive to */
    /* memory access patterns, see documentation Doc/gpu_geometry.tex */
    #define _DECLARE_CONVERT_TO_GPU_FORMAT(_name, _field_type, _site_type, _size)                           \
        void to_gpu_format_##_name(_field_type *out, _field_type *in)                                       \
        {                                                                                                   \
            _site_type *source, *target;                                                                    \
            _CHECK_GEOMETRY_MATCHING(out, in);                                                              \
                                                                                                            \
            int stride = 0;                                                                                 \
            _PIECE_FOR(in->type, ixp)                                                                       \
            {                                                                                               \
                stride = in->type->master_end[ixp] - in->type->master_start[ixp] + 1;                       \
                target  = _FIELD_BLK(out, ixp);                                                             \
                _SITE_FOR(in->type, ixp, ix)                                                                \
                {                                                                                           \
                    source = _FIELD_AT(in, ix);                                                             \
                    int ix_loc = _GPU_IDX_TO_LOCAL(in, ix, ixp);                                            \
                    write_gpu_##_site_type(stride, (*source), target, ix_loc, 0);                           \
                                                                                                            \
                }                                                                                           \
            }                                                                                               \
        }

    #define _DECLARE_CONVERT_TO_CPU_FORMAT(_name, _field_type, _site_type, _size)                           \
        void to_cpu_format_##_name(_field_type *out, _field_type *in)                                       \
        {                                                                                                   \
            _site_type *target, *source;                                                                    \
            _CHECK_GEOMETRY_MATCHING(out, in);                                                              \
                                                                                                            \
            int stride = 0;                                                                                 \
            int ix_loc = 0;\
            _PIECE_FOR(in->type, ixp)                                                                       \
            {                                                                                               \
                stride = in->type->master_end[ixp] - in->type->master_start[ixp] + 1;\
                source = _FIELD_BLK(in, ixp);                                                               \
                _SITE_FOR(in->type, ixp, ix)                                                                \
                {                                                                                           \
                    ix_loc = _GPU_IDX_TO_LOCAL(in, ix, ixp);\
                    target = _FIELD_AT(out, ix);                                                            \
                    read_gpu_##_site_type(stride, (*target), source, ix_loc, 0);                            \
                }                                                                                           \
            }                                                                                               \
        }

#endif

/* ================================================= Empty defs ========================================= */

#ifndef WITH_GPU

    #define _FREE_GPU_FIELD_DATA(_name, _site_type) do {} while(0)
    #define _ALLOC_GPU_FIELD_DATA(_name, _site_type, _size) do {} while(0)
    #define _DECLARE_COPY_TO(_name, _field_type, _site_type, _size)
    #define _DECLARE_COPY_FROM(_name, _field_type, _site_type, _size) 
    #define _DECLARE_CONVERT_TO_GPU_FORMAT(_name, _field_type, _site_type, _size) 
    #define _DECLARE_CONVERT_TO_CPU_FORMAT(_name, _field_type, _site_type, _size)
    
#endif

#ifndef WITH_MPI

    #define _FREE_MPI_FIELD_DATA do {} while(0)
    #define _ALLOC_MPI_FIELD_DATA(_name) do {} while(0)

#endif

/* ============================================== All cases ============================================== */

/* deallocation function declaration */
#define _DECLARE_FREE_FUNC(_name, _field_type, _site_type)                                                  \
    void free_##_name(_field_type *f)                                                                       \
    {                                                                                                       \
        if (f!=NULL) {                                                                                      \
            if (f->ptr!=NULL)                                                                               \
                afree(f->ptr);                                                                              \
            _FREE_GPU_FIELD_DATA(_name, _site_type);                                                        \
            _FREE_MPI_FIELD_DATA;                                                                           \
            afree(f);                                                                                       \
            f = NULL;                                                                                       \
        }                                                                                                   \
    }

/* allocation function declaration */
#define _DECLARE_ALLOC_FUNC(_name, _field_type, _site_type, _size)                                          \
    _field_type *alloc_##_name(unsigned int n, geometry_descriptor *type)                                   \
    { \
        /* Allocate field struct pointer */                                                                 \
        _field_type *f;                                                                                     \
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
            int field_size = _size * type->gsize_gauge * sizeof(*(f->ptr));                                 \
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
        _ALLOC_GPU_FIELD_DATA(_name, _site_type, _size);                                                    \
                                                                                                            \
        /* Allocate buffers for MPI comms, if compiling with MPI */                                         \
        _ALLOC_MPI_FIELD_DATA(_name);                                                                       \
                                                                                                            \
        return f;                                                                                           \
    }

#define _DECLARE_MEMORY_FUNC(_name, _field_type, _site_type, _size)                                         \
    _DECLARE_FREE_FUNC(_name, _field_type, _site_type)                                                      \
    _DECLARE_ALLOC_FUNC(_name,_field_type, _site_type, _size)                                               \
    _DECLARE_COPY_TO(_name, _field_type, _site_type, _size)                                                 \
    _DECLARE_COPY_FROM(_name, _field_type, _site_type, _size)                                               \
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
#undef _ALLOC_GPU_FIELD_DATA
#undef _ALLOC_MPI_FIELD_DATA
#undef _FREE_GPU_FIELD_DATA
#undef _FREE_MPI_FIELD_DATA