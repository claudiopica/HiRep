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

/*
 * Here we use for all macros the parameters:
 *
 *  _name = suffix of the allocation and deallocation functions
 *  _field_type = field type to allocate/deallocate
 *  _site_type = type of elementary objects on the lattice sites
 *  _size = the number of elementary objects per lattice site
 *
 */

/* ================================================= MPI and GPU ========================================= */

#if defined(WITH_MPI) && defined(WITH_GPU)

    /* From the book */
    #define CHECK(call)										                                                \
    {												                                                        \
        const cudaError_t error = call;								                                        \
        if (error != cudaSuccess) 									                                        \
        {												                                                    \
            fprintf(stderr, "Error: %s:%d, ", __FILE__, __LINE__);					                        \
	        fprintf(stderr, "code: %d, reason: %s\n", error, cudaGetErrorString(error));		            \
	        exit(1);										                                                \
        }												                                                    \
    }	

    #define _FREE_GPU_FIELD_DATA(_name, _site_type)                                                         \
        if (f->gpu_ptr != NULL)                                                                             \
        {                                                                                                   \
            cudaError_t err;                                                                                \
            int device_id;                                                                                  \
            cudaGetDevice(&device_id);                                                                      \
            for (int ixp = device_id*2; ixp<2*(device_id+1); ++ixp)                                         \
            {                                                                                               \
                printf("Deallocating piece %d on device %d\n", ixp, device_id);                             \
                CHECK(cudaFree(f->block_handles[ixp][device_id]));                                          \
            }                                                                                               \
            /*TODO: Deallocate handles.*/\
        }

    /* Allocate device memory */
    /* Note: to be used inside function declaration */
    #define _ALLOC_GPU_FIELD_DATA(_name, _site_type, _size)                                                 \
        cudaError_t err;                                                                                    \
        int block_size = 0;                                                                                 \
        int device_id;                                                                                      \
        cudaGetDevice(&device_id);                                                                          \
        /* TODO: The map 2*device_id -> ixp might not be ideal. */                                          \
        for (int ixp = device_id*2; ixp<(2*device_id)+2; ++ixp)                                             \
        {                                                                                                   \
            f->block_handles[ixp] = (_site_type**)malloc(MPI_WORLD_SIZE*sizeof(_site_type*));                \
            block_size = f->type->master_end[ixp] - f->type->master_start[ixp] + 1;                         \
            f->block_handles[ixp][device_id] = _GPU_4FIELD_BLK(f, ixp);                                     \
            printf("Allocating piece %d on device %d\n", ixp, device_id);                                   \
            CHECK(cudaMalloc((void **)&(f->block_handles[ixp][device_id]), (_size)*block_size*sizeof(*(f->gpu_ptr))));\
        }

    //#define _FREE_MPI_FIELD_DATA                                                                            
        /*if (f->comm_req != NULL)                                                                          \
            cudaFree(f->comm_req);*/

    //#define _ALLOC_MPI_FIELD_DATA(_name)                                                                    
       /* if (type->nbuffers_gauge > 0)                                                                     \
        {                                                                                                   \
            cudaError_t err;                                                                                \
            *//* How do we allocate on the right device? *//*                                               \
            err = cudaMalloc((void **)f->comm_req, 2 * type->nbuffers_gauge * sizeof(MPI_Request));         \
            error(err != cudaSuccess, 1, "alloc_" #_name " [" __FILE__ "]",                                 \
                            "Could not allocate buffers for multi-GPU communication.\n");                   \
            for (int ix = 0; ix < 2 * type->nbuffers_gauge; ++ix)                                           \
                f->comm_req[ix] = MPI_REQUEST_NULL;                                                         \
        }                                                                                                   \
        else                                                                                                \
        {                                                                                                   \
            f->comm_req = NULL;                                                                             \
        }*/

        /* Declare function to copy field from host to device */
    #define _DECLARE_COPY_TO(_name, _field_type, _site_type, _size)                                         \
        void copy_to_gpu_##_name(_field_type *f)                                                            \
        {                                                                                                   \
            _field_type *tmp = alloc_##_name(f->type);                                                      \
            to_gpu_format_##_name(tmp, f);                                                                  \
            cudaError_t err;                                                                                \
            int block_size = 0;                                                                             \
            int device_id;                                                                                  \
            cudaGetDevice(&device_id);                                                                      \
            _site_type *block_start_in, *block_start_tmp;                                                   \
            for (int ixp = device_id*2; ixp<(2*device_id)+2; ++ixp)                                         \
            {                                                                                               \
                block_size = f->type->master_end[ixp] - f->type->master_start[ixp] + 1;                     \
                block_start_tmp = _4FIELD_BLK(tmp, tmp->type->master_start[ixp]);                           \
                block_start_in = _GPU_4FIELD_BLK(f, f->type->master_start[ixp]);                            \
                                                                                                            \
                cudaMemcpy(block_start_in, block_start_tmp, block_size, cudaMemcpyHostToDevice);            \
            }                                                                                               \
            cudaDeviceSynchronize();                                                                        \
            free_##_name(tmp);                                                                              \
        }

    /* Declare function to copy field from device to host */
    #define _DECLARE_COPY_FROM(_name, _field_type, _site_type, _size)                                       \
        void copy_from_gpu_##_name(_field_type *f)                                                          \
        {                                                                                                   \
            _field_type *tmp = alloc_##_name(f->type);                                                      \
            cudaError_t err;                                                                                \
            int block_size = 0;                                                                             \
            int device_id;                                                                                  \
            cudaGetDevice(&device_id);                                                                      \
            _site_type *block_start_in, *block_start_tmp;                                                   \
            for (int ixp = device_id*2; ixp<(2*device_id)+2; ++ixp)                                         \
            {                                                                                               \
                block_size = f->type->master_end[ixp] - f->type->master_start[ixp] + 1;                     \
                block_start_tmp = _4FIELD_BLK(tmp, tmp->type->master_start[ixp]);                           \
                block_start_in = _GPU_4FIELD_BLK(f, f->type->master_start[ixp]);                            \
                                                                                                            \
                cudaMemcpy(block_start_tmp, block_start_in, block_size, cudaMemcpyDeviceToHost);            \
            }                                                                                               \
            cudaDeviceSynchronize();                                                                        \
            to_cpu_format_##_name(f, tmp);                                                                  \
            free_##_name(tmp);                                                                              \
        }

#endif

/* ================================================= MPI and CPU ========================================= */

#if defined(WITH_MPI) //&& !defined(WITH_GPU)

        #define _FREE_MPI_FIELD_DATA                                                                        \
            if (f->comm_req != NULL)                                                                        \
                afree(f->comm_req)                                                                          \

        #define _ALLOC_MPI_FIELD_DATA(_name)                                                                \
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
  
#endif

/* ================================================= Single GPU ========================================== */

#if defined(WITH_GPU) && !defined(WITH_MPI)

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
                err = cudaMalloc((void **)&f->gpu_ptr, field_size);                                         \
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

#endif

/* ================================================= Single and Multi-GPU ================================ */

/* These macros work with or without MPI! */
#ifdef WITH_GPU
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
                target = _4FIELD_BLK(out, ixp);                                                             \
                _SITE_FOR(in->type, ixp, ix)                                                                \
                {                                                                                           \
                    ix_loc = _GPU_IDX_TO_LOCAL(in, ix, ixp);                                                \
                    for (int comp = 0; comp < _size; ++comp)                                                \
                    {                                                                                       \
                        source = _4FIELD_AT(in, ix, comp);                                                  \
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
                source = _4FIELD_BLK(in, ixp);                                                              \
                _SITE_FOR(in->type, ixp, ix)                                                                \
                {                                                                                           \
                    ix_loc = _GPU_IDX_TO_LOCAL(in, ix, ixp);                                                \
                    for (int comp = 0; comp < _size; ++comp)                                                \
                    {                                                                                       \
                        target = _4FIELD_AT(out, ix, comp);                                                 \
                        read_gpu_##_site_type(stride, (*target), source, ix_loc, comp);                     \
                    }                                                                                       \
                }                                                                                           \
            }                                                                                               \
        }

#endif

/* ================================================= Empty defs ========================================= */

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
    #define _ALLOC_MPI_FIELD_DATA(_name) do {} while (0)

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
        printf("Executing free on dev %d\n", PID);\
        if (f != NULL)                                                                                      \
        {                                                                                                   \
            if (f->ptr != NULL)                                                                             \
                afree(f->ptr);                                                                              \
            _FREE_GPU_FIELD_DATA(_name, _site_type);/* Every MPI process needs to free -> check, whether its null differently...*/\
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
        _ALLOC_MPI_FIELD_DATA(_name);                                                                       \
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