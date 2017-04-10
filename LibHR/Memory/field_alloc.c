/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *   
* All rights reserved.                                                      * 
\***************************************************************************/

/*******************************************************************************
*
* File field_alloc.c
*
* Functions for fields allocation
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


/* MPI allocation and deallocation code */
#ifdef WITH_MPI

#define _FREE_MPI_CODE if(u->comm_req!=NULL) afree(u->comm_req)

#define _ALLOC_MPI_CODE(_name) \
if (type->nbuffers_gauge>0) {\
f->comm_req=amalloc(2*type->nbuffers_gauge*sizeof(MPI_Request),ALIGN);\
error((f->comm_req)==NULL,1,"alloc_" #_name " [" __FILE__ "]",\
"Could not allocate memory space for field (MPI)");\
for (int ix=0; ix<2*type->nbuffers_gauge; ++ix)\
f->comm_req[ix]=MPI_REQUEST_NULL;\
} else {\
f->comm_req=NULL;\
}do{}while(0)

#else /* WITH_MPI */

#define _FREE_MPI_CODE do{}while(0)
#define _ALLOC_MPI_CODE(_name) do{}while(0)

#endif /* WITH_MPI */


/* GPU allocation and deallocation code */
#ifdef WITH_GPU

#define _FREE_GPU_CODE if(u->gpu_ptr!=NULL) cudaFree(u->gpu_ptr)


#define _ALLOC_GPU_CODE(_name,_size)\
if(alloc_mem_t & GPU_MEM) {\
cudaError_t err;\
err = cudaMalloc((void **) &f->gpu_ptr, _size*type->gsize_gauge*sizeof(*(f->gpu_ptr))); \
error(err!=cudaSuccess,1,"alloc_" #_name " [" __FILE__ "]", \
"Could not allocate GPU memory space for field"); \
} else f->gpu_ptr=NULL

#else /* WITH_GPU */

#define _FREE_GPU_CODE do{}while(0)
#define _ALLOC_GPU_CODE(_name,_size) do{}while(0)

#endif /* WITH_GPU */


/* deallocation function */
#define _DECLARE_FREE_FUNC(_name,_type)\
void free_##_name(_type *u){ \
if (u!=NULL) {\
if (u->ptr!=NULL) afree(u->ptr);\
_FREE_GPU_CODE;\
_FREE_MPI_CODE;\
afree(u);\
}\
}

/* allocation function */
#define _DECLARE_ALLOC_FUNC(_name,_type,_size)\
_type *alloc_##_name(geometry_descriptor *type){ \
_type *f;\
\
f=amalloc(sizeof(*f),ALIGN);\
error(f==NULL,1,"alloc_" #_name " [" __FILE__ "]",\
"Could not allocate memory space for field (structure)");\
f->type=type;\
\
if(alloc_mem_t & CPU_MEM) {\
f->ptr=amalloc(_size*type->gsize_gauge*sizeof(*(f->ptr)),ALIGN);\
error((f->ptr)==NULL,1,"alloc_" #_name " [" __FILE__ "]",\
"Could not allocate memory space for field (data)");\
} else { f->ptr=NULL; }\
\
_ALLOC_GPU_CODE(_name,_size);\
\
_ALLOC_MPI_CODE(_name);\
\
return f;\
}



/*
 _name = suffix of the allocation and deallocation functions
 _type = field type to allocate/deallocate
 _size = the number of elementary objects per lattice site
 */

#define _DECLARE_MEMORY_FUNC(_name,_type,_size) \
_DECLARE_FREE_FUNC(_name,_type);\
_DECLARE_ALLOC_FUNC(_name,_type,_size)

_DECLARE_MEMORY_FUNC(gfield, suNg_field, 4);
_DECLARE_MEMORY_FUNC(gfield_flt, suNg_field_flt, 4);

_DECLARE_MEMORY_FUNC(gfield_f, suNf_field, 4);
_DECLARE_MEMORY_FUNC(gfield_f_flt, suNf_field_flt, 4);
_DECLARE_MEMORY_FUNC(scalar_field, suNg_scalar_field, 1);

_DECLARE_MEMORY_FUNC(avfield, suNg_av_field, 4);
_DECLARE_MEMORY_FUNC(gtransf, suNg_field, 1);

_DECLARE_MEMORY_FUNC(clover_ldl, ldl_field, 1);
_DECLARE_MEMORY_FUNC(clover_term, suNfc_field, 4);
_DECLARE_MEMORY_FUNC(clover_force, suNf_field, 6);
