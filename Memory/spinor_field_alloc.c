/***************************************************************************\
* Copyright (c) 2008, Claudio Pica                                          *   
* All rights reserved.                                                      * 
\***************************************************************************/

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

/* MPI allocation and deallocation code */
#ifdef WITH_MPI

#define _FREE_MPI_CODE if(u->comm_req!=NULL) afree(u->comm_req)

#define _ALLOC_MPI_CODE(_name) \
if (type->nbuffers_spinor>0) {\
f->comm_req=amalloc(n*2*type->nbuffers_spinor*sizeof(MPI_Request),ALIGN);\
error((f->comm_req)==NULL,1,"alloc_" #_name " [" __FILE__ "]",\
"Could not allocate memory space for field (MPI)");\
for (int ix=0; ix<n*2*type->nbuffers_spinor; ++ix)\
f->comm_req[ix]=MPI_REQUEST_NULL;\
for (int i=1; i<n; ++i) f[i].comm_req=f[i-1].comm_req+2*type->nbuffers_spinor;\
} else {\
for (int i=0; i<n; ++i) f[i].comm_req=NULL;\
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
err = cudaMalloc((void **) &(f->gpu_ptr), n*_size*type->gsize_spinor*sizeof(*(f->gpu_ptr))); \
error(err!=cudaSuccess,1,"alloc_" #_name " [" __FILE__ "]", \
"Could not allocate GPU memory space for field"); \
for (int i=1; i<n; ++i) f[i].gpu_ptr=f[i-1].gpu_ptr+type->gsize_spinor*_size;\
} else for (int i=0; i<n; ++i) f[i].gpu_ptr=NULL

#else /* WITH_GPU */

#define _FREE_GPU_CODE do{}while(0)
#define _ALLOC_GPU_CODE(_name,_size) do{}while(0)

#endif /* WITH_GPU */


/* deallocation function */
#define _DECLARE_FREE_FUNC(_name,_type)\
void free_##_name(_type *u){ \
if (u!=NULL) { \
if (u->ptr!=NULL) afree(u->ptr);\
_FREE_GPU_CODE;\
_FREE_MPI_CODE;\
afree(u);\
}\
}

/* allocation function */
#define _DECLARE_ALLOC_FUNC(_name,_type,_size)\
_type *alloc_##_name(unsigned int n, geometry_descriptor *type){ \
_type *f;\
\
if (n==0) return NULL;\
f=amalloc(n*sizeof(*f),ALIGN);\
error(f==NULL,1,"alloc_" #_name " [" __FILE__ "]",\
"Could not allocate memory space for field (structure)");\
for (int i=0; i<n; ++i) f[i].type=type;\
\
if(alloc_mem_t & CPU_MEM) {\
f->ptr=amalloc(n*_size*type->gsize_spinor*sizeof(*(f->ptr)),ALIGN);\
for(int i=1; i<n; ++i) f[i].ptr=f[i-1].ptr+type->gsize_spinor*_size;\
} else { for (int i=0; i<n; ++i) f[i].ptr=NULL; }	      \
\
_ALLOC_GPU_CODE(_name,_size);\
\
_ALLOC_MPI_CODE(_name);\
\
return f;\
}


/*error((f->ptr)==NULL,1,"alloc_" #_name " [" __FILE__ "]",\
"Could not allocate memory space for field (data)");\*/

/*
 _name = suffix of the allocation and deallocation functions
 _type = field type to allocate/deallocate
 _size = the number of elementary objects per lattice site
 */

#define _DECLARE_MEMORY_FUNC(_name,_type,_size) \
_DECLARE_FREE_FUNC(_name,_type);\
_DECLARE_ALLOC_FUNC(_name,_type,_size)


_DECLARE_MEMORY_FUNC(spinor_field_f, spinor_field, 1);
_DECLARE_MEMORY_FUNC(spinor_field_f_flt, spinor_field_flt, 1);

_DECLARE_MEMORY_FUNC(sfield, scalar_field, 1);

